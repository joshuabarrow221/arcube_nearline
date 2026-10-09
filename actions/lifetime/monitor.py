"""Configured watcher integration for pooled detector fits and slow controls.

One serialized state namespace owns ingestion, fit membership and publication.
This favors reproducibility over parallel fits to the same growing window.
Immutable array objects and cached fits retain audit evidence when a file is
replaced. Scientific fit failures are data; extraction/IO failures also make
this command exit unsuccessfully after publishing other available methods.
"""
from datetime import datetime, timedelta, timezone
from fnmatch import fnmatch
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from lifetime_io import output_lock, public_permissions
from monitor_config import load_config, paired_input, file_version
from packet_pooling import pool_samples, digest_json
from purity_sources import atomic_json, aggregate_measurements, load_measurements, prm_observations
from track_shards import write_shard, pool_tracks


def cohort_for(relative, config):
    """Conservative directory cohorts, with explicit operator-reviewed overrides.

    No regex is fitted to observed lifetimes. More than one matching rule is
    an error. Without a rule every Step/trial directory remains distinct; an
    operator can merge known compatible periods only after independent review.
    """
    matches = [r for r in config['cohort_rules'] if fnmatch(relative, r['pattern'])]
    if len(matches) > 1: raise ValueError('ambiguous cohort rules for '+relative)
    cohort = matches[0]['cohort'] if matches else str(Path(relative).parent)
    diagnostic = matches[0].get('diagnostic', False) if matches else any(
        fnmatch(relative, pattern) for pattern in config['diagnostic_patterns'])
    return cohort, diagnostic


def import_packet_sample(descriptor, directory, identity, metadata):
    """Promote one extraction using full acquisition identity, not its basename.

    The historical extractor writes into an input-specific attempt directory,
    so identical basenames in distinct runs cannot collide. The shared pool
    descriptor then uses the FLOW-relative identity and an immutable absolute
    object path. Existing packet gates and object hashes remain unchanged.
    """
    row = json.loads(Path(descriptor).read_text())
    objects = Path(descriptor).parent/row['objects_file']
    row.update(metadata, sample_id=identity, objects_file=str(objects.resolve()))
    atomic_json(Path(directory)/(identity+'.sample.json'), row)
    return row


def ingest(path, config, state):
    """Process a completed FLOW/native generation exactly once, including late packets."""
    pair = paired_input(path, config)
    if pair is None:
        if (Path(path).parent.resolve() != Path(config['snapshot_root']) or
                not Path(path).match(config['snapshot_pattern']) or not Path(path).is_file()):
            raise ValueError('input is outside the configured acquisition/snapshot trees')
        return []  # Snapshot import happens coherently across the whole archive.
    flow, native, relative = pair
    if not flow.is_file(): raise FileNotFoundError('matching FLOW is not available: '+str(flow))
    identity = hashlib.sha256(relative.encode()).hexdigest()[:24]
    folder = state/'inputs'/identity; folder.mkdir(parents=True, exist_ok=True)
    cohort, diagnostic = cohort_for(relative, config)
    from nearline_util import date_from_filename
    timestamp = date_from_filename(str(flow)).isoformat()
    metadata = dict(timestamp=timestamp, cohort=cohort, monitoring_excluded=diagnostic,
                    flow_file=str(flow), relative_flow=relative)
    # Include implementation changes in receipt reuse, but do not include the
    # display style: restyling a graph should never re-extract millions of packets.
    algorithms = ('track_shards.py', 'packet_lifetime.py', 'packet_pedestal.py', 'calibration_pedestal.py')
    code = hashlib.sha256(b''.join(Path(__file__).with_name(p).read_bytes() for p in algorithms)).hexdigest()
    versions = dict(flow=file_version(flow), packet=file_version(native) if config['packet_enabled'] else None,
                    config=config['_processing_revision'], code=code)
    receipt_path = folder/'receipt.json'
    previous = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    if previous.get('versions') == versions and not previous.get('errors'):
        return []
    errors = []
    try:
        write_shard(flow, state/'tracks', identity, metadata)
    except Exception as exc:
        errors.append('track extraction: '+str(exc))
        (state/'tracks'/(identity+'.sample.json')).unlink(missing_ok=True)
    packet_status = 'disabled' if not config['packet_enabled'] else 'waiting for matching native packet'
    # A stale descriptor must not survive replacement of its source generation.
    packet_descriptor = state/'packets'/(identity+'.sample.json')
    packet_descriptor.unlink(missing_ok=True)
    if config['packet_enabled'] and native.is_file():
        try:
            with __import__('h5py').File(flow, 'r') as h:
                pedestal = h['charge/calib_prompt_hits'].attrs.get('pedestal_file', '')
                if isinstance(pedestal, bytes): pedestal = pedestal.decode()
            # Prefer the exact static calibration named in FLOW. An optional
            # calibration_root relocates only its basename, for portable replay.
            if config.get('calibration_root'):
                pedestal = str(Path(os.path.expandvars(config['calibration_root']))/Path(pedestal).name)
            if not pedestal or not Path(pedestal).is_file():
                raise FileNotFoundError('the exact FLOW-named pedestal calibration is unavailable')
            attempt = Path(tempfile.mkdtemp(prefix='packet-', dir=folder))
            extraction = attempt/'samples'; extraction.mkdir()
            command = [sys.executable, '-u', str(Path(__file__).with_name('packet_lifetime.py')),
                str(native), '--flow-file', str(flow), '--ped-source', 'calibration', '--ped', pedestal,
                '--outdir', str(attempt), '--detector-wide', '--charge-mode', 'adc-minus-ped',
                '--successive-min-ticks', '28', '--successive-max-ticks', '28', '--ext-io', str(config['packet_ext_io']),
                '--sample-directory', str(extraction), '--sample-timestamp', timestamp,
                '--cohort', cohort, '--extract-only']
            with (attempt/'analysis.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            files = list(extraction.glob('*.sample.json'))
            if len(files) != 1: raise ValueError('packet extraction did not publish exactly one sample')
            import_packet_sample(files[0], state/'packets', identity, metadata)
            packet_status = 'extracted'
        except Exception as exc:
            errors.append('packet extraction: '+str(exc)); packet_status = 'failed'
    # Reject a source that changed while being read. No partial generation is
    # allowed to enter an otherwise valid growing pool.
    if versions['flow'] != file_version(flow) or (config['packet_enabled'] and versions['packet'] != file_version(native)):
        errors.append('input changed during extraction')
        (state/'tracks'/(identity+'.sample.json')).unlink(missing_ok=True)
        packet_descriptor.unlink(missing_ok=True)
    receipt = dict(input_file=str(flow), relative_flow=relative, versions=versions,
                   packet_status=packet_status, errors=errors, updated_at=datetime.now(timezone.utc).isoformat())
    atomic_json(folder/(digest_json(receipt)+'.receipt.json'), receipt)
    atomic_json(receipt_path, receipt)
    return errors


def publish(payload, config):
    """Render a complete generation before atomically switching stable web URLs.

    Each generation is retained for rollback. Stable PNG, HTML and diagnostic
    aliases point through one '.lifetime-current' directory symlink. Readers
    never see a partly written image or a half-written interactive document.
    Existing regular files are backed up on the first local/site installation.
    """
    from monitor_plot import render
    root = Path(config['publish_root']); root.mkdir(parents=True, exist_ok=True)
    generations = root/'.lifetime-generations'; generations.mkdir(exist_ok=True)
    generation = Path(tempfile.mkdtemp(prefix=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S-'), dir=generations))
    atomic_json(generation/'elifetime_status.json', payload)
    render(payload, generation, config)
    names = ['elifetime_time_series.png', 'elifetime_time_series.png.html',
             'elifetime_diagnostics.html', 'elifetime_status.json']
    for p in [root, generations, generation]:
        try: p.chmod(0o755)
        except PermissionError:
            if os.environ.get('ARCUBE_NEARLINE_LOCAL_OUTPUT') != '1': raise
    for name in names: public_permissions(generation/name)
    link = root/('.current-'+generation.name)
    link.symlink_to(generation.relative_to(root), target_is_directory=True)
    os.replace(link, root/'.lifetime-current')
    for name in names:
        alias = root/name; expected = '.lifetime-current/'+name
        if alias.is_symlink() and os.readlink(alias) == expected: continue
        if alias.exists():
            backup = root/'.previous'/generation.name; backup.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(alias, backup/name)
        temp = root/('.alias-'+name)
        temp.unlink(missing_ok=True); temp.symlink_to(expected); os.replace(temp, alias)
    return str(generation)


def run_monitor(config_path, input_path=None):
    """Ingest, recompute active windows, and publish even with no FLOW available."""
    config = load_config(config_path); state = Path(config['state_root'])
    for name in ('tracks', 'packets', 'inputs', 'track_fits'): (state/name).mkdir(parents=True, exist_ok=True)
    # NERSC uses the shared-filesystem flufl lease. It outlives the documented
    # single-worker job limit (6 hours); local WSL uses explicit kernel flock.
    # Sites should keep the worker wall limit below this lease, not silently
    # switch to a host-local lock on shared storage.
    with output_lock(str(state/'.monitor.lock'), default_timeout=timedelta(hours=6), lifetime=timedelta(hours=24)):
        errors = ingest(input_path, config, state) if input_path else []
        tracks = pool_tracks(state/'tracks', state/'track_fits', 6, config['timezone'])
        packets = (pool_samples(state/'packets', 6, config['timezone'], bootstrap=config['bootstrap'])
                   if config['packet_enabled'] else [])
        for r in packets:
            r['monitoring_excluded'] = any(json.loads((state/'packets'/(m['sample_id']+'.sample.json')).read_text()).get('monitoring_excluded', False) for m in r['members'])
        snapshots = sorted(Path(config['snapshot_root']).glob(config['snapshot_pattern']))
        frame = load_measurements(snapshots)
        quality = json.loads(Path(config['gas_quality_config']).read_text()) if config.get('gas_quality_config') else None
        conversion = json.loads(Path(config['gas_conversion_config']).read_text()) if config.get('gas_conversion_config') else None
        controls = aggregate_measurements([], config['timezone'], quality, measurements=frame, conversion_config=conversion)
        receipts = [json.loads(p.read_text()) for p in (state/'inputs').glob('*/receipt.json')]
        payload = dict(schema_version=6, generated_at=datetime.now(timezone.utc).isoformat(),
            timezone=config['timezone'], track_window_hours=6, packet_window_hours=6, prm_window='local calendar day',
            lifetimes=tracks+packets+controls, raw_prm_measurements=prm_observations(frame),
            inputs=receipts, slow_control_snapshots=[str(p) for p in snapshots],
            input_errors=[dict(input_file=r['input_file'], errors=r['errors']) for r in receipts if r['errors']],
            uncertainty='Conditional fit/SEM errors; shared calibration, selection and event correlations excluded.',
            config_revision=config['_revision'])
        atomic_json(state/'history.json', payload)
        generation = publish(payload, config)
        print('Published lifetime generation:', generation, flush=True)
        if errors: raise RuntimeError('; '.join(errors))
    return payload
