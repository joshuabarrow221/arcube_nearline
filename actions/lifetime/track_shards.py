"""Persistent primary-FLOW segments and fixed six-hour joint track fits.

The stored selection is read once per file generation. Explicit H5Flow links
partition a legacy mixed selection into event-timing categories. Counts and
source-row lineage stay inspectable; neither per-file fitted taus nor published
JSON enter these pools. A failed fit remains a status record with null lifetime.
"""
from collections import defaultdict
from datetime import timedelta
import hashlib
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np
from lifetime import SPLIT_SAMPLES, fit_samples
from track_pooling import reconstruction_signature, json_value
from packet_pooling import digest_json, window_start
from purity_sources import atomic_json


def selection_metadata(flow):
    """Fingerprint recorded selector settings, excluding per-file results.

    Geometry and calibration metadata alone cannot distinguish changed track
    cuts. Keep both group and dataset attributes, including unfamiliar settings
    conservatively, but omit known counts and optional lifetime-fit outputs.
    Older files may record no settings; an empty dictionary documents that
    missing evidence rather than claiming the selector has been verified.
    """
    result = {}
    for sample, path in (*SPLIT_SAMPLES, ('all_mip', 'analysis/rock_muon_segments/data')):
        if path not in flow:
            continue
        result[sample] = {}
        for location, node in [('group', flow[str(Path(path).parent)]), ('dataset', flow[path])]:
            result[sample][location] = json_value({
                key: value for key, value in node.attrs.items()
                if key not in ('selected_track_count', 'selected_segment_count')
                and not key.startswith('electron_lifetime')
            })
    return result


def stored_segments(flow):
    """Read split output, or derive a disjoint timing partition of legacy rows.

    All segments from a selected track inherit their parent event's trigger
    category. ``cosmic`` has zero recorded external triggers; ``beam`` has one
    or more. The names do not identify physical cosmic/beam origin. Later,
    fit_samples concatenates these disjoint arrays into ``all_mip`` (All tracks),
    so the cosmic-enriched and All tracks curves intentionally share segments.
    This partition does not apply different geometry or charge cuts to either
    category and cannot restore candidates omitted by the upstream selection.
    """
    split = {s: flow[p][:] for s, p in SPLIT_SAMPLES if p in flow}
    if split:
        if len(split) != 2: raise ValueError('incomplete pair of FLOW timing selections')
        metadata = selection_metadata(flow)
        comparable = {}
        for sample, expected_timing in [('beam', 'beam'), ('cosmic', 'off_beam')]:
            comparable[sample] = {}
            for location, settings in metadata[sample].items():
                # The two timing cuts must be disjoint when explicitly stored.
                # A mislabeled "all" sample would otherwise double-count the
                # zero-trigger subset when fit_samples forms the union.
                timing = settings.get('event_timing')
                if timing is not None and timing != expected_timing:
                    raise ValueError('unexpected event_timing for '+sample+' selection')
                comparable[sample][location] = {
                    key: value for key, value in settings.items() if key != 'event_timing'
                }
        if comparable['beam'] != comparable['cosmic']:
            raise ValueError('split FLOW selections have different recorded settings')
        return split, 'stored split FLOW selection', {}
    path = 'analysis/rock_muon_segments/data'
    if path not in flow: raise ValueError('FLOW has no stored muon selection; run the selection workflow first')
    segments = flow[path][:]
    track_path = 'analysis/rock_muon_tracks/data'
    ref_path = 'analysis/rock_muon_tracks/ref/analysis/rock_muon_segments/ref'
    if track_path not in flow or ref_path not in flow:
        # Older FLOW can still support the combined estimate. Do not invent a
        # cosmic/external partition when its association evidence is missing.
        return {'all_mip': segments}, 'legacy FLOW selection without timing references', {}
    tracks, refs = flow[track_path][:], flow[ref_path][:]
    attrs = json_value(dict(flow[ref_path].attrs))
    if attrs['dset0'].lstrip('/') != track_path or attrs['dset1'].lstrip('/') != path:
        raise ValueError('unexpected track-to-segment reference orientation')
    if refs.shape != (len(segments), 2) or not np.array_equal(np.sort(refs[:,1]), np.arange(len(segments))):
        raise ValueError('each stored segment must have exactly one track reference')
    if np.any(refs[:,0] < 0) or np.any(refs[:,0] >= len(tracks)): raise ValueError('invalid track reference')
    events = flow['charge/events/data']
    ids = np.unique(tracks['event_id']).astype('i8')
    row_ids = bool(np.all(ids >= 0) and np.all(ids < len(events)))
    rows = (events[ids] if len(ids) else events[:0]) if row_ids else events[:]
    # Read only referenced event chunks; their IDs verify the row assumption.
    if not np.array_equal(rows['id'], ids):
        rows = events[:]
        if len(np.unique(rows['id'])) != len(rows): raise ValueError('duplicate FLOW event IDs')
    lookup = {int(r['id']): int(r['n_ext_trigs']) > 0 for r in rows}
    if not set(map(int, ids)).issubset(lookup): raise ValueError('invalid track event ID')
    # The association is event -> selected track -> segment, not segment row
    # position or a direction cut. Multiple tracks in one event share its flag;
    # an externally triggered event may itself contain a physical cosmic muon.
    external = np.empty(len(segments), dtype=bool)
    external[refs[:,1]] = [lookup[int(tracks[i]['event_id'])] for i in refs[:,0]]
    return {'beam': segments[external], 'cosmic': segments[~external]}, 'legacy FLOW selection partitioned by event n_ext_trigs', dict(
        original_segment_sha256=hashlib.sha256(segments.tobytes()).hexdigest(),
        source_rows={'beam': np.flatnonzero(external).tolist(), 'cosmic': np.flatnonzero(~external).tolist()})


def write_shard(flow_path, directory, identity, metadata):
    """Read primary FLOW and publish its selected arrays without modifying it."""
    with h5py.File(flow_path, 'r') as flow:
        arrays, origin, links = stored_segments(flow)
        signature = reconstruction_signature(flow)
        signature['stored_selections'] = selection_metadata(flow)
    return store_arrays(directory, identity, metadata, arrays, signature, origin, links)


def store_arrays(directory, identity, metadata, arrays, signature, origin, links):
    """Publish immutable arrays followed by their small descriptor.

    Also supports replay of hash-verified primary-FLOW captures. Such callers
    must validate their acquisition receipts before entering here and include
    that provenance in metadata. Published lifetime values are never inputs.
    """
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    if signature['run_info']['crs_ticks'] != .1: raise ValueError('legacy fitter requires 0.1-us ticks')
    signature = dict(signature)  # Do not mutate a caller's provenance record.
    hashes = {s: hashlib.sha256(a.tobytes()).hexdigest() for s,a in arrays.items()}
    signature.update(selection_origin=origin, segment_dtypes={s:str(a.dtype) for s,a in arrays.items()})
    payload_id = digest_json(hashes)
    filename = identity+'-'+payload_id+'.npz'
    target = directory/filename
    if not target.exists():
        fd, temp = tempfile.mkstemp(dir=directory, suffix='.npz'); os.close(fd)
        try:
            np.savez_compressed(temp, **arrays); os.replace(temp, target)
        finally:
            if os.path.exists(temp): os.unlink(temp)
    row = dict(metadata, sample_id=identity, arrays_file=filename, array_hashes=hashes,
               configuration=signature, configuration_id=digest_json(signature), selection_origin=origin, **links)
    atomic_json(directory/(identity+'.sample.json'), row)
    return row


def pool_tracks(directory, output, hours=6, timezone='America/Chicago'):
    """Fit each compatible calendar pool once per exact membership/code version.

    Caller holds the monitor's serialization lock while reading descriptors,
    fitting and publishing. Array hashes are checked even when a fit is cached.
    Thus retries cannot multiply weights or resurrect a retired configuration.
    """
    directory, output = Path(directory), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    for path in sorted(directory.glob('*.sample.json')):
        r = json.loads(path.read_text())
        groups[r['cohort'], r['configuration_id'], window_start(r['timestamp'], hours, timezone)].append(r)
    code = hashlib.sha256(b''.join(Path(__file__).with_name(p).read_bytes() for p in
        ('track_shards.py', 'lifetime.py', 'lifetime_funcs.py'))).hexdigest()
    result = []
    for (cohort, configuration, start), members in sorted(groups.items()):
        fingerprint = digest_json(dict(members=members, code=code, hours=hours, timezone=timezone))
        folder = output/fingerprint; folder.mkdir(exist_ok=True)
        arrays, origins = defaultdict(list), defaultdict(list)
        for i, member in enumerate(members):
            with np.load(directory/member['arrays_file'], allow_pickle=False) as source:
                for sample, expected in member['array_hashes'].items():
                    array = source[sample]
                    if hashlib.sha256(array.tobytes()).hexdigest() != expected: raise ValueError('track shard hash mismatch')
                    arrays[sample].append(array); origins[sample].append(np.full(len(array), i, dtype='i4'))
        cache = folder/'results.json'
        if cache.exists():
            result.extend(json.loads(cache.read_text())); continue
        pooled = folder/'segments.h5'
        with h5py.File(pooled, 'w') as h:
            for sample, pieces in arrays.items():
                path = dict(SPLIT_SAMPLES).get(sample, 'analysis/rock_muon_segments/data')
                h.create_dataset(path, data=np.concatenate(pieces))
                h.create_dataset(str(Path(path).parent)+'/member_index', data=np.concatenate(origins[sample]))
            h.attrs.update(complete=True, pool_members=json.dumps(members))
        end = (start.tz_localize(None)+timedelta(hours=hours)).tz_localize(timezone)
        fitted = fit_samples(str(pooled), str(folder/'fit.png'), timestamp=(start+(end-start)/2).isoformat())
        for row in fitted:
            row.update(method='track_pool', source=cohort+'/'+configuration[:12], cohort=cohort,
                window_hours=hours, period_start=start.isoformat(), period_end=end.isoformat(),
                n_files=len(members), members=members, configuration_id=configuration,
                first_observed_at=min(m['timestamp'] for m in members), last_observed_at=max(m['timestamp'] for m in members),
                monitoring_excluded=any(m.get('monitoring_excluded', False) for m in members),
                aggregation='selected FLOW segments pooled before fitting; no averaging of file lifetimes',
                selection_origin=members[0]['selection_origin'], fit_code_sha256=code)
        atomic_json(cache, fitted); result.extend(fitted)
    return result
