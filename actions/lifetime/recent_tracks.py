#!/usr/bin/env python3
"""Select and jointly fit FLOW segments from a trailing acquisition-time window.

File acquisition start times define membership; segment ``t`` is drift time,
never wall-clock time. Source files are opened read-only. Each invocation gets
a new output directory so a failed rerun cannot leave an apparently current
plot or accidentally ingest its own previous pool.
"""
import argparse
from contextlib import ExitStack, contextmanager
from datetime import datetime, timedelta, timezone
import glob
import hashlib
import json
import logging
from pathlib import Path
import re
import sys
import tempfile

import h5py
import numpy as np

from lifetime import LEGACY_SAMPLE, SPLIT_SAMPLES, fit_samples, valid_segments
from purity_sources import atomic_json, aware_time
from track_pooling import json_value, reconstruction_signature
from track_shards import recorded_settings, selection_metadata


QUALITY_IDENTITY = ('selector_sha256', 'adapter_sha256', 'quality_code_sha256',
                    'geometry_quality_code_sha256', 'quality_config', 'tick_period_us', 'selection_origin')
REQUIRED_FIELDS = ('dx', 'dQ', 'nhits', 't')


class _RepeatedFitWarnings(logging.Filter):
    """Show each pylandau message once while retaining its full frequency."""
    def __init__(self):
        super().__init__()
        self.records = []
        self._messages = {}

    def filter(self, record):
        message = record.getMessage()
        key = (record.levelno, message)
        if key in self._messages:
            self._messages[key]['occurrences'] += 1
            return False
        entry = dict(logger=record.name, level=record.levelname,
                     message=message, occurrences=1)
        self._messages[key] = entry
        self.records.append(entry)
        return True


@contextmanager
def fit_warning_receipt():
    """Temporarily deduplicate this dependency's logs, including failed fits.

    The first occurrence keeps its original logger, severity and handlers.
    Python/scipy warnings and all other loggers are unaffected.
    """
    logger = logging.getLogger('pylandau')
    capture = _RepeatedFitWarnings()
    logger.addFilter(capture)
    try:
        yield capture.records
    finally:
        logger.removeFilter(capture)


def discover_inputs(inputs, output=None):
    """Resolve files, recursive directories and quoted globs exactly once."""
    found = {}
    output = Path(output).resolve() if output is not None else None
    for specification in inputs:
        expanded = glob.glob(str(specification), recursive=True)
        if not expanded:
            raise FileNotFoundError(f'no inputs match {specification}')
        for entry in expanded:
            entry = Path(entry)
            candidates = (p for p in entry.rglob('*') if p.is_file()
                          and p.suffix.lower() in ('.h5', '.hdf5', '.flow')) if entry.is_dir() else [entry]
            for path in candidates:
                path = path.resolve()
                if output is not None and path.is_relative_to(output):
                    continue
                stat = path.stat()
                # Resolve symlinks and hard links, but retain distinct files
                # with repeated file-local segment IDs.
                identity = (stat.st_dev, stat.st_ino) if stat.st_ino else str(path)
                found.setdefault(identity, path)
    return sorted(found.values())


def filename_time(name):
    """Parse the recorded acquisition offset, including DST folds correctly."""
    match = re.search(r'(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})_(CST|CDT|UTC|CET)(?=[^A-Za-z]|$)', Path(str(name)).name)
    if not match:
        raise ValueError(f'no acquisition timestamp in filename: {name}')
    offsets = {'CST': '-0600', 'CDT': '-0500', 'UTC': '+0000', 'CET': '+0100'}
    return aware_time(datetime.strptime(match[1] + offsets[match[2]], '%Y_%m_%d_%H_%M_%S%z'))


def source_identity(path, source):
    """A sidecar and its source must not acquire two weights in one pool."""
    name = json_value(source.attrs.get('source_file', source.attrs.get('source_flow_url', str(path))))
    if str(name).startswith(('http://', 'https://')):
        return str(name)
    name = Path(name)
    if not name.is_absolute():
        name = path.parent/name
    return str(name.resolve())


def acquisition_time(path, source):
    identity = source_identity(path, source)
    try:
        return filename_time(identity), identity
    except ValueError:
        return filename_time(path), str(path)


def segment_paths(source, requested=None):
    if requested:
        return [requested] if requested in source else []
    split = [path for _, path in SPLIT_SAMPLES if path in source]
    if split:
        if len(split) != 2:
            raise ValueError('incomplete pair of stored timing selections; specify --dataset explicitly to fit one category')
        # Reuse validation of the split labels/settings without requiring
        # event references or copying legacy arrays simply to discover paths.
        settings = selection_metadata(source)
        comparable = []
        for sample, expected in (('beam', 'beam'), ('cosmic', 'off_beam')):
            locations = {}
            for location, attrs in settings[sample].items():
                if attrs.get('event_timing', expected) != expected:
                    raise ValueError(f'unexpected event_timing for {sample} selection')
                locations[location] = {k: v for k, v in attrs.items() if k != 'event_timing'}
            comparable.append(locations)
        if comparable[0] != comparable[1]:
            raise ValueError('split FLOW selections have different recorded settings')
        return split
    return [LEGACY_SAMPLE[1]] if LEGACY_SAMPLE[1] in source else []


def pool_signature(source, selected, source_path, stack, allow_missing=False):
    """Compare recorded reconstruction/selection and require 0.1-us ticks."""
    original = source
    if 'geometry_info' not in original and 'source_file' in source.attrs:
        recorded = Path(source_identity(source_path, source))
        if recorded.is_file() and recorded.resolve() != source_path.resolve():
            original = stack.enter_context(h5py.File(recorded, 'r'))
    try:
        reconstruction = reconstruction_signature(original)
    except KeyError as exc:
        if not allow_missing:
            raise ValueError('missing reconstruction metadata; use full FLOW inputs or explicitly '
                             'acknowledge legacy metadata with --allow-missing-metadata') from exc
        # Preserve every available relevant attribute; missing metadata must
        # not erase a *known* disagreement or a nonstandard tick period.
        reconstruction = {name: json_value(dict(original[name].attrs)) if name in original else None
                          for name in ('geometry_info', 'lar_info', 'charge/calib_prompt_hits', 'run_info')}
    ticks = reconstruction.get('run_info') or {}
    tick_period = ticks.get('crs_ticks', selected.attrs.get('tick_period_us'))
    if tick_period is None and not allow_missing:
        raise ValueError('missing drift tick period')
    if tick_period is not None and float(tick_period) != .1:
        raise ValueError('the existing lifetime fitter requires 0.1-us drift ticks')
    return dict(reconstruction=reconstruction, selections=selection_metadata(selected),
                quality={key: json_value(selected.attrs[key]) for key in QUALITY_IDENTITY if key in selected.attrs},
                missing_metadata_acknowledged=allow_missing)


def gather_recent(paths, start, end, directory, selection='quality-stored', dataset=None,
                  quality_config=None, allow_missing_metadata=False):
    """Gather raw selected arrays with parallel row lineage and a receipt.

    Bad files are recorded and cause the caller to stop before fitting. Missing
    requested datasets are explicit skips, matching the notebook gather loop.
    Invalid numerical segment rows are removed once after concatenation.
    """
    arrays, members, skipped, errors = [], [], [], []
    member_rows, dataset_rows, source_rows = [], [], []
    signature, dtype, seen_sources, payloads = None, None, set(), set()
    for input_index, path in enumerate(paths):
        path = Path(path).resolve()
        selected_path = None
        try:
            with ExitStack() as stack:
                source = stack.enter_context(h5py.File(path, 'r'))
                if 'pool_members' in source.attrs:
                    skipped.append(dict(path=str(path), reason='already pooled output'))
                    continue
                stamp, time_source = acquisition_time(path, source)
                if not start <= stamp < end:
                    skipped.append(dict(path=str(path), reason='outside acquisition-time window', timestamp=stamp.isoformat()))
                    continue
                identity = source_identity(path, source)
                if identity in seen_sources:
                    skipped.append(dict(path=str(path), reason='duplicate source acquisition', source=identity))
                    continue
                if not bool(source.attrs.get('complete', True)):
                    raise ValueError('incomplete selected-segment input')
                selected, selected_path = source, path
                if selection != 'stored':
                    required = ('charge/events/data', 'charge/calib_prompt_hits/data', 'geometry_info')
                    if any(key not in source for key in required):
                        raise ValueError('fresh quality selection needs full FLOW hits/events/geometry; '
                                         '--selection stored explicitly opts into unverified upstream selections')
                    if selection == 'reselect':
                        from track_selection import select_muons
                    else:
                        from stored_quality import select_stored_muons as select_muons
                    # A completed sidecar can still fail compatibility checks.
                    # Keep its evidence instead of reusing its name for the
                    # next successfully gathered source.
                    selected_path = directory/'selections'/f'{input_index:05d}.segments.h5'
                    select_muons(str(path), str(selected_path), quality_config=quality_config)
                    selected = stack.enter_context(h5py.File(selected_path, 'r'))
                    if not bool(selected.attrs.get('complete', False)):
                        raise ValueError('fresh selection did not complete')
                keys = segment_paths(selected, dataset)
                if not keys:
                    skipped.append(dict(path=str(path), reason='missing requested segment dataset' if dataset else 'no stored segment dataset'))
                    continue
                current = pool_signature(source, selected, path, stack, allow_missing_metadata)
                current['selected_dataset_settings'] = {key: recorded_settings(selected, key) for key in keys}
                if selection == 'quality-stored':
                    # Vetted sidecars can have only counters on their output
                    # groups. The original selector settings still determine
                    # which tracks were available before the quality veto.
                    current['upstream_selections'] = selection_metadata(source)
                if signature is not None and current != signature:
                    raise ValueError('incompatible reconstruction or selection metadata; fit separate acquisition cohorts')
                pieces, records = [], []
                for key in keys:
                    node = selected[key]
                    if not isinstance(node, h5py.Dataset) or node.ndim != 1:
                        raise ValueError(f'{key} must be a one-dimensional segment dataset')
                    data = node[:]
                    if not set(REQUIRED_FIELDS).issubset(data.dtype.names or ()):
                        raise ValueError(f'{key} lacks required fields {REQUIRED_FIELDS}')
                    if dtype is not None and data.dtype != dtype:
                        raise ValueError('incompatible segment dtypes')
                    dtype = data.dtype
                    pieces.append(data)
                    records.append(dict(path=key, n_segments=len(data), sha256=hashlib.sha256(data.tobytes()).hexdigest()))
                # Renamed copies have different path identities. Reject exact
                # nonempty payload repeats rather than silently double weight.
                digest = hashlib.sha256(json.dumps(records, sort_keys=True).encode()).hexdigest()
                if sum(len(a) for a in pieces) and digest in payloads:
                    raise ValueError('duplicate nonempty segment payload, possibly a renamed input copy')
                payloads.add(digest)
                index = len(members)
                for key_index, data in enumerate(pieces):
                    arrays.append(data)
                    member_rows.append(np.full(len(data), index, dtype='i4'))
                    dataset_rows.append(np.full(len(data), key_index, dtype='i4'))
                    source_rows.append(np.arange(len(data), dtype='i8'))
                members.append(dict(input_file=str(path), source=identity, selection_file=str(selected_path),
                                    timestamp=stamp.isoformat(), timestamp_source=time_source, datasets=records,
                                    n_raw_segments=sum(len(a) for a in pieces)))
                if 'quality/events' in selected:
                    event_table_path = ('quality/events' if isinstance(selected['quality/events'], h5py.Dataset)
                                        else 'quality/events/data')
                    diagnostic_path = ('quality/events/diagnostics' if 'quality/events/diagnostics' in selected
                                       else event_table_path)
                    members[-1]['quality_diagnostics'] = dict(file=str(selected_path), dataset=diagnostic_path,
                                                              n_events=len(selected[event_table_path]))
                seen_sources.add(identity)
                signature = current
                print(f'+{members[-1]["n_raw_segments"]:6d} segments  {path.name}', flush=True)
        except Exception as exc:
            failure = dict(path=str(path), error=f'{type(exc).__name__}: {exc}')
            if selected_path is not None and selected_path != path and selected_path.exists():
                failure['selection_file'] = str(selected_path)
            errors.append(failure)
            print(f'Input error: {path}: {failure["error"]}', file=sys.stderr, flush=True)
    report = dict(members=members, skipped=skipped, errors=errors, reconstruction_signature=signature,
                  n_files=len(members), n_raw_segments=sum(len(a) for a in arrays),
                  n_valid_segments=0, n_invalid_segments=0)
    if not arrays:
        return None, None, report
    raw = np.concatenate(arrays)
    # Keep row provenance by applying the same numerical mask as the fitter.
    valid_rows = np.flatnonzero(np.logical_and.reduce(
        [np.isfinite(raw[name]) for name in REQUIRED_FIELDS]
        + [raw[name] > 0 for name in ('dx', 'dQ', 'nhits')]))
    segments = valid_segments(raw)
    lineage = dict(member_index=np.concatenate(member_rows)[valid_rows],
                   source_dataset_index=np.concatenate(dataset_rows)[valid_rows],
                   source_row=np.concatenate(source_rows)[valid_rows])
    report.update(n_valid_segments=len(segments), n_invalid_segments=len(raw)-len(segments))
    return segments, lineage, report


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('inputs', nargs='+', help='FLOW files, directories (recursive), or quoted glob patterns')
    result.add_argument('--output', required=True, type=Path, help='parent directory for a new, isolated result directory')
    result.add_argument('--hours', type=float, default=6., help='trailing elapsed hours (default: 6)')
    result.add_argument('--end', help='exclusive window end as ISO timestamp with UTC offset (default: current time)')
    result.add_argument('--selection', choices=('quality-stored', 'reselect', 'stored'), default='quality-stored',
                        help='hit-vet stored segments (default), reselect from hits, or explicitly trust stored upstream segments')
    result.add_argument('--reselect', action='store_const', const='reselect', dest='selection', help=argparse.SUPPRESS)
    result.add_argument('--dataset', help='one HDF5 segment dataset; only with --selection stored')
    result.add_argument('--quality-config', type=Path, help='JSON object passed to hit quality selection')
    result.add_argument('--allow-missing-metadata', action='store_true',
                        help='acknowledge incomplete legacy metadata and assume 0.1-us ticks only if absent')
    return result


def main(argv=None):
    options = parser().parse_args(argv)
    if not np.isfinite(options.hours) or options.hours <= 0:
        parser().error('--hours must be finite and positive')
    if options.dataset and options.selection != 'stored':
        parser().error('--dataset requires --selection stored')
    if options.quality_config and options.selection == 'stored':
        parser().error('--quality-config requires hit quality selection')
    try:
        end = aware_time(options.end or datetime.now(timezone.utc))
        start = end-timedelta(hours=options.hours)
        paths = discover_inputs(options.inputs, options.output)
        config = json.loads(options.quality_config.read_text()) if options.quality_config else None
        if config is not None and not isinstance(config, dict):
            raise ValueError('--quality-config must contain a JSON object')
    except (ValueError, OSError, OverflowError) as exc:
        parser().error(str(exc))
    options.output.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix='recent-', dir=options.output)).resolve()
    receipt = dict(period_start=start.isoformat(), period_end=end.isoformat(), window_hours=options.hours,
                   membership='file acquisition start in [period_start, period_end); all selected segments from each included file',
                   selection=options.selection,
                   noise_quality=('fresh hit quality selection' if options.selection == 'reselect' else
                                  'stored segments vetted through original hit references' if options.selection == 'quality-stored' else
                                  'unverified upstream selection; stored segments cannot establish hit-level noise quality'),
                   metadata_limitation='equal recorded metadata does not verify stable live hardware settings',
                   drift_time_unit='0.1 us tick; 0..1960 ticks corresponds to 0..196 us',
                   fit_status='unavailable', lifetime_us=None, error_us=None, fit_warnings=[],
                   output_directory=str(directory))
    segments, lineage, gathered = gather_recent(paths, start, end, directory, options.selection,
        options.dataset.lstrip('/') if options.dataset else None, config, options.allow_missing_metadata)
    receipt.update(gathered)
    status = 2
    if receipt['errors']:
        receipt['fit_message'] = 'input errors prevent a partial-input lifetime fit; inspect errors and rerun corrected inputs'
    elif segments is None or not len(segments):
        receipt['fit_message'] = 'no valid segments in the requested acquisition-time window'
    else:
        pooled = directory/'segments.h5'
        with h5py.File(pooled, 'w') as output:
            output.create_dataset(LEGACY_SAMPLE[1], data=segments)
            group = output[str(Path(LEGACY_SAMPLE[1]).parent)]
            for name, values in lineage.items():
                group.create_dataset(name, data=values)
            output.attrs.update(complete=True, pool_members=json.dumps(receipt['members']),
                                reconstruction_signature=json.dumps(receipt['reconstruction_signature'], sort_keys=True),
                                period_start=start.isoformat(), period_end=end.isoformat(), selection=options.selection)
            output.attrs.update(receipt['reconstruction_signature']['quality'])
        with fit_warning_receipt() as warnings:
            receipt['fit_warnings'] = warnings
            result = fit_samples(str(pooled), str(directory/'fit.png'), timestamp=(start+(end-start)/2).isoformat())[0]
        receipt.update(fit=result, pooled_segments=str(pooled), fit_status=result['fit_status'],
                       lifetime_us=result['lifetime_us'], error_us=result['error_us'])
        status = 0 if result['fit_status'] == 'ok' else 2
    atomic_json(directory/'summary.json', receipt)
    message = receipt.get('fit_message') or receipt.get('fit', {}).get('fit_message')
    if message:
        print(f'Fit {receipt["fit_status"]}: {message}', file=sys.stderr, flush=True)
    error_count = len(receipt['errors'])
    print(f'{receipt["n_files"]} files, {receipt.get("n_valid_segments", 0)} valid segments; '
          f'fit {receipt["fit_status"]}; {error_count} input error{"" if error_count == 1 else "s"}')
    print(f'Receipt: {directory / "summary.json"}')
    return status


if __name__ == '__main__':
    raise SystemExit(main())
