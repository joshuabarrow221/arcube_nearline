"""Veto noisy stored FLOW segments using their original event/hit references.

Surviving segment rows are copied byte-for-byte: removing individual hits would
change charge and path length differently and bias the lifetime estimator.
This adapter cannot recover tracks omitted by the upstream through-going cut.
"""
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

from packet_pedestal import take
from track_quality import TrackQualityConfig, quality_hit_mask


HITS = 'charge/calib_prompt_hits'
EVENTS = 'charge/events'
BOUNDS_TOLERANCE_CM = 2.1


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def _json(value):
    if isinstance(value, dict):
        return {key: _json(item) for key, item in value.items()}
    if isinstance(value, (np.ndarray, np.generic)):
        return _json(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_json(item) for item in value]
    return value.decode() if isinstance(value, bytes) else value


@dataclass
class QualityContext:
    """Geometry and timing loaded once; evaluate the original hits of one event."""
    config: TrackQualityConfig
    tile_data: np.ndarray
    tile_keys: np.ndarray
    module_bounds: np.ndarray
    metadata: dict

    def evaluate(self, hits):
        # H5Flow LUT: row zero is the default, first key has unit stride.
        keys = np.column_stack((hits['io_group'], hits['io_channel'])).astype('i8')
        valid_keys = np.all((keys >= self.tile_keys[:, 0]) &
                            (keys <= self.tile_keys[:, 1]), axis=1)
        lengths = self.tile_keys[:, 1] - self.tile_keys[:, 0] + 1
        indices = 1 + keys[:, 0] - self.tile_keys[0, 0]
        indices += (keys[:, 1] - self.tile_keys[1, 0]) * lengths[0]
        indices[~valid_keys] = 0
        tiles = self.tile_data['data'][indices].astype('i8')
        mapped = valid_keys & self.tile_data['filled'][indices] & (tiles >= 0)
        # Explicit 2x2 readout convention, also used by FLOW's drift velocity
        # lookup: adjacent io_groups are the two TPCs of one module.
        modules = (hits['io_group'].astype('i8') - 1) // 2
        mapped &= (modules >= 0) & (modules < len(self.module_bounds))
        keep = np.zeros(len(hits), dtype=bool)
        local_keep, diagnostics = quality_hit_mask(
            hits[mapped], self.config, tile_ids=tiles[mapped], module_ids=modules[mapped])
        keep[mapped] = local_keep
        xyz = np.column_stack([hits[name] for name in ('x', 'y', 'z')])
        inside = np.zeros(len(hits), dtype=bool)
        bounds = self.module_bounds[modules[mapped]]
        inside[mapped] = np.all(np.isfinite(xyz[mapped]) &
            (xyz[mapped] >= bounds[:, 0] - BOUNDS_TOLERANCE_CM) &
            (xyz[mapped] <= bounds[:, 1] + BOUNDS_TOLERANCE_CM), axis=1)
        keep &= inside
        diagnostics.update(n_input_hits=len(hits), n_kept_hits=int(keep.sum()),
            n_unmapped_hits=int((~mapped).sum()),
            n_outside_module_hits=int((mapped & ~inside).sum()),
            bounds_tolerance_cm=BOUNDS_TOLERANCE_CM)
        return keep, diagnostics


def build_quality_context(source, quality_config=None):
    """Use source clock/geometry metadata, rejecting conflicting overrides."""
    values = (quality_config.to_dict() if isinstance(quality_config, TrackQualityConfig)
              else dict(quality_config or {}))
    run = source['run_info'].attrs
    geometry = source['geometry_info'].attrs
    authoritative = dict(tick_us=float(run['crs_ticks']),
                         drift_axis=_text(geometry.get('drift_direction', 'x')))
    if 'rollover_ticks' in run:
        authoritative['rollover_ticks'] = float(run['rollover_ticks'])
    for key, value in authoritative.items():
        if key in values and values[key] is not None and values[key] != value:
            raise ValueError(f'quality {key} conflicts with FLOW metadata')
        values[key] = value
    config = TrackQualityConfig(**values)
    group = source['geometry_info/tile_id']
    keys = np.asarray(group.attrs['meta'][0]['min_max_keys'], dtype='i8')
    tiles = group['data'][:]
    if keys.shape != (2, 2) or np.any(keys[:, 1] < keys[:, 0]):
        raise ValueError('invalid FLOW tile lookup geometry')
    expected = 1 + int(np.prod(keys[:, 1] - keys[:, 0] + 1))
    if len(tiles) != expected or not {'data', 'filled'} <= set(tiles.dtype.names or ()):
        raise ValueError('invalid FLOW tile lookup payload')
    bounds = np.asarray(geometry['module_RO_bounds'], dtype=float)
    if (bounds.shape != (4, 2, 3) or not np.isfinite(bounds).all()
            or np.any(bounds[:, 0] >= bounds[:, 1])):
        raise ValueError('quality selection requires four 2x2 module_RO_bounds in xyz cm')
    quality_path = Path(__file__).with_name('track_quality.py')
    metadata = dict(quality_config=json.dumps(config.to_dict(), sort_keys=True),
        quality_code_sha256=hashlib.sha256(quality_path.read_bytes()).hexdigest(),
        geometry_quality_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        tick_period_us=config.tick_us,
        geometry_bounds_policy='per-module xyz bounds plus 2.1 cm tolerance',
        module_mapping='2x2 paired io_groups: module=(io_group-1)//2')
    return QualityContext(config, tiles, keys, bounds, metadata)


def _reference(source, first, second):
    path = f'{first}/ref/{second}/ref'
    if path not in source:
        raise ValueError(f'quality selection requires reference {path}')
    node = source[path]
    if (_text(node.attrs.get('dset0', '')).lstrip('/') != first + '/data'
            or _text(node.attrs.get('dset1', '')).lstrip('/') != second + '/data'):
        raise ValueError(f'unexpected reference orientation: {path}')
    if len(node.shape) != 2 or node.shape[1] != 2 or node.dtype.kind not in 'iu':
        raise ValueError(f'invalid reference array: {path}')
    return node


def _parents(refs, child_count, parent_count, label):
    if (refs.shape != (child_count, 2) or
            not np.array_equal(np.sort(refs[:, 1]), np.arange(child_count))):
        raise ValueError(f'each {label} must have exactly one parent reference')
    if np.any(refs[:, 0] < 0) or np.any(refs[:, 0] >= parent_count):
        raise ValueError(f'invalid parent reference for {label}')
    result = np.empty(child_count, dtype='i8')
    result[refs[:, 1]] = refs[:, 0]
    return result


def _settings(source, path):
    return {location: _json({key: value for key, value in node.attrs.items()
            if key not in ('selected_track_count', 'selected_segment_count')
            and not key.startswith('electron_lifetime')})
        for location, node in [('group', source[path]), ('data', source[path + '/data'])]}


def select_stored_muons(input_file, output_file, quality_config=None):
    """Copy clean stored segments to a complete, separate timing-split sidecar.

    Every retained segment must have full hit membership and one parent track
    and event. One rejected member hit vetoes the whole segment. Event-wide
    occupancy is evaluated before filtering, including hits outside the track.
    ``events_processed`` denotes coverage of the source's stored selection;
    the quality/events table lists the selected parent events actually read.
    """
    input_file, output_file = Path(input_file), Path(output_file)
    if input_file.resolve() == output_file.resolve():
        raise ValueError('quality sidecar must not replace source FLOW')
    arrays, origins, track_counts, audits = {'beam': [], 'cosmic': []}, {'beam': [], 'cosmic': []}, {'beam': 0, 'cosmic': 0}, []
    with h5py.File(input_file, 'r') as source:
        context = build_quality_context(source, quality_config)
        events, hits = source[EVENTS + '/data'], source[HITS + '/data']
        event_hits_ref = _reference(source, EVENTS, HITS)
        region_path = f'{EVENTS}/ref/{HITS}/ref_region'
        if region_path not in source or len(source[region_path]) != len(events):
            raise ValueError('quality selection requires complete event-to-hit ref_region')
        regions = source[region_path]
        layouts = []
        split = ['analysis/beam_rock_muon_segments', 'analysis/cosmic_muon_segments']
        if any(path + '/data' in source for path in split):
            if not all(path + '/data' in source for path in split):
                raise ValueError('incomplete pair of stored timing selections')
            layouts = [(path, sample) for path, sample in zip(split, ['beam', 'cosmic'])]
        elif 'analysis/rock_muon_segments/data' in source:
            layouts = [('analysis/rock_muon_segments', None)]
        else:
            raise ValueError('FLOW has no stored muon segments')
        selector_settings, output_dtype = {}, None
        # Cache only event masks shared across split selections; a stored segment
        # carries far fewer rows than the complete FLOW hit collection.
        event_cache = {}
        for segment_path, declared_sample in layouts:
            track_path = segment_path.replace('_segments', '_tracks')
            if track_path + '/data' not in source:
                raise ValueError('quality selection requires original parent tracks')
            segments, tracks = source[segment_path + '/data'][:], source[track_path + '/data'][:]
            if output_dtype is not None and segments.dtype != output_dtype:
                raise ValueError('stored timing selections have different segment dtypes')
            output_dtype = segments.dtype
            selector_settings[segment_path] = dict(segments=_settings(source, segment_path),
                tracks=_settings(source, track_path), track_dtype=_json(tracks.dtype.descr))
            parent_tracks = _parents(_reference(source, track_path, segment_path)[:], len(segments), len(tracks), 'segment')
            parent_events = _parents(_reference(source, EVENTS, track_path)[:], len(tracks), len(events), 'track')
            event_rows = np.unique(parent_events)
            records = take(events, event_rows)
            event_lookup = {int(row): record for row, record in zip(event_rows, records)}
            for track_row, event_row in enumerate(parent_events):
                if int(tracks[track_row]['event_id']) != int(event_lookup[int(event_row)]['id']):
                    raise ValueError('stored track event_id disagrees with event reference')
            links = _reference(source, segment_path, HITS)[:].astype('i8')
            if (np.any(links < 0) or np.any(links[:, 0] >= len(segments))
                    or np.any(links[:, 1] >= len(hits))):
                raise ValueError('invalid segment-to-hit reference')
            if len(np.unique(links, axis=0)) != len(links):
                raise ValueError('duplicate segment-to-hit reference')
            counts = np.bincount(links[:, 0], minlength=len(segments))
            if not np.array_equal(counts, segments['nhits']) or np.any(counts == 0):
                raise ValueError('stored segment nhits disagrees with complete hit references')
            order = np.argsort(links[:, 0], kind='stable')
            links = links[order]
            offsets = np.r_[0, np.cumsum(counts)]
            keep = np.ones(len(segments), dtype=bool)
            samples = np.empty(len(segments), dtype='U6')
            for event_row in event_rows:
                event_row = int(event_row)
                event = event_lookup[event_row]
                sample = 'beam' if event['n_ext_trigs'] > 0 else 'cosmic'
                if declared_sample is not None and sample != declared_sample:
                    raise ValueError('stored timing selection disagrees with parent trigger category')
                if event_row not in event_cache:
                    region = regions[event_row]
                    start, stop = int(region['start']), int(region['stop'])
                    if start < 0 or stop < start or stop > len(event_hits_ref):
                        raise ValueError('invalid event-to-hit reference region')
                    event_links = event_hits_ref[start:stop]
                    hit_rows = event_links[event_links[:, 0] == event_row, 1].astype('i8')
                    if (np.any(hit_rows < 0) or np.any(hit_rows >= len(hits))
                            or len(np.unique(hit_rows)) != len(hit_rows)):
                        raise ValueError('invalid or duplicate event-to-hit reference')
                    hit_rows.sort()
                    event_hits = take(hits, hit_rows)
                    hit_keep, diagnostics = context.evaluate(event_hits)
                    event_cache[event_row] = (hit_rows, hit_keep, event_hits['io_group'], diagnostics)
                hit_rows, hit_keep, hit_groups, diagnostics = event_cache[event_row]
                segment_rows = np.flatnonzero(parent_events[parent_tracks] == event_row)
                samples[segment_rows] = sample
                for segment_row in segment_rows:
                    members = links[offsets[segment_row]:offsets[segment_row + 1], 1]
                    positions = np.searchsorted(hit_rows, members)
                    if (np.any(positions >= len(hit_rows)) or
                            not np.array_equal(hit_rows[positions], members)):
                        raise ValueError('segment hit lies outside its parent event references')
                    if ('io_group' in segments.dtype.names and
                            np.any(hit_groups[positions] != segments[segment_row]['io_group'])):
                        raise ValueError('segment io_group disagrees with its hit references')
                    keep[segment_row] = bool(np.all(hit_keep[positions]))
                audits.append(dict(diagnostics, event_row=event_row, event_id=int(event['id']),
                    source_segments=segment_path, n_segments=len(segment_rows),
                    n_kept_segments=int(keep[segment_rows].sum())))
            for sample in ('beam', 'cosmic'):
                selected = np.flatnonzero(keep & (samples == sample))
                arrays[sample].append(segments[selected])
                origins[sample].extend((segment_path, int(row)) for row in selected)
                track_counts[sample] += len(np.unique(parent_tracks[selected]))
        n_events = len(events)
        metadata = context.metadata
    settings_json = json.dumps(selector_settings, sort_keys=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=output_file.parent, suffix='.h5.tmp')
    os.close(fd)
    try:
        with h5py.File(temporary, 'w') as output:
            output.attrs.update(metadata, complete=True, source_file=str(input_file.resolve()),
                events_processed=n_events,
                selector_sha256=hashlib.sha256(settings_json.encode()).hexdigest(),
                selector_identity_kind='recorded stored-selection metadata; deployed source code unknown',
                adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                selection_origin='stored FLOW segments vetoed using original event/hit electronics evidence')
            for sample, name in [('beam', 'beam_rock_muon_segments'), ('cosmic', 'cosmic_muon_segments')]:
                group = output.create_group('analysis/' + name)
                array = np.concatenate(arrays[sample]) if arrays[sample] else np.empty(0, dtype=output_dtype)
                group.create_dataset('data', data=array)
                group.create_dataset('source_row', data=np.array([row for _, row in origins[sample]], dtype='i8'))
                group.create_dataset('source_dataset', data=np.array([path for path, _ in origins[sample]], dtype=h5py.string_dtype()))
                group.attrs.update(selected_track_count=track_counts[sample], selected_segment_count=len(array),
                    event_timing='beam' if sample == 'beam' else 'off_beam',
                    source_selector_settings=settings_json,
                    quality_selection='whole-segment veto; surviving dQ/dx and t unchanged')
            fields = ['event_row', 'event_id', 'n_input_hits', 'n_kept_hits', 'n_segments', 'n_kept_segments', 'n_unmapped_hits', 'n_outside_module_hits']
            table = np.array([tuple(row[field] for field in fields) for row in audits], dtype=[(field, 'i8') for field in fields])
            output.create_dataset('quality/events/data', data=table)
            output.create_dataset('quality/events/diagnostics', data=np.array([json.dumps(row, sort_keys=True) for row in audits], dtype=h5py.string_dtype()))
        os.replace(temporary, output_file)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return str(output_file)
