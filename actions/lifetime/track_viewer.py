#!/usr/bin/env python3
"""Export receipted lifetime parent tracks to an offline browser display.

The summary is an index, not a hit store. Run this beside the original FLOW,
selection sidecars, and pooled segments on NERSC, then copy the HTML anywhere.
No fits or quality cuts are rerun. All source files are opened read-only.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np

from packet_pedestal import take
from stored_quality import EVENTS, HITS, _parents, _reference


POOL = 'analysis/rock_muon_segments'
FIELDS = ('dx', 'dQ', 'nhits', 't')


def _text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def _clean(value):
    """Strict JSON: missing/nonfinite physical values stay null, never zero."""
    if isinstance(value, dict):
        return {str(k): _clean(v) for k, v in value.items()}
    if isinstance(value, (np.ndarray, np.generic)):
        return _clean(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_clean(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value.decode() if isinstance(value, bytes) else value


def _valid(rows):
    if not set(FIELDS) <= set(rows.dtype.names or ()):
        raise ValueError('segment dataset lacks dx, dQ, nhits or t')
    return np.logical_and.reduce([np.isfinite(rows[f]) for f in FIELDS] +
                                 [rows[f] > 0 for f in FIELDS[:3]])


def _resolve(value, receipt, path_maps):
    original = str(value)
    mappings = []
    for item in path_maps:
        if isinstance(item, str):
            old, sep, new = item.partition('=')
            if not sep or not old or not new:
                raise ValueError('--path-map must be OLD_PREFIX=NEW_PREFIX')
        else:
            old, new = item
        mappings.append((str(old).rstrip('/'), str(new).rstrip('/')))
    for old, new in sorted(mappings, key=lambda pair: len(pair[0]), reverse=True):
        if original == old or original.startswith(old + '/'):
            original = new + original[len(old):]
            break
    path = Path(original)
    if not path.is_absolute():
        path = receipt.parent / path
    if not path.is_file():
        raise FileNotFoundError(f'Missing referenced file: {path}. Summary JSON alone does not '
            'contain track geometry; export on NERSC with its FLOW/sidecar/pool files, '
            'or use --path-map OLD_PREFIX=NEW_PREFIX for relocated files.')
    return path.resolve()


def _layouts(source):
    split = ('analysis/beam_rock_muon_segments', 'analysis/cosmic_muon_segments')
    if any(p + '/data' in source for p in split):
        if not all(p + '/data' in source for p in split):
            raise ValueError('incomplete original split segment selection')
        return list(split)
    if POOL + '/data' in source:
        return [POOL]
    raise ValueError('original FLOW has no stored segment selection')


def _member_inventory(member, index, mode, receipt, path_maps):
    source_path = _resolve(member['input_file'], receipt, path_maps)
    selected_path = _resolve(member['selection_file'], receipt, path_maps)
    with ExitStack() as stack:
        source = stack.enter_context(h5py.File(source_path, 'r'))
        selected = source if source_path == selected_path else stack.enter_context(h5py.File(selected_path, 'r'))
        if not bool(source.attrs.get('complete', True)) or not bool(selected.attrs.get('complete', True)):
            raise ValueError('incomplete FLOW or selection sidecar')
        if EVENTS + '/data' not in source or HITS + '/data' not in source:
            raise ValueError('member.input_file must be the full original FLOW with events and hits')
        events = source[EVENTS + '/data'][:]
        layout_names = (_layouts(source) if mode == 'quality-stored' else
                        [str(Path(d['path']).parent) for d in member['datasets']])
        layouts = {}
        for name in layout_names:
            track_name = name.replace('_segments', '_tracks')
            segments = source[name + '/data'][:]
            tracks = source[track_name + '/data'][:]
            parents = _parents(_reference(source, track_name, name)[:], len(segments), len(tracks), 'segment')
            track_events = _parents(_reference(source, EVENTS, track_name)[:], len(tracks), len(events), 'track')
            if 'event_id' in (tracks.dtype.names or ()) and not np.array_equal(tracks['event_id'], events['id'][track_events]):
                raise ValueError('track event_id disagrees with original event references')
            valid = _valid(segments)
            status = np.where(valid, 'quality_rejected' if mode == 'quality-stored' else 'not_selected', 'invalid').astype('U20')
            layouts[name] = dict(segments=segments, tracks=tracks, parents=parents,
                track_events=track_events, status=status, selected=np.zeros(len(segments), dtype=bool), track_name=track_name)
        pool_pieces, pool_rows = [], []
        mapped = set()
        for dataset_index, record in enumerate(member['datasets']):
            key = record['path']
            rows = selected[key][:]
            if len(rows) != record['n_segments'] or hashlib.sha256(rows.tobytes()).hexdigest() != record['sha256']:
                raise ValueError(f'selection dataset hash/count disagrees with receipt: {selected_path}:{key}')
            group = selected[str(Path(key).parent)]
            if mode == 'quality-stored':
                if 'source_row' not in group or 'source_dataset' not in group:
                    raise ValueError('quality sidecar lacks original segment row lineage')
                originals = group['source_row'][:]
                names = [_text(v).removesuffix('/data').lstrip('/') for v in group['source_dataset'][:]]
            else:
                originals = np.arange(len(rows))
                names = [str(Path(key).parent)] * len(rows)
            if len(originals) != len(rows) or len(names) != len(rows) or originals.dtype.kind not in 'iu':
                raise ValueError('invalid original segment lineage shape/type')
            for i, (name, original_row) in enumerate(zip(names, originals)):
                original_row = int(original_row)
                if name not in layouts or original_row < 0 or original_row >= len(layouts[name]['segments']):
                    raise ValueError('invalid original segment lineage location')
                original = layouts[name]['segments'][original_row]
                if rows.dtype != original.dtype or rows[i].tobytes() != original.tobytes():
                    raise ValueError('selected segment differs from its original FLOW row')
                identity = (name, original_row)
                if identity in mapped:
                    raise ValueError('duplicate original segment lineage')
                mapped.add(identity)
                layouts[name]['selected'][original_row] = True
            valid = _valid(rows)
            for i in np.flatnonzero(valid):
                layouts[names[i]]['status'][int(originals[i])] = 'retained'
            pool_pieces.append(rows[valid])
            pool_rows.append((dataset_index, np.flatnonzero(valid)))
        bounds = np.asarray(source['geometry_info'].attrs['module_RO_bounds'], dtype=float)
        if bounds.shape != (4, 2, 3) or not np.isfinite(bounds).all() or np.any(bounds[:, 0] >= bounds[:, 1]):
            raise ValueError('invalid four-module geometry bounds in xyz cm')
        tick_us = float(source['run_info'].attrs['crs_ticks'])
        if not np.isfinite(tick_us) or tick_us <= 0:
            raise ValueError('invalid source drift tick period')
        diagnostics = []
        if 'quality/events/diagnostics' in selected:
            diagnostics = [json.loads(_text(v)) for v in selected['quality/events/diagnostics'][:]]
        geometry = dict(id=f'member-{index}', bounds=bounds.tolist(),
                        drift_axis=_text(source['geometry_info'].attrs.get('drift_direction', 'x')))
    return dict(index=index, source_path=source_path, selected_path=selected_path,
        layouts=layouts, events=events, pool_pieces=pool_pieces, pool_rows=pool_rows,
        diagnostics=diagnostics, geometry=geometry, tick_us=tick_us)


def _verify_pool(summary, inventories, receipt, path_maps):
    if not summary.get('pooled_segments'):
        raise ValueError('receipt has no pooled_segments; there is no completed pooled sample to visualize')
    path = _resolve(summary['pooled_segments'], receipt, path_maps)
    pieces, members, datasets, rows = [], [], [], []
    for inv in inventories:
        for data, (dataset_index, selected_rows) in zip(inv['pool_pieces'], inv['pool_rows']):
            pieces.append(data)
            members.extend([inv['index']] * len(data))
            datasets.extend([dataset_index] * len(data))
            rows.extend(selected_rows)
    expected = np.concatenate(pieces) if pieces else None
    if expected is None or len(expected) != summary['n_valid_segments']:
        raise ValueError('receipt valid-segment count disagrees with selected inputs')
    with h5py.File(path, 'r') as pool:
        if not bool(pool.attrs.get('complete', False)):
            raise ValueError('incomplete pooled segment file')
        actual = pool[POOL + '/data'][:]
        if expected.dtype != actual.dtype or expected.shape != actual.shape or expected.tobytes() != actual.tobytes():
            raise ValueError('pooled segment payload disagrees with selected inputs')
        for field, values in [('member_index', members), ('source_dataset_index', datasets), ('source_row', rows)]:
            if not np.array_equal(pool[f'{POOL}/{field}'][:], values):
                raise ValueError(f'pooled {field} lineage disagrees with selected inputs')
    return str(path)


def _xyz(row, suffix):
    names = [f'{axis}_{suffix}' for axis in 'xyz']
    if not set(names) <= set(row.dtype.names or ()):
        return None
    values = [float(row[name]) for name in names]
    return values if np.isfinite(values).all() else None


def _value(row, name, default=None):
    return row[name].item() if name in (row.dtype.names or ()) else default


def _track_payload(inv, name, track_row, max_hits, seed, source):
    layout = inv['layouts'][name]
    segment_rows = np.flatnonzero(layout['parents'] == track_row)
    track = layout['tracks'][track_row]
    event_row = int(layout['track_events'][track_row])
    event = inv['events'][event_row]
    key = f"{inv['index']}:{layout['track_name']}:{track_row}"
    segments = []
    for row in segment_rows:
        segment = layout['segments'][row]
        valid = bool(_valid(np.array([segment]))[0])
        segments.append(dict(row=int(row), id=_value(segment, 'rock_segment_id', int(row)),
            start=_xyz(segment, 'start'), end=_xyz(segment, 'end'), mid=_xyz(segment, 'mid'),
            status=str(layout['status'][row]), quality_selected=bool(layout['selected'][row]), numerical_valid=valid,
            io_group=_value(segment, 'io_group'), t_us=float(segment['t']) * inv['tick_us'],
            dqdx=float(segment['dQ'] / segment['dx']) if valid else None,
            nhits=_value(segment, 'nhits'), dx=_value(segment, 'dx'), dQ=_value(segment, 'dQ')))
    output_hits = {field: [] for field in ('x', 'y', 'z', 'q', 'io_group', 'segment_row', 'status', 'hit_row',
                                         'io_channel', 'chip_id', 'channel_id', 't_drift_us', 'ts_pps_us')}
    hit_total = sum(int(s['nhits']) for s in segments)
    nonfinite_hits = 0
    if max_hits:
        refs = _reference(source, name, HITS)
        region_key = f'{name}/ref/{HITS}/ref_region'
        # Only the stored segment reference table, never complete event hits.
        if region_key in source:
            pieces = []
            for row in segment_rows:
                region = source[region_key][row]
                start, stop = int(region['start']), int(region['stop'])
                if start < 0 or stop < start or stop > len(refs):
                    raise ValueError('invalid segment-to-hit reference region')
                part = refs[start:stop]
                pieces.append(part[part[:, 0] == row])
            links = np.concatenate(pieces) if pieces else np.empty((0, 2), dtype='i8')
        else:
            links = refs[:]
            links = links[np.isin(links[:, 0], segment_rows)]
        if (np.any(links < 0) or np.any(links[:, 1] >= len(source[HITS + '/data']))
                or len(np.unique(links, axis=0)) != len(links)):
            raise ValueError('invalid or duplicate segment-to-hit references')
        counts = {int(row): int(np.sum(links[:, 0] == row)) for row in segment_rows}
        if any(counts[s['row']] != s['nhits'] for s in segments):
            raise ValueError('segment nhits disagrees with original hit references')
        event_refs = _reference(source, EVENTS, HITS)
        event_region_key = f'{EVENTS}/ref/{HITS}/ref_region'
        if event_region_key not in source or len(source[event_region_key]) != len(inv['events']):
            raise ValueError('missing or incomplete event-to-hit reference regions')
        region = source[event_region_key][event_row]
        start, stop = int(region['start']), int(region['stop'])
        if start < 0 or stop < start or stop > len(event_refs):
            raise ValueError('invalid event-to-hit reference region')
        event_links = event_refs[start:stop]
        event_members = np.sort(event_links[event_links[:, 0] == event_row, 1])
        if (np.any(event_members < 0) or np.any(event_members >= len(source[HITS + '/data']))
                or len(np.unique(event_members)) != len(event_members)):
            raise ValueError('invalid or duplicate event-to-hit references')
        positions = np.searchsorted(event_members, links[:, 1])
        if np.any(positions >= len(event_members)) or not np.array_equal(event_members[positions], links[:, 1]):
            raise ValueError('segment hit lies outside its parent event references')
        hit_total = len(links)
        if len(links) > max_hits:
            entropy = int.from_bytes(hashlib.sha256(f'{seed}:{key}'.encode()).digest()[:8], 'little')
            rng = np.random.default_rng(entropy)
            links = links[np.sort(rng.choice(len(links), max_hits, replace=False))]
        hits = take(source[HITS + '/data'], links[:, 1])
        if not {'x', 'y', 'z'} <= set(hits.dtype.names or ()):
            raise ValueError('original prompt hits lack xyz coordinates')
        if 'io_group' in layout['segments'].dtype.names and ('io_group' not in hits.dtype.names or
                np.any(hits['io_group'] != layout['segments']['io_group'][links[:, 0]])):
            raise ValueError('segment io_group disagrees with its displayed hit references')
        finite = np.isfinite(np.column_stack([hits[a] for a in 'xyz'])).all(axis=1)
        nonfinite_hits = int((~finite).sum())
        hits, links = hits[finite], links[finite]
        for axis in 'xyz':
            output_hits[axis] = hits[axis].tolist()
        output_hits.update(q=hits['Q'].tolist() if 'Q' in hits.dtype.names else [None] * len(hits),
            io_group=hits['io_group'].tolist() if 'io_group' in hits.dtype.names else [None] * len(hits),
            segment_row=links[:, 0].tolist(), status=layout['status'][links[:, 0]].tolist(), hit_row=links[:, 1].tolist())
        for field in ('io_channel', 'chip_id', 'channel_id'):
            output_hits[field] = hits[field].tolist() if field in hits.dtype.names else [None] * len(hits)
        for field in ('t_drift', 'ts_pps'):
            output_hits[field + '_us'] = (hits[field].astype(float) * inv['tick_us']).tolist() if field in hits.dtype.names else [None] * len(hits)
    timestamp = _value(event, 'unix_ts')
    event_time = datetime.fromtimestamp(timestamp, timezone.utc).isoformat() if timestamp is not None else None
    statuses = layout['status'][segment_rows]
    return dict(key=key, file=inv['source_path'].name, source_path=str(inv['source_path']),
        member_index=inv['index'], track_dataset=layout['track_name'], track_row=int(track_row),
        track_id=_value(track, 'rock_muon_id', int(track_row)), event_id=int(event['id']), event_row=event_row,
        sample='beam' if _value(event, 'n_ext_trigs', 0) > 0 else 'cosmic', event_time=event_time,
        geometry_id=inv['geometry']['id'], n_segments=len(segments), n_retained=int(np.sum(statuses == 'retained')),
        n_quality_rejected=int(np.sum(statuses == 'quality_rejected')), n_invalid=int(np.sum(statuses == 'invalid')),
        length_cm=_value(track, 'length'), io_groups=sorted({s['io_group'] for s in segments if s['io_group'] is not None}),
        quality_diagnostics=[d for d in inv['diagnostics'] if d.get('event_row') == event_row and d.get('source_segments') == name],
        segments=segments, hits=output_hits, hit_total=hit_total, hit_displayed=len(output_hits['x']),
        nonfinite_display_hits_omitted=nonfinite_hits)


def build_payload(summary_path, max_tracks=100, seed=0, population='retained', max_hits_per_track=2000, path_maps=()):
    """Uniform deterministic parent sampling after verifying saved pool lineage."""
    if not 1 <= max_tracks <= 1000 or not 0 <= max_hits_per_track <= 100000:
        raise ValueError('max_tracks must be 1..1000 and max_hits_per_track 0..100000')
    if seed < 0 or population not in ('retained', 'rejected', 'all'):
        raise ValueError('seed must be nonnegative; population must be retained, rejected or all')
    receipt = Path(summary_path).resolve()
    summary = json.loads(receipt.read_text())
    mode = summary.get('selection')
    if mode not in ('stored', 'quality-stored'):
        raise ValueError('viewer supports stored and quality-stored receipts; reselect has different track lineage')
    if summary.get('errors'):
        raise ValueError('receipt contains input errors; no complete pooled sample was fitted')
    if not summary.get('members'):
        raise ValueError('receipt contains no input members')
    inventories = []
    for i, member in enumerate(summary['members']):
        inventories.append(_member_inventory(member, i, mode, receipt, path_maps))
        print(f"Verified {i + 1}/{len(summary['members'])}: {Path(member['input_file']).name}", file=sys.stderr, flush=True)
    pooled_path = _verify_pool(summary, inventories, receipt, path_maps)
    candidates, total = [], 0
    for inv in inventories:
        for name, layout in inv['layouts'].items():
            for row in np.unique(layout['parents']):
                status = layout['status'][layout['parents'] == row]
                retained = np.any(status == 'retained')
                rejected = not retained and np.any(status == 'quality_rejected')
                total += 1
                if population == 'all' or (population == 'retained' and retained) or (population == 'rejected' and rejected):
                    candidates.append((inv['index'], name, int(row)))
    rng = np.random.default_rng(seed)
    picks = np.sort(rng.choice(len(candidates), min(max_tracks, len(candidates)), replace=False))
    chosen = [candidates[i] for i in picks]
    tracks = []
    for inv in inventories:
        local = [item for item in chosen if item[0] == inv['index']]
        if not local:
            continue
        with h5py.File(inv['source_path'], 'r') as source:
            for _, name, row in local:
                tracks.append(_track_payload(inv, name, row, max_hits_per_track, seed, source))
    keys = ('selection', 'period_start', 'period_end', 'n_files', 'n_valid_segments', 'fit_status',
            'lifetime_us', 'error_us', 'membership', 'noise_quality')
    metadata = {key: summary.get(key) for key in keys}
    metadata.update(receipt=str(receipt), pooled_segments=pooled_path,
        fit_message=summary.get('fit_message') or summary.get('fit', {}).get('fit_message'),
        fit_diagnostics=summary.get('fit', {}).get('fit_diagnostics'),
        quality=summary.get('reconstruction_signature', {}).get('quality', {}))
    return _clean(dict(schema_version=1, summary=metadata,
        sampling=dict(seed=seed, total_parent_tracks=total, eligible_parent_tracks=len(candidates),
            exported_tracks=len(tracks), population=population, max_hits_per_track=max_hits_per_track,
            method='uniform parent-track sample without replacement; uniform referenced-hit display subsample'),
        geometries=[inv['geometry'] for inv in inventories if any(t['member_index'] == inv['index'] for t in tracks)],
        tracks=tracks))


def write_viewer(payload, output, json_output=None):
    """Embed data and Plotly in one portable HTML document, with no CDN."""
    from plotly.offline import get_plotlyjs
    output = Path(output)
    if output.suffix.lower() != '.html':
        raise ValueError('viewer output must end in .html')
    encoded = json.dumps(_clean(payload), allow_nan=False, separators=(',', ':'))
    safe = encoded.replace('&', '\\u0026').replace('<', '\\u003c').replace('>', '\\u003e').replace('\u2028', '\\u2028').replace('\u2029', '\\u2029')
    template = Path(__file__).with_suffix('.html').read_text()
    # Replace the data last: user strings must never be interpreted as template tokens.
    html = template.replace('__PLOTLY_JS__', get_plotlyjs()).replace('__TRACK_VIEWER_DATA__', safe)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding='utf-8')
    if json_output:
        json_output = Path(json_output)
        json_output.parent.mkdir(parents=True, exist_ok=True)
        json_output.write_text(encoded + '\n', encoding='utf-8')
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary', type=Path, help='recent_tracks summary.json (including failed attenuation fits)')
    parser.add_argument('--output', required=True, type=Path, help='self-contained .html file')
    parser.add_argument('--json-output', type=Path, help='optional viewer payload for loading into another exported viewer')
    parser.add_argument('--max-tracks', type=int, default=100, help='parent tracks exported (default 100; browser displays 15 at once)')
    parser.add_argument('--max-hits-per-track', type=int, default=2000, help='display hit cap per parent; 0 exports segment geometry only')
    parser.add_argument('--population', choices=('retained', 'rejected', 'all'), default='retained', help='retained: >=1 pooled segment; rejected: no retained and >=1 quality-vetoed; all: original parents')
    parser.add_argument('--seed', type=int, default=0, help='deterministic parent/hit display sampling seed')
    parser.add_argument('--path-map', action='append', default=[], metavar='OLD=NEW', help='replace a recorded file-path prefix; repeatable')
    args = parser.parse_args(argv)
    try:
        receipt = args.summary.resolve()
        recorded = json.loads(receipt.read_text())
        protected = {receipt}
        recorded_paths = [recorded.get('pooled_segments')]
        for member in recorded.get('members', []):
            recorded_paths.extend([member.get('input_file'), member.get('selection_file')])
        protected.update(_resolve(path, receipt, args.path_map) for path in recorded_paths if path)
        outputs = [args.output.resolve()] + ([args.json_output.resolve()] if args.json_output else [])
        if len(set(outputs)) != len(outputs) or any(path in protected for path in outputs):
            raise ValueError('output paths must be distinct from every input summary, FLOW, sidecar and pool file')
        if args.json_output and args.json_output.suffix.lower() != '.json':
            raise ValueError('--json-output must end in .json')
        payload = build_payload(args.summary, args.max_tracks, args.seed, args.population, args.max_hits_per_track, args.path_map)
        result = write_viewer(payload, args.output, args.json_output)
    except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
        parser.exit(2, f'Track viewer export failed: {exc}\n')
    print(f"Exported {len(payload['tracks'])} parent tracks: {result}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
