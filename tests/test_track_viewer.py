"""Viewer provenance and display-sampling checks; not physics validation."""
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
from stored_quality import select_stored_muons
from test_stored_quality import fixture_flow, reference
import track_viewer as viewer


SEGMENTS = 'analysis/rock_muon_segments'
TRACKS = 'analysis/rock_muon_tracks'
POOL_KEY = SEGMENTS + '/data'


def make_flow(path, charge_offset=0., invalid=False, split=False, fully_reject=False):
    """Two parents with IDs unlike rows; shuffled refs include a hot segment."""
    fixture_flow(path)
    with h5py.File(path, 'r+') as f:
        old = f[POOL_KEY][:]
        data = np.zeros(len(old), dtype=old.dtype.descr + [
            (f'{axis}_{where}', 'f8') for where in ('start', 'end', 'mid')
            for axis in ('x', 'y', 'z')])
        for name in old.dtype.names:
            data[name] = old[name]
        data['dQ'] += charge_offset
        data['x_start'] = [20., 10., 20.]
        data['x_end'] = [20., 11., 20.]
        data['x_mid'] = (data['x_start'] + data['x_end']) / 2
        if invalid:
            data['dQ'][2] = np.nan
        del f[POOL_KEY]
        f.create_dataset(POOL_KEY, data=data)
        old_tracks = f[TRACKS + '/data'][:]
        tracks = np.zeros(len(old_tracks), dtype=old_tracks.dtype.descr + [
            ('rock_muon_id', 'i8'), ('length', 'f8')])
        tracks['event_id'] = old_tracks['event_id']
        tracks['rock_muon_id'] = [300, 900]
        tracks['length'] = [100., 120.]
        del f[TRACKS + '/data']
        f.create_dataset(TRACKS + '/data', data=tracks)
        if fully_reject:
            hit = f['charge/calib_prompt_hits/data'][12]
            hit['channel_id'] = 0
            f['charge/calib_prompt_hits/data'][12] = hit
        if split:
            for sample, event_row, track_row, selected, memberships in [
                ('beam_rock_muon', 1, 0, [0], [[0, 13], [0, 14]]),
                ('cosmic_muon', 0, 1, [1, 2], [[0, 0], [0, 1], [1, 12]]),
            ]:
                track_group = f'analysis/{sample}_tracks'
                segment_group = f'analysis/{sample}_segments'
                f.create_dataset(track_group + '/data', data=tracks[[track_row]])
                f.create_dataset(segment_group + '/data', data=data[selected])
                reference(f, 'charge/events', track_group, [[event_row, 0]])
                reference(f, track_group, segment_group, [[0, i] for i in range(len(selected))])
                reference(f, segment_group, 'charge/calib_prompt_hits', memberships)
    return path


def make_receipt(tmp_path, selection='quality-stored', n_files=1,
                 invalid=False, fit_status='ok', split=False, fully_reject=False):
    """Construct the exact persisted pooling contract, without running fits."""
    members, pieces, member_rows, dataset_rows, source_rows = [], [], [], [], []
    for index in range(n_files):
        source = make_flow(tmp_path / f'packet-{index}-2026_11_03_0{index + 1}_00_00_CST.FLOW.hdf5',
                           charge_offset=index, invalid=invalid, split=split, fully_reject=fully_reject)
        if selection == 'quality-stored':
            selected = tmp_path / f'{index:05d}.segments.h5'
            select_stored_muons(source, selected)
            keys = ['analysis/beam_rock_muon_segments/data',
                    'analysis/cosmic_muon_segments/data']
        else:
            selected = source
            keys = (['analysis/beam_rock_muon_segments/data',
                     'analysis/cosmic_muon_segments/data'] if split else [POOL_KEY])
        records = []
        with h5py.File(selected, 'r') as f:
            for dataset_index, key in enumerate(keys):
                data = f[key][:]
                records.append(dict(path=key, n_segments=len(data),
                                    sha256=hashlib.sha256(data.tobytes()).hexdigest()))
                pieces.append(data)
                member_rows.append(np.full(len(data), index, dtype='i4'))
                dataset_rows.append(np.full(len(data), dataset_index, dtype='i4'))
                source_rows.append(np.arange(len(data), dtype='i8'))
        members.append(dict(input_file=str(source), source=str(source),
                            selection_file=str(selected), datasets=records,
                            n_raw_segments=sum(record['n_segments'] for record in records),
                            timestamp='2026-11-03T07:00:00+00:00'))
        if selection == 'quality-stored':
            members[-1]['quality_diagnostics'] = dict(
                file=str(selected), dataset='quality/events/diagnostics', n_events=2)
    raw = np.concatenate(pieces)
    valid = np.logical_and.reduce([np.isfinite(raw[name]) for name in ('dx', 'dQ', 't', 'nhits')]
                                  + [raw[name] > 0 for name in ('dx', 'dQ', 'nhits')])
    pooled = tmp_path / 'segments.h5'
    with h5py.File(pooled, 'w') as f:
        f.create_dataset(POOL_KEY, data=raw[valid])
        for key, values in [('member_index', member_rows),
                            ('source_dataset_index', dataset_rows), ('source_row', source_rows)]:
            f.create_dataset(SEGMENTS + '/' + key, data=np.concatenate(values)[valid])
        f.attrs.update(complete=True, pool_members=json.dumps(members), selection=selection,
                       period_start='2026-11-03T06:00:00+00:00',
                       period_end='2026-11-03T12:00:00+00:00')
    receipt = dict(selection=selection, members=members, errors=[], skipped=[],
                   n_files=n_files, n_raw_segments=len(raw), n_valid_segments=int(valid.sum()),
                   n_invalid_segments=int((~valid).sum()), pooled_segments=str(pooled),
                   output_directory=str(tmp_path), fit_status=fit_status,
                   lifetime_us=321.061 if fit_status == 'ok' else None,
                   error_us=6. if fit_status == 'ok' else None,
                   period_start='2026-11-03T06:00:00+00:00',
                   period_end='2026-11-03T12:00:00+00:00')
    path = tmp_path / 'summary.json'
    path.write_text(json.dumps(receipt))
    return path, receipt


def test_quality_lineage_uses_original_rows_and_parent_references(tmp_path):
    path, _ = make_receipt(tmp_path)
    payload = viewer.build_payload(path)
    tracks = {track['track_row']: track for track in payload['tracks']}
    assert len(tracks) == 2
    assert tracks[0]['event_id'] == 7 and tracks[0]['track_id'] == 300
    assert tracks[1]['event_id'] == 42 and tracks[1]['track_id'] == 900
    segments = {segment['row']: segment for segment in tracks[1]['segments']}
    assert set(segments) == {1, 2}  # Cosmic sidecar row zero maps to original row two.
    assert segments[1]['status'] == 'quality_rejected'
    assert segments[2]['status'] == 'retained'
    assert segments[2]['t_us'] == pytest.approx(100.)
    assert segments[2]['dqdx'] == pytest.approx(10.)
    assert set(tracks[1]['hits']['hit_row']) == {0, 1, 12}


def test_repeated_local_ids_remain_distinct_across_files(tmp_path):
    path, _ = make_receipt(tmp_path, n_files=2)
    payload = viewer.build_payload(path)
    assert len(payload['tracks']) == 4
    assert len({track['key'] for track in payload['tracks']}) == 4
    assert sorted(track['track_id'] for track in payload['tracks']) == [300, 300, 900, 900]


@pytest.mark.parametrize('selection', ['quality-stored', 'stored'])
def test_split_source_datasets_use_their_own_rows_without_legacy_duplicates(tmp_path, selection):
    path, _ = make_receipt(tmp_path, selection=selection, split=True)
    payload = viewer.build_payload(path)
    tracks = {track['track_dataset']: track for track in payload['tracks']}
    assert set(tracks) == {'analysis/beam_rock_muon_tracks', 'analysis/cosmic_muon_tracks'}
    cosmic = tracks['analysis/cosmic_muon_tracks']
    assert cosmic['event_id'] == 42 and cosmic['track_id'] == 900 and cosmic['track_row'] == 0
    assert {segment['row'] for segment in cosmic['segments']} == {0, 1}
    assert cosmic['n_retained'] == (1 if selection == 'quality-stored' else 2)
    assert sum(track['n_segments'] for track in tracks.values()) == 3


def test_rejected_population_selects_parents_with_no_retained_segments(tmp_path):
    path, _ = make_receipt(tmp_path, fully_reject=True)
    payload = viewer.build_payload(path, population='rejected')
    assert len(payload['tracks']) == 1
    track = payload['tracks'][0]
    assert track['track_id'] == 900 and track['n_retained'] == 0
    assert track['n_quality_rejected'] == 2


def test_path_map_restores_relocated_pool_sidecars_and_flow(tmp_path):
    original, moved = tmp_path / 'original', tmp_path / 'moved'
    original.mkdir()
    make_receipt(original)
    original.rename(moved)
    payload = viewer.build_payload(moved / 'summary.json', path_maps=[f'{original}={moved}'])
    assert len(payload['tracks']) == 2
    assert all(Path(track['source_path']).parent == moved for track in payload['tracks'])


def test_numerical_invalidity_is_separate_from_quality_veto(tmp_path):
    path, _ = make_receipt(tmp_path, invalid=True)
    retained = viewer.build_payload(path)
    assert [track['track_row'] for track in retained['tracks']] == [0]
    payload = viewer.build_payload(path, population='all')
    track = next(track for track in payload['tracks'] if track['track_row'] == 1)
    assert {segment['row']: segment['status'] for segment in track['segments']} == {
        1: 'quality_rejected', 2: 'invalid'}


@pytest.mark.parametrize('alteration', ['selected_charge', 'original_charge', 'mapping',
                                        'pooled_charge', 'pooled_lineage', 'reported_count'])
def test_changed_provenance_is_rejected(tmp_path, alteration):
    path, receipt = make_receipt(tmp_path)
    member = receipt['members'][0]
    if alteration == 'reported_count':
        receipt['n_valid_segments'] += 1
        path.write_text(json.dumps(receipt))
    elif alteration in ('selected_charge', 'mapping'):
        with h5py.File(member['selection_file'], 'r+') as f:
            group = f['analysis/cosmic_muon_segments']
            if alteration == 'mapping':
                group['source_row'][0] = 1
            else:
                row = group['data'][0]
                row['dQ'] += 1
                group['data'][0] = row
    elif alteration == 'original_charge':
        with h5py.File(member['input_file'], 'r+') as f:
            row = f[POOL_KEY][2]
            row['dQ'] += 1
            f[POOL_KEY][2] = row
    else:
        with h5py.File(receipt['pooled_segments'], 'r+') as f:
            if alteration == 'pooled_lineage':
                f[SEGMENTS + '/source_row'][0] = 999
            else:
                row = f[POOL_KEY][0]
                row['dQ'] += 1
                f[POOL_KEY][0] = row
    with pytest.raises((ValueError, RuntimeError)):
        viewer.build_payload(path)


@pytest.mark.parametrize('target', ['source', 'sidecar', 'pool'])
def test_cli_cannot_overwrite_referenced_inputs(tmp_path, target):
    path, receipt = make_receipt(tmp_path)
    member = receipt['members'][0]
    protected = Path({'source': member['input_file'], 'sidecar': member['selection_file'],
                      'pool': receipt['pooled_segments']}[target])
    before = protected.read_bytes()
    with pytest.raises(SystemExit) as result:
        viewer.main([str(path), '--output', str(tmp_path / 'viewer.html'),
                     '--json-output', str(protected)])
    assert result.value.code == 2
    assert protected.read_bytes() == before
    assert not (tmp_path / 'viewer.html').exists()


def test_failed_stored_fit_still_has_viewable_tracks(tmp_path):
    path, _ = make_receipt(tmp_path, selection='stored', fit_status='failed')
    payload = viewer.build_payload(path)
    assert payload['summary']['fit_status'] == 'failed'
    assert len(payload['tracks']) == 2
    assert all(segment['status'] == 'retained'
               for track in payload['tracks'] for segment in track['segments'])


@pytest.mark.parametrize('problem', ['input_error', 'reselect', 'missing_pool', 'missing_source', 'missing_reference'])
def test_unsupported_or_incomplete_evidence_is_rejected(tmp_path, problem):
    path, receipt = make_receipt(tmp_path)
    if problem == 'input_error':
        receipt['errors'] = [{'path': 'failed.FLOW.hdf5', 'error': 'unreadable'}]
    elif problem == 'reselect':
        receipt['selection'] = 'reselect'
    elif problem == 'missing_pool':
        Path(receipt['pooled_segments']).unlink()
    elif problem == 'missing_source':
        Path(receipt['members'][0]['input_file']).unlink()
    else:
        with h5py.File(receipt['members'][0]['input_file'], 'r+') as f:
            del f[TRACKS + '/ref/' + SEGMENTS + '/ref']
    path.write_text(json.dumps(receipt))
    with pytest.raises((ValueError, RuntimeError, OSError, KeyError)):
        viewer.build_payload(path)


@pytest.mark.parametrize('problem', ['wrong_event_hit', 'wrong_hit_io_group'])
def test_hit_associations_must_match_parent_event_and_segment_io_group(tmp_path, problem):
    path, receipt = make_receipt(tmp_path)
    with h5py.File(receipt['members'][0]['input_file'], 'r+') as f:
        if problem == 'wrong_event_hit':
            refs = f[SEGMENTS + '/ref/charge/calib_prompt_hits/ref']
            rows = refs[:]
            rows[(rows[:, 0] == 2), 1] = 15  # Valid row and count, wrong parent event.
            refs[:] = rows
        else:
            hit = f['charge/calib_prompt_hits/data'][12]
            hit['io_group'] = 3  # Original segment remains in IO group 1.
            f['charge/calib_prompt_hits/data'][12] = hit
    with pytest.raises((ValueError, RuntimeError)):
        viewer.build_payload(path, max_hits_per_track=1 if problem == 'wrong_event_hit' else 2000)


def test_sampling_is_deterministic_and_hit_limit_is_disclosed(tmp_path):
    path, _ = make_receipt(tmp_path, n_files=2)
    one = viewer.build_payload(path, max_tracks=2, seed=17, max_hits_per_track=1)
    two = viewer.build_payload(path, max_tracks=2, seed=17, max_hits_per_track=1)
    assert [track['key'] for track in one['tracks']] == [track['key'] for track in two['tracks']]
    assert len(one['tracks']) == 2
    for track in one['tracks']:
        assert track['hit_total'] > track['hit_displayed'] == 1
        assert len(track['hits']['x']) == 1


def test_html_embedding_cannot_close_script_and_export_is_self_contained(tmp_path):
    path, _ = make_receipt(tmp_path)
    payload = viewer.build_payload(path)
    marker = '</script><script>window.viewer_injected=true</script>'
    payload['summary']['fit_message'] = marker
    output, json_output = tmp_path / 'viewer.html', tmp_path / 'viewer.json'
    viewer.write_viewer(payload, output, json_output=json_output)
    html = output.read_text()
    assert marker not in html
    assert '__TRACK_VIEWER_DATA__' not in html and '__PLOTLY_JS__' not in html
    assert '__GIF_ENCODER_JS__' not in html and '__TRACK_ANALYSIS_JS__' not in html
    assert 'window.viewer_injected=true' in html
    assert json.loads(json_output.read_text())['summary']['fit_message'] == marker
    assert '<script src=' not in html.lower()


def test_charge_and_quality_evidence_do_not_depend_on_display_sample(tmp_path):
    path, _ = make_receipt(tmp_path, n_files=2)
    one = viewer.build_payload(path, max_tracks=1, seed=7, max_hits_per_track=0)
    all_tracks = viewer.build_payload(path, max_tracks=100, population='all')
    assert len(one['tracks']) == 1 and len(all_tracks['tracks']) == 4
    assert one['charge_scale'] == all_tracks['charge_scale']
    assert one['charge_scale']['direct']['n_segments'] == 4
    assert one['quality_evidence'] == all_tracks['quality_evidence']
    assert one['quality_evidence']['n_quality_rejected'] == 2
    skipped = viewer.build_payload(path, fit_charge=False)
    assert skipped['charge_scale']['direct']['status'] == 'unavailable'
    assert skipped['charge_scale']['direct']['intercept'] is None


def test_event_evidence_is_not_weighted_by_parent_multiplicity():
    event = dict(event_row=4, hot_pixel_hits=0, coherent_noise_hits=100,
                 version='archived', ambiguous_clock_hits=None,
                 burst_groups_flagged=dict(chip=2, tile=1, module=0))
    inventories = [dict(index=0, layouts={}, diagnostics=[event, event.copy()]),
                   dict(index=1, layouts={}, diagnostics=[event.copy()])]
    evidence = viewer._quality_evidence(inventories)
    assert evidence['n_events'] == 2  # File-local event rows remain distinct.
    assert evidence['counters']['coherent_noise_hits'] == 200
    assert evidence['burst_groups_flagged'] == dict(chip=4, tile=2, module=0)
    assert evidence['counters']['hot_pixel_hits'] == 0
    assert evidence['counters']['ambiguous_clock_hits'] is None
    assert evidence['counters']['n_unmapped_hits'] is None


def test_mixed_tick_periods_do_not_get_a_single_conversion_in_direct_fit(tmp_path):
    path, receipt = make_receipt(tmp_path, n_files=2)
    with h5py.File(receipt['members'][1]['input_file'], 'r+') as flow:
        flow['run_info'].attrs['crs_ticks'] = .2
    payload = viewer.build_payload(path)
    assert payload['charge_scale']['direct']['status'] == 'unavailable'
    assert 'different drift tick' in payload['charge_scale']['direct']['fit_message']
    assert len(payload['tracks']) == 4
