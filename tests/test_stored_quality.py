"""Synthetic association/charge-preservation checks, not detector validation."""
from pathlib import Path
import hashlib
import json
import sys
import weakref

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
from stored_quality import QualityContext, build_quality_context, select_stored_muons


HIT_DTYPE = [(name, 'f8') for name in ('x', 'y', 'z', 'Q', 'E', 'ts_pps', 't_drift')]
HIT_DTYPE += [(name, 'i8') for name in ('id', 'io_group', 'io_channel', 'chip_id', 'channel_id')]


def reference(source, first, second, rows):
    node = source.create_dataset(f'{first}/ref/{second}/ref', data=np.asarray(rows, dtype='i8').reshape(-1, 2))
    node.attrs.update(dset0='/' + first + '/data', dset1='/' + second + '/data')
    return node


def fixture_flow(path):
    hits = np.zeros(16, dtype=HIT_DTYPE)
    hits['id'] = np.arange(16) + 100  # references are rows, not these IDs
    hits['io_group'] = hits['io_channel'] = hits['chip_id'] = 1
    hits['Q'] = hits['E'] = 1
    hits['x'] = 20
    hits['x'][:12] = np.linspace(10, 20, 12)
    hits['ts_pps'] = np.arange(16) * 10
    hits['channel_id'][12:] = [1, 2, 3, 4]
    segments = np.array([(10, 2., 40., 2, 1000., 1), (11, 2., 80., 2, 1000., 1), (12, 2., 20., 1, 1000., 1)],
        dtype=[('rock_segment_id', 'i4'), ('dx', 'f8'), ('dQ', 'f8'), ('nhits', 'i4'), ('t', 'f8'), ('io_group', 'i4')])
    with h5py.File(path, 'w') as source:
        source.create_dataset('charge/events/data', data=np.array([(42, 0, 13), (7, 1, 3)],
            dtype=[('id', 'i8'), ('n_ext_trigs', 'i8'), ('nhit', 'i8')]))
        source.create_dataset('charge/calib_prompt_hits/data', data=hits)
        reference(source, 'charge/events', 'charge/calib_prompt_hits', [[0, i] for i in range(13)] + [[1, i] for i in range(13, 16)])
        source.create_dataset('charge/events/ref/charge/calib_prompt_hits/ref_region',
            data=np.array([(0, 13), (13, 16)], dtype=[('start', 'i8'), ('stop', 'i8')]))
        source.create_dataset('analysis/rock_muon_tracks/data', data=np.array([(7,), (42,)], dtype=[('event_id', 'i8')]))
        source.create_dataset('analysis/rock_muon_segments/data', data=segments)
        reference(source, 'charge/events', 'analysis/rock_muon_tracks', [[1, 0], [0, 1]])
        reference(source, 'analysis/rock_muon_tracks', 'analysis/rock_muon_segments', [[1, 2], [0, 0], [1, 1]])
        reference(source, 'analysis/rock_muon_segments', 'charge/calib_prompt_hits', [[2, 12], [1, 1], [0, 13], [1, 0], [0, 14]])
        source.create_group('run_info').attrs.update(crs_ticks=.1, rollover_ticks=10000000)
        geometry = source.create_group('geometry_info')
        geometry.attrs.update(drift_direction='x', module_RO_bounds=np.tile([[0., -70., -70.], [70., 70., 70.]], (4, 1, 1)))
        tile = source.create_group('geometry_info/tile_id')
        meta = np.zeros(1, dtype=[('min_max_keys', 'i8', (2, 2)), ('default', 'i4')])
        meta['min_max_keys'] = [[1, 8], [1, 4]]
        meta['default'] = -1
        tile.attrs['meta'] = meta
        data = np.zeros(33, dtype=[('data', 'i4'), ('filled', '?')])
        data['data'][0] = -1
        data['data'][1:] = 1
        data['filled'][1:] = True
        tile.create_dataset('data', data=data)
    return segments


def test_noise_veto_preserves_charge_original_rows_and_event_lineage(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    before = hashlib.sha256(source.read_bytes()).hexdigest()
    select_stored_muons(source, output)
    assert hashlib.sha256(source.read_bytes()).hexdigest() == before
    with h5py.File(output, 'r') as result:
        np.testing.assert_array_equal(result['analysis/beam_rock_muon_segments/data'][:], original[[0]])
        np.testing.assert_array_equal(result['analysis/cosmic_muon_segments/data'][:], original[[2]])
        np.testing.assert_array_equal(result['analysis/cosmic_muon_segments/source_row'][:], [2])
        assert result.attrs['complete']
        assert result.attrs['events_processed'] == 2
        assert json.loads(result.attrs['quality_config'])['tick_us'] == .1
        audit = result['quality/events/data'][:]
        noisy = audit[audit['event_id'] == 42][0]
        assert noisy['n_input_hits'] == 13
        assert noisy['n_kept_hits'] == 1
        assert noisy['n_segments'] == 2 and noisy['n_kept_segments'] == 1


@pytest.mark.parametrize('failure', ['missing', 'orientation', 'duplicate', 'membership', 'wrong_event', 'incomplete', 'wrong_iog'])
def test_incomplete_or_inconsistent_references_fail_closed(tmp_path, failure):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        key = 'analysis/rock_muon_segments/ref/charge/calib_prompt_hits/ref'
        if failure == 'missing':
            del h[key]
        elif failure == 'orientation':
            h[key].attrs['dset0'] = '/charge/calib_prompt_hits/data'
        elif failure == 'duplicate':
            h[key][1] = h[key][0]
        elif failure == 'membership':
            h[key][0] = [2, 15]  # hit belongs to another event
        elif failure == 'wrong_event':
            h['analysis/rock_muon_tracks/data'][0] = (42,)
        elif failure == 'incomplete':
            rows = h['analysis/rock_muon_segments/data'][:]
            rows['nhits'][0] = 3
            h['analysis/rock_muon_segments/data'][:] = rows
        else:
            rows = h['analysis/rock_muon_segments/data'][:]
            rows['io_group'][0] = 2
            h['analysis/rock_muon_segments/data'][:] = rows
    with pytest.raises(ValueError):
        select_stored_muons(source, output)
    assert not output.exists()


def test_bounds_reject_whole_segment_and_report_separately(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        rows = h['charge/calib_prompt_hits/data'][:]
        rows['x'][13] = 100.
        h['charge/calib_prompt_hits/data'][:] = rows
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as result:
        assert len(result['analysis/beam_rock_muon_segments/data']) == 0
        assert result['quality/events/data']['n_outside_module_hits'].sum() == 1


def test_ambiguous_clock_does_not_abort_file_or_reject_other_groups(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        rows = h['charge/calib_prompt_hits/data'][:]
        rows['ts_pps'][:12] = np.linspace(0, 9_900_000, 12)
        rows['io_group'][12] = 3  # Clean independent clock in the same event.
        h['charge/calib_prompt_hits/data'][:] = rows
        original['io_group'][2] = 3
        h['analysis/rock_muon_segments/data'][:] = original
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as result:
        assert result.attrs['complete']
        np.testing.assert_array_equal(result['analysis/beam_rock_muon_segments/data'][:], original[[0]])
        np.testing.assert_array_equal(result['analysis/cosmic_muon_segments/data'][:], original[[2]])
        audits = [json.loads(row) for row in result['quality/events/diagnostics'][:]]
        noisy = next(row for row in audits if row['event_id'] == 42)
        assert noisy['ambiguous_clock_hits'] == 12
        assert noisy['ambiguous_clock_groups'] == 1
        assert noisy['n_kept_hits'] == 1


def test_source_clock_conflict_and_source_overwrite_fail(tmp_path):
    source = tmp_path / 'source.h5'
    fixture_flow(source)
    with h5py.File(source, 'r') as h, pytest.raises(ValueError, match='conflicts'):
        build_quality_context(h, {'tick_us': 1.})
    with pytest.raises(ValueError, match='replace source'):
        select_stored_muons(source, source)


def test_outside_hits_still_supply_noise_evidence_for_inside_segment(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        rows = h['charge/calib_prompt_hits/data'][:]
        rows['x'][:12] = np.linspace(70., 90., 12)
        h['charge/calib_prompt_hits/data'][:] = rows
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as h:
        # Segment 11 uses the two in-bounds hot hits. The other ten hits still
        # establish repeated activity despite their unphysical drift positions.
        assert list(h['analysis/cosmic_muon_segments/data']['rock_segment_id']) == [12]
        diagnostics = [json.loads(row) for row in h['quality/events/diagnostics'].asstr()[:]]
        noisy = next(row for row in diagnostics if row['event_id'] == 42)
        assert noisy['hot_pixel_hits'] == 12
        assert noisy['n_outside_module_hits'] == 10


def test_empty_stored_selection_is_valid_and_complete(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        for path in ['analysis/rock_muon_tracks/data', 'analysis/rock_muon_segments/data',
                     'charge/events/ref/analysis/rock_muon_tracks/ref',
                     'analysis/rock_muon_tracks/ref/analysis/rock_muon_segments/ref',
                     'analysis/rock_muon_segments/ref/charge/calib_prompt_hits/ref']:
            array, attrs = h[path][:0], dict(h[path].attrs)
            del h[path]
            h.create_dataset(path, data=array).attrs.update(attrs)
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as h:
        assert len(h['analysis/cosmic_muon_segments/data']) == 0
        assert len(h['quality/events/data']) == 0
        assert h.attrs['complete']


def test_explicit_split_selection_retains_original_rows_and_categories(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        for sample, event_row, event_id, selected, memberships in [
            ('beam_rock_muon', 1, 7, [0], [[0, 13], [0, 14]]),
            ('cosmic_muon', 0, 42, [1, 2], [[0, 0], [0, 1], [1, 12]]),
        ]:
            track, segment = f'analysis/{sample}_tracks', f'analysis/{sample}_segments'
            h.create_dataset(track + '/data', data=np.array([(event_id,)], dtype=[('event_id', 'i8')]))
            h.create_dataset(segment + '/data', data=original[selected])
            reference(h, 'charge/events', track, [[event_row, 0]])
            reference(h, track, segment, [[0, i] for i in range(len(selected))])
            reference(h, segment, 'charge/calib_prompt_hits', memberships)
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as h:
        np.testing.assert_array_equal(h['analysis/beam_rock_muon_segments/data'][:], original[[0]])
        np.testing.assert_array_equal(h['analysis/cosmic_muon_segments/data'][:], original[[2]])
        assert list(h['analysis/cosmic_muon_segments/source_row'][:]) == [1]
        assert h['analysis/cosmic_muon_segments/source_dataset'].asstr()[0] == 'analysis/cosmic_muon_segments'


@pytest.mark.parametrize('split', [False, True])
def test_event_hit_buffers_are_released_before_reading_next_event(tmp_path, monkeypatch, split):
    """A field view must not retain every previous event's structured hit array."""
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    if split:
        with h5py.File(source, 'r+') as h:
            for sample, event_row, event_id, selected, memberships in [
                ('beam_rock_muon', 1, 7, [0], [[0, 13], [0, 14]]),
                ('cosmic_muon', 0, 42, [1, 2], [[0, 0], [0, 1], [1, 12]]),
            ]:
                track, segment = f'analysis/{sample}_tracks', f'analysis/{sample}_segments'
                h.create_dataset(track + '/data', data=np.array([(event_id,)], dtype=[('event_id', 'i8')]))
                h.create_dataset(segment + '/data', data=original[selected])
                reference(h, 'charge/events', track, [[event_row, 0]])
                reference(h, track, segment, [[0, i] for i in range(len(selected))])
                reference(h, segment, 'charge/calib_prompt_hits', memberships)

    evaluate = QualityContext.evaluate
    previous_arrays, event_sizes = [], []

    def tracked_evaluate(context, hits):
        # Reachability directly tests retention, independently of allocator/RSS
        # behavior and without depending on the event-cache implementation.
        assert all(reference() is None for reference in previous_arrays)
        previous_arrays.append(weakref.ref(hits))
        event_sizes.append(len(hits))
        return evaluate(context, hits)

    monkeypatch.setattr(QualityContext, 'evaluate', tracked_evaluate)
    select_stored_muons(source, output)
    assert sorted(event_sizes) == [3, 13]
    assert all(reference() is None for reference in previous_arrays)
    with h5py.File(output, 'r') as result:
        np.testing.assert_array_equal(result['analysis/beam_rock_muon_segments/data'][:], original[[0]])
        np.testing.assert_array_equal(result['analysis/cosmic_muon_segments/data'][:], original[[2]])


def test_multiple_parents_share_one_full_event_evaluation(tmp_path, monkeypatch):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    with h5py.File(source, 'r+') as h:
        tracks = 'analysis/rock_muon_tracks'
        segments = 'analysis/rock_muon_segments'
        del h[tracks + '/data']
        # Event 42 now has two parents; event 7 also has an empty parent.
        h.create_dataset(tracks + '/data', data=np.array([(7,), (42,), (42,), (7,)],
            dtype=[('event_id', 'i8')]))
        del h[f'charge/events/ref/{tracks}/ref']
        del h[f'{tracks}/ref/{segments}/ref']
        reference(h, 'charge/events', tracks, [[0, 2], [1, 0], [1, 3], [0, 1]])
        reference(h, tracks, segments, [[2, 2], [0, 0], [1, 1]])
    evaluate = QualityContext.evaluate
    sizes = []

    def counted_evaluate(context, hits):
        sizes.append(len(hits))
        return evaluate(context, hits)

    monkeypatch.setattr(QualityContext, 'evaluate', counted_evaluate)
    select_stored_muons(source, output)
    assert sorted(sizes) == [3, 13]
    with h5py.File(output, 'r') as result:
        np.testing.assert_array_equal(result['analysis/cosmic_muon_segments/data'][:], original[[2]])
        audits = [json.loads(row) for row in result['quality/events/diagnostics'].asstr()[:]]
        assert len(audits) == 2
        noisy = next(row for row in audits if row['event_id'] == 42)
        # Both parent tracks see the other ten unselected hot-pixel hits.
        assert noisy['hot_pixel_hits'] == 12
        assert noisy['n_segments'] == 2 and noisy['n_kept_segments'] == 1
        assert result['analysis/cosmic_muon_segments'].attrs['selected_track_count'] == 1


def test_scattered_event_hit_rows_preserve_quality_and_segment_values(tmp_path):
    source, output = tmp_path / 'source.h5', tmp_path / 'selected.h5'
    original = fixture_flow(source)
    permutation = np.array([0, 13, 1, 14, 2, 15, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    inverse = np.argsort(permutation)
    with h5py.File(source, 'r+') as h:
        hits = h['charge/calib_prompt_hits/data'][:]
        h['charge/calib_prompt_hits/data'][:] = hits[permutation]
        for parent in ('charge/events', 'analysis/rock_muon_segments'):
            node = h[parent + '/ref/charge/calib_prompt_hits/ref']
            refs = node[:]
            refs[:, 1] = inverse[refs[:, 1]]
            node[:] = refs
    select_stored_muons(source, output)
    with h5py.File(output, 'r') as h:
        np.testing.assert_array_equal(h['analysis/beam_rock_muon_segments/data'][:], original[[0]])
        np.testing.assert_array_equal(h['analysis/cosmic_muon_segments/data'][:], original[[2]])
        audits = [json.loads(row) for row in h['quality/events/diagnostics'].asstr()[:]]
        assert next(row for row in audits if row['event_id'] == 42)['hot_pixel_hits'] == 12
