"""Synthetic rolling-window membership, lineage and CLI failure checks."""
import json
import logging
from pathlib import Path
import sys
from types import ModuleType

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'actions/lifetime'))
import recent_tracks as recent
from lifetime import LEGACY_SAMPLE, SPLIT_SAMPLES
from purity_sources import aware_time


DTYPE = [('rock_segment_id', 'i4'), ('dx', 'f8'), ('dQ', 'f8'), ('nhits', 'f8'), ('t', 'f8')]


def fixture(directory, stamp='2026_11_03_01_00_00_CST', rows=None, split=False, ticks=.1):
    path = directory/f'packet-{stamp}.FLOW.hdf5'
    data = np.array(rows if rows is not None else [(1, 1., 20., 4., 200.)], dtype=DTYPE)
    with h5py.File(path, 'w') as source:
        for name in ('geometry_info', 'lar_info', 'charge/calib_prompt_hits'):
            source.create_group(name)
        source.create_dataset('charge/events/data', data=np.zeros(1, dtype=[('id', 'i4')]))
        source.create_group('run_info')
        source['run_info'].attrs.update(crs_ticks=ticks, e_field=.05, is_mc=False,
                                       charge_thresholds='medm', data_packet_type=0)
        for label, key in SPLIT_SAMPLES if split else [LEGACY_SAMPLE]:
            source.create_dataset(key, data=data[:0] if label == 'beam' else data)
    return path


def gather(paths, tmp_path, **kwargs):
    return recent.gather_recent(paths, aware_time('2026-11-03T06:00:00Z'),
        aware_time('2026-11-03T12:00:00Z'), tmp_path, selection='stored', **kwargs)


def test_rolling_window_uses_acquisition_not_mtime_and_is_half_open(tmp_path):
    before = fixture(tmp_path, '2026_11_02_23_59_59_CST')
    first = fixture(tmp_path, '2026_11_03_00_00_00_CST')
    last = fixture(tmp_path, '2026_11_03_05_59_59_CST', rows=[(1, 1, 21, 4, 200)])
    after = fixture(tmp_path, '2026_11_03_06_00_00_CST')
    data, lineage, receipt = gather([before, first, last, after], tmp_path)
    assert receipt['n_files'] == 2
    assert len(receipt['skipped']) == 2 and not receipt['errors']
    np.testing.assert_array_equal(data['dQ'], [20, 21])
    np.testing.assert_array_equal(lineage['member_index'], [0, 1])
    assert recent.filename_time('packet-2026_11_01_01_30_00_CST.FLOW.h5') - recent.filename_time('packet-2026_11_01_01_30_00_CDT.FLOW.h5') == __import__('pandas').Timedelta(hours=1)


def test_missing_empty_invalid_and_repeated_local_ids_preserve_lineage(tmp_path):
    missing = fixture(tmp_path, '2026_11_03_00_00_00_CST')
    with h5py.File(missing, 'r+') as source:
        del source[LEGACY_SAMPLE[1]]
    empty = fixture(tmp_path, '2026_11_03_01_00_00_CST', rows=[])
    first = fixture(tmp_path, '2026_11_03_02_00_00_CST', rows=[(1, 1, 20, 4, 200), (2, 0, 20, 4, 100)])
    second = fixture(tmp_path, '2026_11_03_03_00_00_CST', rows=[(1, 1, 25, 4, 300), (2, 1, 25, 4, np.nan)])
    data, lineage, receipt = gather([missing, empty, first, second], tmp_path)
    assert receipt['n_files'] == 3 and receipt['n_raw_segments'] == 4
    assert receipt['n_valid_segments'] == 2 and receipt['n_invalid_segments'] == 2
    assert receipt['skipped'][0]['reason'] == 'no stored segment dataset'
    np.testing.assert_array_equal(data['rock_segment_id'], [1, 1])
    np.testing.assert_array_equal(lineage['member_index'], [1, 2])
    np.testing.assert_array_equal(lineage['source_row'], [0, 0])


def test_split_categories_enter_union_once_even_when_legacy_also_exists(tmp_path):
    path = fixture(tmp_path, split=True)
    with h5py.File(path, 'r+') as source:
        source.create_dataset(LEGACY_SAMPLE[1], data=source[SPLIT_SAMPLES[1][1]][:])
    data, lineage, receipt = gather([path], tmp_path)
    assert len(data) == 1 and not receipt['errors']
    assert len(receipt['members'][0]['datasets']) == 2
    np.testing.assert_array_equal(lineage['source_dataset_index'], [1])


def test_discovery_deduplicates_overlapping_globs_and_omits_output(tmp_path):
    path = fixture(tmp_path)
    output = tmp_path/'output'; output.mkdir()
    fixture(output)
    assert recent.discover_inputs([tmp_path, path, str(tmp_path/'*.hdf5')], output) == [path.resolve()]


@pytest.mark.parametrize('failure', ['ticks', 'calibration', 'dtype', 'incomplete', 'copy'])
def test_incompatible_or_duplicated_inputs_are_errors(tmp_path, failure):
    first = fixture(tmp_path)
    second = fixture(tmp_path, '2026_11_03_02_00_00_CST', rows=[(1, 1, 21, 4, 200)])
    with h5py.File(second, 'r+') as source:
        if failure == 'ticks': source['run_info'].attrs['crs_ticks'] = 1.
        if failure == 'calibration': source['charge/calib_prompt_hits'].attrs['gain_file'] = 'different'
        if failure == 'incomplete': source.attrs['complete'] = False
        if failure in ('dtype', 'copy'):
            data = source[LEGACY_SAMPLE[1]][:]
            del source[LEGACY_SAMPLE[1]]
            if failure == 'dtype': data = data.astype([(n, 'f4') for n in data.dtype.names])
            else: data['dQ'] = 20.
            source.create_dataset(LEGACY_SAMPLE[1], data=data)
    _, _, receipt = gather([first, second], tmp_path)
    assert len(receipt['errors']) == 1


def test_missing_metadata_requires_explicit_acknowledgement_but_bad_ticks_still_fail(tmp_path):
    path = fixture(tmp_path)
    with h5py.File(path, 'r+') as source:
        del source['geometry_info']
    assert gather([path], tmp_path)[2]['errors']
    assert not gather([path], tmp_path, allow_missing_metadata=True)[2]['errors']
    with h5py.File(path, 'r+') as source:
        source['run_info'].attrs['crs_ticks'] = 1.
    assert gather([path], tmp_path, allow_missing_metadata=True)[2]['errors']


def test_custom_dataset_metadata_is_compared(tmp_path):
    first = fixture(tmp_path)
    second = fixture(tmp_path, '2026_11_03_02_00_00_CST', rows=[(1, 1, 21, 4, 200)])
    for i, path in enumerate((first, second)):
        with h5py.File(path, 'r+') as source:
            source.create_dataset('custom/data', data=source[LEGACY_SAMPLE[1]][:])
            source['custom'].attrs['selection_version'] = i
    assert gather([first, second], tmp_path, dataset='custom/data')[2]['errors']


def test_cli_reuses_joint_fitter_and_reports_stored_quality(monkeypatch, tmp_path):
    path = fixture(tmp_path)
    calls = []
    def fit(pool, plot, timestamp):
        with h5py.File(pool, 'r') as source:
            calls.append(source[LEGACY_SAMPLE[1]][:])
        for _ in range(3):
            logging.getLogger('pylandau').warning('synthetic fit warning')
        return [dict(fit_status='ok', lifetime_us=2000., error_us=100.)]
    monkeypatch.setattr(recent, 'fit_samples', fit)
    output = tmp_path/'results'
    assert recent.main([str(path), '--selection', 'stored', '--end', '2026-11-03T12:00:00Z', '--output', str(output)]) == 0
    receipt = json.loads(next(output.glob('*/summary.json')).read_text())
    assert len(calls) == 1 and len(calls[0]) == 1
    assert receipt['window_hours'] == 6 and receipt['n_files'] == 1
    assert receipt['noise_quality'].startswith('unverified')
    assert receipt['lifetime_us'] == 2000.
    assert receipt['fit_warnings'] == [dict(logger='pylandau', level='WARNING',
                                           message='synthetic fit warning', occurrences=3)]


def test_fit_warning_counts_first_display_and_restoration_after_failure(caplog):
    import warnings
    logger = logging.getLogger('pylandau')
    original_filters = list(logger.filters)
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError, match='fit failed'):
        with recent.fit_warning_receipt() as records:
            logger.warning('repeated %s', 'message')
            logger.warning('repeated message')
            logger.error('distinct message')
            for _ in range(2):
                logging.getLogger('other-fitter').warning('untouched logger')
            with pytest.warns(RuntimeWarning, match='scipy-style warning'):
                warnings.warn('scipy-style warning', RuntimeWarning)
            raise RuntimeError('fit failed')
    assert records == [dict(logger='pylandau', level='WARNING', message='repeated message', occurrences=2),
                       dict(logger='pylandau', level='ERROR', message='distinct message', occurrences=1)]
    assert [r.getMessage() for r in caplog.records].count('repeated message') == 1
    assert [r.getMessage() for r in caplog.records].count('untouched logger') == 2
    assert logger.filters == original_filters
    logger.warning('repeated message')
    assert [r.getMessage() for r in caplog.records].count('repeated message') == 2


@pytest.mark.parametrize('diagnostic_layout', ['dataset', 'group', 'data_only'])
def test_pool_retains_existing_quality_provenance_for_shared_fitter(monkeypatch, tmp_path, diagnostic_layout):
    path = fixture(tmp_path)
    diagnostic_path = {'dataset': 'quality/events', 'group': 'quality/events/diagnostics',
                       'data_only': 'quality/events/data'}[diagnostic_layout]
    with h5py.File(path, 'r+') as source:
        source.attrs.update(quality_config='{"window_us": 50}', quality_code_sha256='quality-version',
                            selection_origin='previous whole-segment quality veto')
        source.create_dataset(diagnostic_path, data=[b'{"rejected_segments": 3}'] * 3)
        if diagnostic_layout == 'group':
            source.create_dataset('quality/events/data', data=[0, 1, 2])
    def fit(pool, plot, timestamp):
        with h5py.File(pool, 'r') as source:
            assert source.attrs['quality_config'] == '{"window_us": 50}'
            assert source.attrs['quality_code_sha256'] == 'quality-version'
            assert source.attrs['selection_origin'] == 'previous whole-segment quality veto'
        return [dict(fit_status='failed', lifetime_us=None, error_us=None)]
    monkeypatch.setattr(recent, 'fit_samples', fit)
    output = tmp_path/'results'
    assert recent.main([str(path), '--selection', 'stored', '--end', '2026-11-03T12:00:00Z', '--output', str(output)]) == 2
    receipt = json.loads(next(output.glob('*/summary.json')).read_text())
    assert receipt['members'][0]['quality_diagnostics'] == dict(file=str(path.resolve()), dataset=diagnostic_path, n_events=3)
    assert receipt['lifetime_us'] is None


def test_no_data_or_bad_input_never_invokes_fitter(monkeypatch, tmp_path):
    path = fixture(tmp_path, rows=[])
    monkeypatch.setattr(recent, 'fit_samples', lambda *a, **k: pytest.fail('fitter must not run'))
    output = tmp_path/'results'
    args = [str(path), '--selection', 'stored', '--end', '2026-11-03T12:00:00Z', '--output', str(output)]
    assert recent.main(args) == 2
    receipt = json.loads(next(output.glob('*/summary.json')).read_text())
    assert receipt['lifetime_us'] is None and receipt['n_valid_segments'] == 0
    with pytest.raises(SystemExit):
        recent.main(args + ['--end', '2026-11-03T12:00:00'])


def test_nonempty_pool_error_is_reported_and_blocks_fit(monkeypatch, tmp_path, capsys):
    good = fixture(tmp_path)
    bad = fixture(tmp_path, '2026_11_03_02_00_00_CST', rows=[(1, 1, 21, 4, 200)])
    with h5py.File(bad, 'r+') as source:
        source.attrs['complete'] = False
    monkeypatch.setattr(recent, 'fit_samples', lambda *a, **k: pytest.fail('partial-input fitter must not run'))
    output = tmp_path/'results'
    assert recent.main([str(good), str(bad), '--selection', 'stored',
                        '--end', '2026-11-03T12:00:00Z', '--output', str(output)]) == 2
    terminal = capsys.readouterr()
    receipt = json.loads(next(output.glob('*/summary.json')).read_text())
    assert receipt['n_files'] == 1 and receipt['n_valid_segments'] == 1
    assert receipt['fit_status'] == 'unavailable' and receipt['lifetime_us'] is None
    assert receipt['errors'] == [dict(path=str(bad.resolve()),
                                    error='ValueError: incomplete selected-segment input')]
    assert str(bad.resolve()) in terminal.err
    assert 'ValueError: incomplete selected-segment input' in terminal.err
    assert receipt['fit_message'] in terminal.err
    assert '1 input error' in terminal.out
    assert not list(output.glob('*/segments.h5'))


def test_failed_selection_sidecar_is_retained_without_reusing_its_path(monkeypatch, tmp_path):
    paths = [fixture(tmp_path, f'2026_11_03_0{hour}_00_00_CST', rows=[(1, 1, 20+hour, 4, 200)])
             for hour in (1, 2, 3)]
    for index, path in enumerate(paths):
        with h5py.File(path, 'r+') as source:
            source.create_dataset('charge/calib_prompt_hits/data', data=[0])
            source[str(Path(LEGACY_SAMPLE[1]).parent)].attrs['selector_version'] = int(index == 1)
    def select(source_path, output_path, quality_config):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_path, 'w') as output, h5py.File(source_path, 'r') as source:
            output.create_dataset(LEGACY_SAMPLE[1], data=source[LEGACY_SAMPLE[1]][:])
            output.attrs.update(complete=True, quality_code_sha256='quality', tick_period_us=.1)
    module = ModuleType('stored_quality')
    module.select_stored_muons = select
    monkeypatch.setitem(sys.modules, 'stored_quality', module)
    data, _, receipt = recent.gather_recent(paths, aware_time('2026-11-03T06:00:00Z'),
        aware_time('2026-11-03T12:00:00Z'), tmp_path)
    assert len(data) == 2 and len(receipt['errors']) == 1
    failure = receipt['errors'][0]
    assert failure['path'] == str(paths[1].resolve())
    assert Path(failure['selection_file']).name == '00001.segments.h5'
    assert [Path(member['selection_file']).name for member in receipt['members']] == [
        '00000.segments.h5', '00002.segments.h5']
    with h5py.File(failure['selection_file'], 'r') as failed_sidecar:
        np.testing.assert_array_equal(failed_sidecar[LEGACY_SAMPLE[1]]['dQ'], [22.])


def test_quality_stored_is_default_and_requires_hits_before_selection(monkeypatch, tmp_path):
    path = fixture(tmp_path)
    monkeypatch.setattr(recent, 'fit_samples', lambda *a, **k: pytest.fail('fitter must not run'))
    output = tmp_path/'results'
    assert recent.main([str(path), '--end', '2026-11-03T12:00:00Z', '--output', str(output)]) == 2
    receipt = json.loads(next(output.glob('*/summary.json')).read_text())
    assert receipt['selection'] == 'quality-stored'
    assert 'full FLOW hits' in receipt['errors'][0]['error']


def test_reselect_passes_quality_configuration_and_preserves_identity(monkeypatch, tmp_path):
    import track_selection
    path = fixture(tmp_path)
    with h5py.File(path, 'r+') as source:
        source.create_dataset('charge/calib_prompt_hits/data', data=[0])
    received = []
    def select(source_path, output_path, quality_config):
        received.append(quality_config)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_path, 'w') as output, h5py.File(source_path, 'r') as source:
            output.create_dataset(LEGACY_SAMPLE[1], data=source[LEGACY_SAMPLE[1]][:])
            output.attrs.update(complete=True, selector_sha256='selector', quality_code_sha256='quality',
                                quality_config=json.dumps(quality_config), adapter_sha256='adapter', tick_period_us=.1)
    monkeypatch.setattr(track_selection, 'select_muons', select)
    data, _, receipt = recent.gather_recent([path], aware_time('2026-11-03T06:00:00Z'),
        aware_time('2026-11-03T12:00:00Z'), tmp_path, selection='reselect', quality_config={'window_us': 50})
    assert len(data) == 1 and not receipt['errors']
    assert received == [{'window_us': 50}]
    assert receipt['reconstruction_signature']['quality']['quality_code_sha256'] == 'quality'


def test_quality_stored_dispatch_and_upstream_settings_survive_vetting(monkeypatch, tmp_path):
    first = fixture(tmp_path)
    second = fixture(tmp_path, '2026_11_03_02_00_00_CST', rows=[(1, 1, 21, 4, 200)])
    for index, path in enumerate((first, second)):
        with h5py.File(path, 'r+') as source:
            source.create_dataset('charge/calib_prompt_hits/data', data=[0])
            source[str(Path(LEGACY_SAMPLE[1]).parent)].attrs['selector_version'] = index
    received = []
    def select(source_path, output_path, quality_config):
        received.append(quality_config)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(output_path, 'w') as output, h5py.File(source_path, 'r') as source:
            output.create_dataset(LEGACY_SAMPLE[1], data=source[LEGACY_SAMPLE[1]][:])
            output.attrs.update(complete=True, quality_code_sha256='quality', tick_period_us=.1)
    module = ModuleType('stored_quality')
    module.select_stored_muons = select
    monkeypatch.setitem(sys.modules, 'stored_quality', module)
    data, _, receipt = recent.gather_recent([first, second], aware_time('2026-11-03T06:00:00Z'),
        aware_time('2026-11-03T12:00:00Z'), tmp_path, quality_config={'window_us': 50})
    assert len(data) == 1 and len(receipt['errors']) == 1
    assert 'incompatible reconstruction or selection' in receipt['errors'][0]['error']
    assert received == [{'window_us': 50}, {'window_us': 50}]
    assert receipt['reconstruction_signature']['upstream_selections']['all_mip']['group']['selector_version'] == 0
