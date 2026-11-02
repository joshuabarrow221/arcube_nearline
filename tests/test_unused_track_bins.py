"""Regression tests for the unused-slice failure found in November FLOW.

Synthetic peak values isolate the control-flow requirement: changes outside
the established attenuation range cannot affect fitted data, initialization or
the result. The final attenuation solver and covariance checks remain real.
These fixtures are not detector observations or physical model validation.
"""
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
import lifetime_funcs as LT


EDGES = np.linspace(0, 1960, 20)
MIDPOINTS = (EDGES[:-1] + EDGES[1:]) / 2
ACTIVE = list(range(4, 18))
EXCLUDED = [0, 1, 2, 3, 18]


def arguments(excluded_count=0, sparse_bins=()):
    """Make supported fitted slices and unusable excluded histogram entries."""
    hits, charges, times = [], [], []
    for index, midpoint in enumerate(MIDPOINTS):
        if index in ACTIVE:
            count = 4 if index in sparse_bins else 7
            hits.extend([3.] * count)
            charges.extend([20.] * count)
        else:
            count = excluded_count
            # Outside both histogram axes: fitting any excluded slice would
            # fail even if it had enough segment rows to meet the count floor.
            hits.extend([1000.] * count)
            charges.extend([1000.] * count)
        times.extend([midpoint] * count)
    return dict(nhits=np.array(hits), dqdx=np.array(charges),
                time_drifted=np.array(times), time_bins=EDGES,
                dqdx_bins=np.linspace(0, 250, 70),
                nhits_bins=np.linspace(0, 30, 42), wanted_title='Synthetic slice invariance')


def controlled_peaks(monkeypatch, peak_bins=None, failure_bins=(), invalid_error_bins=()):
    """Return known active-slice peaks while recording every peak-fit call."""
    calls = []
    peak_bins = ACTIVE if peak_bins is None else list(peak_bins)
    expected = 12. * np.exp(-MIDPOINTS[peak_bins] / 10. / 1500.)

    def peak_fit(histogram, centers):
        # A call beyond the included slices is itself a failed regression.
        index = len(calls)
        assert index < len(peak_bins)
        calls.append((histogram.copy(), centers.copy()))
        if peak_bins[index] in failure_bins:
            raise ValueError('Synthetic peak evaluation failure')
        error = np.inf if peak_bins[index] in invalid_error_bins else .2
        return expected[index], 1., 1., 1., error, expected[index]

    monkeypatch.setattr(LT, 'langau_fit', peak_fit)
    return calls, expected


@pytest.mark.parametrize('excluded_count', [0, 3, 100])
def test_empty_sparse_or_pathological_unused_slices_do_not_change_fit(monkeypatch, excluded_count):
    calls, expected = controlled_peaks(monkeypatch)
    tau, error, figure = LT.langau_lifetime(**arguments(excluded_count))
    try:
        diagnostics = figure.lifetime_diagnostics
        assert len(calls) == 14
        assert tau == pytest.approx(1500., abs=.001)
        assert error > 0
        assert diagnostics['fitted_bins'] == ACTIVE
        assert diagnostics['nominal_fitted_bins'] == ACTIVE
        assert diagnostics['usable_drift_bin_count'] == 14
        assert diagnostics['minimum_usable_drift_bins'] == 10
        assert diagnostics['excluded_bins'] == EXCLUDED
        assert diagnostics['estimator_version'] == LT.TRACK_ESTIMATOR_VERSION
        assert diagnostics['attenuation_initial_guess'] == [float(expected.max()), 2000]
        np.testing.assert_allclose(np.array(diagnostics['mpv'])[ACTIVE].astype(float), expected)
        for index in EXCLUDED:
            assert diagnostics['mpv'][index] is None
            assert diagnostics['mpv_error'][index] is None
            assert diagnostics['segments_per_bin'][index] == excluded_count
        np.testing.assert_allclose(diagnostics['drift_time_us'], MIDPOINTS / 10.)
        json.dumps(diagnostics, allow_nan=False)
        # The attenuation curve is confined to the actually fitted drift range.
        model_line = next(line for line in figure.axes[0].lines
                          if len(line.get_xdata()) == 1000)
        np.testing.assert_allclose(model_line.get_xdata()[[0, -1]], MIDPOINTS[[4, 17]] / 10.)
    finally:
        plt.close(figure)


@pytest.mark.parametrize('sparse_bin', [4, 10, 17])
def test_one_sparse_nominal_slice_is_recorded_without_changing_other_peaks(monkeypatch, sparse_bin):
    used = [i for i in ACTIVE if i != sparse_bin]
    controlled_peaks(monkeypatch, peak_bins=used)
    tau, error, figure = LT.langau_lifetime(**arguments(excluded_count=100, sparse_bins=[sparse_bin]))
    try:
        d = figure.lifetime_diagnostics
        assert tau == pytest.approx(1500., abs=.001)
        assert d['fitted_bins'] == used
        assert d['ndf'] == 11
        assert d['mpv'][sparse_bin] is None
        assert d['skipped_bins'][0]['bin'] == sparse_bin
        assert d['skipped_bins'][0]['status'] == 'insufficient_histogram_support'
    finally:
        plt.close(figure)


def test_exactly_ten_usable_bins_pass_with_eight_degrees_of_freedom(monkeypatch):
    used, skipped = ACTIVE[:10], ACTIVE[10:]
    controlled_peaks(monkeypatch, peak_bins=used)
    tau, error, figure = LT.langau_lifetime(**arguments(sparse_bins=skipped))
    try:
        d = figure.lifetime_diagnostics
        assert tau == pytest.approx(1500., abs=.001)
        assert error > 0
        assert d['fitted_bins'] == used
        assert d['usable_drift_bin_count'] == 10
        assert d['ndf'] == 8
        assert [r['bin'] for r in d['skipped_bins']] == skipped
        json.dumps(d, allow_nan=False)
    finally:
        plt.close(figure)


def test_nine_usable_bins_fail_before_attenuation_and_retain_reasons(monkeypatch):
    used, skipped = ACTIVE[:9], ACTIVE[9:]
    controlled_peaks(monkeypatch, peak_bins=used)
    def forbidden(*args, **kwargs):
        pytest.fail('Attenuation solver must not run with only nine usable bins')
    monkeypatch.setattr(LT, 'curve_fit', forbidden)
    figures_before = set(plt.get_fignums())
    with pytest.raises(LT.TrackFitError, match='9; require at least 10') as raised:
        LT.langau_lifetime(**arguments(sparse_bins=skipped))
    assert raised.value.fit_quality == 'insufficient_usable_drift_bins'
    d = raised.value.fit_diagnostics
    assert d['fitted_bins'] == used
    assert d['nominal_fitted_bins'] == ACTIVE
    assert len(d['skipped_bins']) == 5
    assert d['usable_drift_bin_count'] == 9
    json.dumps(d, allow_nan=False)
    assert set(plt.get_fignums()) == figures_before


def test_peak_failures_and_invalid_errors_are_selected_before_attenuation(monkeypatch):
    failed, invalid = [4, 8], [10, 17]
    calls, _ = controlled_peaks(monkeypatch, failure_bins=failed, invalid_error_bins=invalid)
    used = [i for i in ACTIVE if i not in failed + invalid]
    tau, error, figure = LT.langau_lifetime(**arguments())
    try:
        d = figure.lifetime_diagnostics
        assert len(calls) == 14
        assert d['fitted_bins'] == used
        assert d['usable_drift_bin_count'] == 10
        assert d['ndf'] == 8
        assert tau == pytest.approx(1500., abs=.001)
        expected_amplitude = float(12. * np.exp(-MIDPOINTS[used[0]] / 10. / 1500.))
        assert d['attenuation_initial_guess'] == [expected_amplitude, 2000]
        statuses = {r['bin']: r['status'] for r in d['skipped_bins']}
        assert statuses == {4: 'peak_fit_failed', 8: 'peak_fit_failed',
                            10: 'invalid_peak_or_uncertainty', 17: 'invalid_peak_or_uncertainty'}
        for i in failed + invalid:
            assert d['mpv'][i] is None and d['mpv_error'][i] is None
        json.dumps(d, allow_nan=False)
    finally:
        plt.close(figure)


def test_peak_runtime_error_preserves_actual_cause(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError('Cannot calculate MPV, check function parameters')
    monkeypatch.setattr(LT, 'curve_fit', fail)
    with pytest.raises(ValueError, match='Cannot calculate MPV') as raised:
        LT.langau_fit(np.ones(8), np.arange(1., 9.))
    assert isinstance(raised.value.__cause__, RuntimeError)
    assert 'did not converge' not in str(raised.value)


def test_invalid_included_peak_error_is_not_ignored(monkeypatch):
    monkeypatch.setattr(LT, 'langau_fit', lambda *args: (10., 1., 1., 1., np.inf, 10.))
    with pytest.raises(LT.TrackFitError, match='0; require at least 10') as raised:
        LT.langau_lifetime(**arguments())
    d = raised.value.fit_diagnostics
    assert len(d['skipped_bins']) == 14
    assert all(r['status'] == 'invalid_peak_or_uncertainty' for r in d['skipped_bins'])
    json.dumps(d, allow_nan=False)
