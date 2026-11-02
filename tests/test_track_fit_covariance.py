"""Numerical regression fixtures, never detector observations or published data.

The first two matrices reproduce the pseudoinverse degeneracy found in the
November 2 pooled fit. A positive diagonal alone is not evidence of a usable
uncertainty; covariance validation must also be independent of parameter units.
"""
import json
from pathlib import Path
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'actions/lifetime'))
import lifetime
import lifetime_funcs as LT


DEGENERATE_COVARIANCES = [
    [[0.011143732110037387, 8.345086032371066e-21],
     [8.345086032371066e-21, 6.249294240028265e-39]],
    [[0.011143732221789599, 9.606142980336582e-21],
     [9.606142980336582e-21, 8.280707138514735e-39]],
]


@pytest.mark.parametrize('covariance', DEGENERATE_COVARIANCES)
def test_positive_variances_do_not_hide_rank_one_covariance(covariance):
    check = LT.attenuation_covariance_diagnostics(covariance)
    assert not check['valid']
    assert 'singular' in check['reason']
    assert check['eigenvalues'][0] <= check['rank_tolerance']
    assert np.all(np.diag(covariance) > 0)


@pytest.mark.parametrize('scales', [[1., 1.], [1e-15, 1e15], [1e15, 1e-15]])
def test_validity_does_not_depend_on_parameter_units(scales):
    scales = np.asarray(scales)
    covariance = np.array([[1., .25], [.25, 1.]]) * np.outer(scales, scales)
    check = LT.attenuation_covariance_diagnostics(covariance)
    assert check['valid']
    np.testing.assert_allclose(check['eigenvalues'], [.75, 1.25])


@pytest.mark.parametrize('covariance', [
    [[1., 0.], [0., 0.]], [[1., 0.], [0., -1.]],
    [[1., 2.], [2., 1.]], [[1., .3], [.2, 1.]],
    [[1., 0.], [0., np.inf]], [[1., np.nan], [np.nan, 1.]],
])
def test_invalid_covariance_checks_are_json_safe(covariance):
    check = LT.attenuation_covariance_diagnostics(covariance)
    assert not check['valid']
    json.dumps(check, allow_nan=False)


def mock_slice_and_attenuation(monkeypatch, covariance, tau=2000.):
    """Exercise the real post-fit path without redoing nineteen slice fits."""
    monkeypatch.setattr(LT, 'langau_fit', lambda *a: (10., 1., 1., 1., .1, 10.))
    monkeypatch.setattr(LT, 'curve_fit',
                        lambda *a, **k: (np.array([10., tau]), np.array(covariance)))
    edges = np.linspace(0, 1960, 20)
    return dict(nhits=np.ones(19 * 5), dqdx=np.full(19 * 5, 10.),
                time_drifted=np.repeat((edges[:-1] + edges[1:]) / 2, 5),
                time_bins=edges, dqdx_bins=np.linspace(0, 250, 70),
                nhits_bins=np.linspace(0, 30, 42), wanted_title='synthetic covariance check')


@pytest.mark.parametrize('covariance', DEGENERATE_COVARIANCES + [
    [[1., 0.], [0., np.inf]], [[1., np.nan], [np.nan, 1.]],
])
def test_failed_candidate_retains_unrounded_evidence_without_a_plot(monkeypatch, covariance):
    arguments = mock_slice_and_attenuation(monkeypatch, covariance, tau=2.8777208821862e10)
    figures_before = set(plt.get_fignums())
    with pytest.raises(LT.InvalidAttenuationCovariance) as raised:
        LT.langau_lifetime(**arguments)
    diagnostics = raised.value.fit_diagnostics
    assert diagnostics['candidate_parameters']['lifetime_us'] == 2.8777208821862e10
    assert not diagnostics['attenuation_covariance_check']['valid']
    assert len(diagnostics['fitted_bins']) == 14
    json.dumps(diagnostics, allow_nan=False)
    assert set(plt.get_fignums()) == figures_before


def test_valid_small_error_keeps_legacy_rounding(monkeypatch):
    arguments = mock_slice_and_attenuation(monkeypatch, [[1., 0.], [0., .01]])
    tau, error, figure = LT.langau_lifetime(**arguments)
    try:
        assert tau == 2000. and error == 0
        assert figure.lifetime_diagnostics['candidate_error_us_unrounded'] == .1
        assert figure.lifetime_diagnostics['attenuation_covariance_check']['valid']
    finally:
        plt.close(figure)


def test_fit_samples_keeps_failure_payload_and_null_measurement(tmp_path, monkeypatch):
    """The monitor consumes fit_samples, so assert its publication contract."""
    segment = np.array([(1., 2., 10., 1.)],
                       dtype=[('dx', 'f8'), ('dQ', 'f8'), ('t', 'f8'), ('nhits', 'f8')])
    source = tmp_path / 'segments.h5'
    with h5py.File(source, 'w') as stream:
        stream['analysis/rock_muon_segments/data'] = segment
    diagnostics = {'candidate_parameters': {'amplitude': 10., 'lifetime_us': 2.8e10},
                   'candidate_error_us_unrounded': 8e-20,
                   'attenuation_covariance': DEGENERATE_COVARIANCES[0],
                   'attenuation_covariance_check':
                       LT.attenuation_covariance_diagnostics(DEGENERATE_COVARIANCES[0])}
    def fail(**kwargs):
        raise LT.InvalidAttenuationCovariance('Invalid attenuation covariance', diagnostics)
    monkeypatch.setattr(LT, 'langau_lifetime', fail)
    output = tmp_path / 'fit.png'
    rows = lifetime.fit_samples(source, output, timestamp='2026-11-02T00:00:00Z')
    assert len(rows) == 1
    row = rows[0]
    assert row['fit_status'] == 'failed'
    assert row['fit_quality'] == 'invalid_attenuation_covariance'
    assert row['lifetime_us'] is None and row['error_us'] is None
    assert row['fit_diagnostics'] == diagnostics
    assert not output.exists()
    json.dumps(rows, allow_nan=False)
