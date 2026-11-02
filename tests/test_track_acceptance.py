"""Track policy boundaries and aggregation invariants, using synthetic fits.

These checks validate the requested monitoring policy, not its physical
coverage or the validity of the historical proxy estimator.
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


@pytest.mark.parametrize('candidate,accepted', [(19999.9999, True), (20000., True),
                                               (20000.0001, False), (25000., False)])
def test_ceiling_uses_unrounded_fit_and_keeps_rejection_evidence(tmp_path, monkeypatch, candidate, accepted):
    source = tmp_path / 'input.h5'
    with h5py.File(source, 'w') as handle:
        handle[lifetime.LEGACY_SAMPLE[1]] = np.array([(1., 10., 100., 2.)],
            dtype=[('dx', 'f8'), ('dQ', 'f8'), ('t', 'f8'), ('nhits', 'f8')])
    def fitted(**kwargs):
        figure = plt.figure()
        figure.lifetime_diagnostics = dict(
            candidate_parameters=dict(amplitude=10., lifetime_us=candidate),
            candidate_error_us_unrounded=2.25, chi2_ndf=120.)
        return round(candidate, 3), 2., figure
    monkeypatch.setattr(LT, 'langau_lifetime', fitted)
    output = tmp_path / 'fit.png'
    before = set(plt.get_fignums())
    row = lifetime.fit_samples(source, output, timestamp='2026-11-02T12:00:00Z')[0]
    assert lifetime.valid_lifetime(row) == accepted
    assert output.exists() == accepted
    assert set(plt.get_fignums()) == before
    # Poor model agreement remains visible, without introducing an unrequested
    # chi-square cutoff or scaling its formal uncertainty to force agreement.
    assert row['fit_diagnostics']['chi2_ndf'] == 120.
    if accepted:
        assert row['fit_status'] == 'ok'
    else:
        assert row['fit_status'] == 'rejected'
        assert row['fit_quality'] == 'above_lifetime_limit'
        assert row['lifetime_us'] is None and row['error_us'] is None
        assert row['candidate_lifetime_us'] == candidate
        assert row['candidate_error_us'] == 2.25
    json.dumps(row, allow_nan=False)


def test_old_track_rows_cannot_bypass_ceiling_other_methods_unchanged():
    for method in ('track', 'track_pool'):
        assert not lifetime.valid_lifetime(dict(method=method, fit_status='ok', lifetime_us=20001., error_us=1.))
        assert lifetime.valid_lifetime(dict(method=method, fit_status='ok', lifetime_us=20000., error_us=1.))
    for method in ('prm', 'gas', 'gas_o2', 'packet', 'packet_pool'):
        assert lifetime.valid_lifetime(dict(method=method, fit_status='ok', lifetime_us=20001., error_us=1.))


def track(name, value, error=10., cohort='a', configuration='one', hour=13):
    return dict(input_file=name, timestamp=f'2026-11-02T{hour:02d}:00:00Z',
        lifetime_us=value, error_us=error, method='track', sample='all_mip',
        fit_status='ok', calculation_source='flow_segments', cohort=cohort,
        configuration_id=configuration)


def test_means_exclude_over_ceiling_and_keep_cohort_configuration_boundaries():
    rows = [track('one', 1000., 30.), track('two', 1400., 40.),
            track('over', 25000.), track('different-cohort', 4000., cohort='b'),
            track('different-config', 5000., configuration='two')]
    result = lifetime.aggregate_track_lifetimes(rows)
    assert len(result) == 3
    first = next(r for r in result if r['cohort']=='a' and r['configuration_id']=='one')
    assert first['n_measurements'] == 2 and first['n_rejected'] == 1
    assert first['lifetime_us'] == 1200.
    assert first['propagated_fit_error_us'] == 25.
    assert first['between_file_sem_us'] == pytest.approx(200.)
    assert first['error_us'] == pytest.approx(200.)


def test_raw_mean_centroids_use_only_their_own_compatible_members():
    rows = [track('early', 1000., cohort='a', hour=13),
            track('late', 2000., cohort='b', hour=17)]
    groups = lifetime.overlay_groups(rows, raw_points=True)
    means = next(iter(groups.values()))['means']
    expected = {'a': '2026-11-02T13:00:00+00:00', 'b': '2026-11-02T17:00:00+00:00'}
    assert {r['cohort']: r['timestamp'] for r in means} == expected
