"""Unit tests for bpscan. Run:  python3 -m pytest -q tests/"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import numi, stats, target_model as tm, ifbeam, detector_rates as dr  # noqa: E402


def test_spill_structure():
    assert abs(numi.BUCKET_NS - 18.83) < 0.01
    assert abs(numi.BATCH_NS - 1581.8) < 0.5
    assert abs(numi.SPILL_NS - 9491) < 2


def test_batch_index_and_guards():
    e = numi.nominal_batch_edges_ns(100.0)
    t = np.array([50, 100, 120, e[1] - 1, e[1] + 1, e[6] - 1, e[6] + 5])
    assert list(numi.batch_index(t, e)) == [0, 1, 1, 1, 2, 6, 0]
    # with a 30 ns guard, times within 30 ns of any edge are rejected
    assert list(numi.batch_index(t, e, guard_ns=30)) == [0, 0, 0, 0, 0, 0, 0]
    assert numi.batch_index([e[2] + 500], e, 30)[0] == 3


def test_extrapolation_matches_cpp_formula():
    # a straight line through (z121, 1) and (ztgt, 2) evaluated at z=0
    x = numi.extrapolate(1.0, numi.Z_HP121, 2.0, numi.Z_HPTGT)
    slope = (2 - 1) / (numi.Z_HPTGT - numi.Z_HP121)
    assert np.isclose(x, 2.0 + slope * (0 - numi.Z_HPTGT))


def test_intercept_model():
    assert np.isclose(tm.intercept_fraction(0, 1.0, 1e6), 1.0)
    assert np.isclose(tm.intercept_fraction(0, 1.0, 2.0), 0.6827, atol=1e-4)
    # symmetric offsets about the centre give zero difference
    assert abs(tm.two_position_difference(-0.15, 0.15, 1.43, 7.4)) < 1e-12
    # curvature matches a finite difference
    k = tm.curvature(1.43, 7.4)
    h = 0.05
    f = lambda x: tm.intercept_fraction(x, 1.43, 7.4) / tm.intercept_fraction(0, 1.43, 7.4)
    assert np.isclose(1 - f(h), k * h * h, rtol=1e-2)


def test_required_counts_formula_against_toys():
    rng = np.random.default_rng(0)
    n_tot, phi = 60000, 1 / 6
    r = []
    for _ in range(4000):
        n1 = rng.poisson(n_tot * phi); n2 = rng.poisson(n_tot * (1 - phi))
        r.append(stats.ratio_and_error(n1, phi, n2, 1 - phi)[0])
    expected = np.sqrt(1 / (n_tot * phi) + 1 / (n_tot * (1 - phi)))
    assert np.isclose(np.std(r), expected, rtol=0.05)
    assert np.isclose(stats.required_total_counts(0.01, 3), 7.2 * 9 / 1e-4)


def test_exact_ratio_test_null():
    p, (lo, hi) = stats.exact_ratio_test(100, 1.0, 500, 5.0)
    assert p > 0.5 and lo < 1 < hi


def test_template_scan_unbiased():
    R = lambda x: tm.relative_rate(x, 1.43, 7.4, eps_miss=0.2)
    out = stats.simulate_scan(np.arange(-4.5, 4.6, 1.5), 500, R, 0.3, n_toys=200, fit="template", rng=1)
    assert abs(out["bias"]) < 0.03 and out["rms"] < 0.15


def test_fit_batch_edges_recovers_offset():
    rng = np.random.default_rng(2)
    t = 1234.0 + rng.uniform(0, numi.SPILL_NS, 20000) + rng.normal(0, 10, 20000)
    t1, _ = dr.fit_batch_edges(np.concatenate([t, rng.uniform(-5e3, 2e4, 500)]), 1100.0)
    assert abs(t1 - 1234.0) < 30


def test_match_spills():
    idx = dr.match_spills([10.0, 11.21, 50.0], [9.9, 11.2, 12.4])
    assert list(idx) == [0, 1, -1]


def _write_synthetic(stem, n=50, dx1=-0.3):
    rows, q = [], []
    rng = np.random.default_rng(3)
    for i in range(n):
        t = 1.72e9 + 1.2 * i
        jit = rng.normal(0, 0.05)
        for dev in ifbeam.ARRAY_DEVICES:
            for k in range(7):
                if "I" in dev.split(":")[1][:2]:          # HITGT / VITGT intensities
                    v = 190.0
                else:
                    v = 0.1 + jit + (dx1 if (k == 1 and dev.startswith("E:H")) else 0.0)
                rows.append((dev, t + 0.001, k, v))
        q.append(dict(time_s=t, spillpot=5e13, hornI=-200.0, posx=0.0, posy=0.0, widthx=1.4,
                      widthy=1.4, goodbeam=1))
    pd.DataFrame(rows, columns=["device", "time_s", "index", "value"]).to_csv(f"{stem}.csv", index=False)
    pd.DataFrame(q).to_csv(f"{stem}_quality.csv", index=False)


def test_ifbeam_loader(tmp_path):
    stem = tmp_path / "get_data_synth"
    _write_synthetic(stem)
    b = ifbeam.load_batches(stem)
    assert len(b) == 50 * 6
    assert np.allclose(b.groupby("spill")["frac"].sum(), 1.0)
    d = ifbeam.batch1_minus_rest(b)
    # both BPM stations carry the same -0.3 mm batch-1 shift -> -0.3 mm at target
    assert np.allclose(d.dx, -0.3, atol=1e-9) and np.allclose(d.dy, 0.0, atol=1e-9)


def test_mx2_slices_from_chunk():
    """Slices >= 1 only; light-weighted time; tracks attached to their slice."""
    ak = pytest.importorskip("awkward")
    from bpscan import mx2_tracks as mx
    ch = dict(ev_gps_time_sec=ak.Array([100, 101]), ev_gps_time_usec=ak.Array([5, 7]),
              ev_trigger_type=ak.Array([16, 16]), ev_gate=ak.Array([1, 2]), n_slices=ak.Array([2, 1]),
              clus_id_time=ak.Array([[10., 20., 500., 9000.], [30.]]),
              clus_id_time_slice=ak.Array([[1, 1, 2, 0], [1]]),
              clus_id_pe=ak.Array([[1., 3., 2., 5.], [1.]]), clus_id_z=ak.Array([[1., 5., 3., 0.], [2.]]),
              n_tracks=ak.Array([2, 0]), trk_time_slice=ak.Array([[1, 1], []]),
              trk_hits=ak.Array([[10, 20], []]), trk_nodes=ak.Array([[5, 9], []]),
              trk_theta=ak.Array([[0.1, 0.3], []]))
    g, s = mx.slices_from_chunk(ch, 0)
    assert len(g) == 2 and abs(g.utime_ms[0] - 100000.005) < 1e-6
    assert len(s) == 3                                   # slice 0 dropped
    first = s[(s.event == 0) & (s.slice == 1)].iloc[0]
    assert abs(first.t_rel - 17.5) < 1e-9 and first.ntrk == 2 and first.trk_hits_max == 20
    assert first.zmax - first.zmin == 4.0
    assert s[(s.event == 0) & (s.slice == 2)].ntrk.iloc[0] == 0
