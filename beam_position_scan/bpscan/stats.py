"""Poisson sensitivity tools for beam-position rate comparisons.

All rates are *per POT*: comparisons between batches or scan points must be
made on counts normalised by the POT delivered in that batch / at that point.
POT uncertainties (toroid ~1%, per-batch BPM-intensity share ~ few 0.1%) are
systematic and are carried separately via `pot_syst_frac`.

Key formulas (large-N Gaussian limit, verified by toy MC in tests/):

  Ratio r = (N1/P1)/(N2/P2).  sigma_r / r = sqrt(1/N1 + 1/N2 (+ pot syst^2)).

  Detecting a fractional difference delta at Z sigma between a sample with a
  fraction phi of the POT (batch 1: phi=1/6) and the rest (1-phi):
      N_total >= Z^2 / delta^2 * (1/phi + 1/(1-phi))
  For phi = 1/6:  N_total >= 7.2 Z^2 / delta^2.
"""
import numpy as np
from scipy import stats as st


def ratio_and_error(n1, pot1, n2, pot2, pot_syst_frac=0.0):
    """Per-POT rate ratio r=(n1/pot1)/(n2/pot2) and its 1-sigma error."""
    n1, n2 = float(n1), float(n2)
    r = (n1 / pot1) / (n2 / pot2)
    rel = np.sqrt(1.0 / n1 + 1.0 / n2 + pot_syst_frac**2)
    return r, r * rel


def exact_ratio_test(n1, pot1, n2, pot2):
    """Exact conditional test of equal per-POT rates.

    Given n = n1+n2, under H0 n1 ~ Binomial(n, p0) with p0 = pot1/(pot1+pot2).
    Returns (two-sided p-value, Clopper-Pearson 68% CI on the rate ratio).
    """
    n = int(n1 + n2)
    p0 = pot1 / (pot1 + pot2)
    pval = st.binomtest(int(n1), n, p0).pvalue
    ci = st.binomtest(int(n1), n, p0).proportion_ci(confidence_level=0.6827, method="exact")
    to_ratio = lambda p: (p / (1 - p)) * (pot2 / pot1)
    return pval, (to_ratio(ci.low), to_ratio(ci.high))


def required_total_counts(delta, z=3.0, phi=1.0 / 6.0, pot_syst_frac=0.0):
    """Total counts needed to see a fractional rate difference `delta` at `z` sigma
    when a fraction phi of the POT is in sample 1. Returns inf if the POT
    systematic alone exceeds delta/z."""
    stat_budget = (delta / z) ** 2 - pot_syst_frac**2
    if stat_budget <= 0:
        return np.inf
    return (1.0 / phi + 1.0 / (1.0 - phi)) / stat_budget


def min_detectable_delta(n_total, z=3.0, phi=1.0 / 6.0, pot_syst_frac=0.0):
    """Smallest fractional difference visible at z sigma with n_total counts."""
    return z * np.sqrt((1.0 / phi + 1.0 / (1.0 - phi)) / n_total + pot_syst_frac**2)


def counts_expected(rate_per_hour_at_ref, ref_power_kw, power_kw, hours, live_frac=1.0):
    """Linear scaling of a measured rate to another beam power and exposure.

    Detector rates scale with POT, i.e. with protons/spill x spills/hour. A
    lower *power* reached by fewer protons/pulse at fixed cycle scales linearly
    too. live_frac covers DAQ dead time and quality-cut losses.
    """
    return rate_per_hour_at_ref * (power_kw / ref_power_kw) * hours * live_frac


def simulate_scan(positions_mm, counts_at_opt, response, x_true, n_toys=2000, rng=None,
                  fit="parabola", pot_syst_frac=0.0):
    """Toy-MC precision of a position scan.

    positions_mm   scan points (beam centroid setpoints) [mm]
    counts_at_opt  expected counts per point if the beam were at the optimum
    response       callable R(x) -> relative rate (1 at optimum), e.g. a lambda
                   wrapping target_model.relative_rate
    x_true         true optimum offset relative to the scan coordinate [mm]
    fit            'parabola' : weighted quadratic fit, optimum = vertex
                   'centroid' : rate-weighted mean of positions (robust, biased
                                 for asymmetric scans)
                   'template' : Poisson maximum likelihood of the *known* shape
                                 response(x - xc) with free normalisation and
                                 free centre xc (grid search). Appropriate when
                                 sigma/width are measured (multiwire) and the
                                 scan reaches the fin edges.
    Returns dict(bias, rms, frac_fail) for the fitted optimum.
    """
    rng = np.random.default_rng(rng)
    x = np.asarray(positions_mm, float)
    mu = counts_at_opt * response(x - x_true)
    est = np.full(n_toys, np.nan)
    for i in range(n_toys):
        pot_scale = 1.0 + pot_syst_frac * rng.standard_normal(x.size)
        n = rng.poisson(mu * pot_scale)
        y = n / counts_at_opt
        w = 1.0 / np.maximum(n, 1) * counts_at_opt**2
        if fit == "template":
            grid = np.arange(x_true - 4.0, x_true + 4.0, 0.01)
            shapes = response(x[None, :] - grid[:, None])          # (G, P)
            A = n.sum() / shapes.sum(axis=1)                        # MLE norm
            lam = np.maximum(A[:, None] * shapes, 1e-12)
            ll = (n[None, :] * np.log(lam) - lam).sum(axis=1)
            est[i] = grid[np.argmax(ll)]
        elif fit == "parabola":
            c2, c1, c0 = np.polyfit(x, y, 2, w=np.sqrt(w))
            if c2 < 0:
                est[i] = -c1 / (2 * c2)
        else:
            est[i] = np.sum(x * n) / max(n.sum(), 1)
    ok = np.isfinite(est) & (np.abs(est - x_true) < 3 * np.ptp(x))
    e = est[ok]
    return dict(bias=float(np.mean(e) - x_true) if e.size else np.nan,
                rms=float(np.std(e)) if e.size else np.nan,
                frac_fail=float(1 - ok.mean()))
