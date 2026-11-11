"""Geometric model of a Gaussian proton beam on a finite-width target fin.

Only the *interception* effect is modelled here: the fraction of protons whose
transverse position falls inside the fin. This is the one piece of the
position dependence that can be computed from first principles without a full
beamline simulation (g4numi). Horn-focusing changes from an off-centre
production point are NOT in this model; they are handled as an external,
literature-based slope (see report section on systematics) and must be
replaced by g4numi derivatives before any quantitative scan decision.

Conventions
-----------
x0      beam centroid relative to the fin centre [mm]
sigma   Gaussian beam RMS in that plane [mm]
width   full fin width in that plane [mm]  (NOvA-era ME fins: 7.4 mm wide,
        63 mm tall; the 1 MW-era fin width is not confirmed here -> pass it)
"""
import numpy as np
from scipy.special import erf
from scipy.stats import norm


def intercept_fraction(x0, sigma, width):
    """Fraction of a 1-D Gaussian beam of RMS sigma centred at x0 inside |x|<width/2."""
    a = 0.5 * width
    x0 = np.asarray(x0, float)
    return norm.cdf((a - x0) / sigma) - norm.cdf((-a - x0) / sigma)


def curvature(sigma, width):
    """k [1/mm^2] such that f(x0)/f(0) ~= 1 - k x0^2 near the centre, i.e. -f''(0)/(2 f(0)).

    f''(0) = -2 (a/sigma^3) phi(a/sigma) with a = width/2.
    """
    a = 0.5 * width
    f0 = intercept_fraction(0.0, sigma, width)
    fpp = -2.0 * a / sigma**3 * norm.pdf(a / sigma)
    return -0.5 * fpp / f0


def relative_rate(x0, sigma, width, eps_miss=0.0, slope_per_mm=0.0, x_opt_focus=0.0,
                  quad_per_mm2=0.0):
    """Relative neutrino-equivalent rate vs beam position (1.0 at the optimum of the
    interception term when the other terms vanish).

    R(x0) = [f(x0) + eps_miss (1 - f(x0))] / [f(0) + eps_miss (1 - f(0))]
            * (1 - quad_per_mm2 * (x0 - x_opt_focus)^2) * (1 + slope_per_mm * x0)

    eps_miss      relative neutrino yield of protons that miss the fin (they still
                  hit downstream material, much less efficiently focused). 0..~0.3.
    quad_per_mm2  optional extra quadratic (focusing) term, e.g. from g4numi.
    slope_per_mm  optional linear term (pure placeholder; leave 0 unless measured).
    """
    f = intercept_fraction(x0, sigma, width)
    f0 = intercept_fraction(0.0, sigma, width)
    num = f + eps_miss * (1.0 - f)
    den = f0 + eps_miss * (1.0 - f0)
    x0 = np.asarray(x0, float)
    return num / den * (1.0 - quad_per_mm2 * (x0 - x_opt_focus) ** 2) * (1.0 + slope_per_mm * x0)


def two_position_difference(x_a, x_b, sigma, width, **kw):
    """Fractional rate difference R(x_a)/R(x_b) - 1 predicted by relative_rate()."""
    return relative_rate(x_a, sigma, width, **kw) / relative_rate(x_b, sigma, width, **kw) - 1.0
