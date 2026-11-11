"""NuMI spill time structure and BPM geometry.

Time structure (Main Injector, slip-stacked operation)
------------------------------------------------------
* MI RF frequency ~53.10 MHz  ->  bucket spacing ~18.83 ns.
* The MI circumference holds 588 buckets (~11.1 us). It is divided into 7
  Booster-batch slots of 84 buckets (~1.58 us each). For NuMI, 6 slots are
  filled ("batches 1-6") and the 7th is left empty as the abort-kicker gap.
* A NuMI spill is therefore ~9.5 us of beam made of six ~1.58 us batches.

These are nominal numbers. The *measured* batch boundaries in any detector
must be calibrated from the data (e.g. the distribution of LRS flash times
relative to the RWM edge shows six plateaus); see detector_rates.fit_batch_edges.
The nominal values are used only as a starting point and for guard-band sizing.

BPM geometry
------------
Copied verbatim from beam_quality/get_data.cpp (itself a port of NOvA's
IFDBSpillInfo BpmProjection). z positions in feet, origin between the two
Budal monitors of the target; only ratios matter for linear extrapolation.
"""
import numpy as np

RF_HZ = 53.103e6                 # Main Injector RF at extraction [Hz]
BUCKET_NS = 1e9 / RF_HZ          # ~18.83 ns
BUCKETS_PER_BATCH = 84
N_SLOTS = 7                      # 7 Booster-batch slots in the MI ring
N_BATCHES = 6                    # filled slots for NuMI
BATCH_NS = BUCKETS_PER_BATCH * BUCKET_NS   # ~1581.8 ns
SPILL_NS = N_BATCHES * BATCH_NS            # ~9491 ns

# Same constants as beam_quality/get_data.cpp (kZhp121, ...). Units: feet.
Z_HP121, Z_VP121 = -68.04458, -66.99283
Z_HPTGT, Z_VPTGT = -31.25508, -30.16533
Z_TARGET = 0.0


def extrapolate(t1, z1, t2, z2, z3=Z_TARGET):
    """Linear extrapolation of a transverse coordinate measured at z1, z2 to z3.

    Identical to extrapolatePosition() in get_data.cpp. Works on numpy arrays.
    """
    t1 = np.asarray(t1, float)
    t2 = np.asarray(t2, float)
    return t1 + (t2 - t1) * (z3 - z1) / (z2 - z1)


def nominal_batch_edges_ns(t_first_ns=0.0):
    """Nominal start times of batches 1..6 and the end of batch 6 [ns].

    Parameters
    ----------
    t_first_ns : float
        Time (in the detector clock, relative to whatever reference the caller
        uses, e.g. the RWM leading edge) at which batch 1 begins. Must be
        calibrated per detector; it absorbs cable delays and the ~3.4 us
        time of flight over ~1 km.

    Returns
    -------
    ndarray of shape (7,)
    """
    return t_first_ns + BATCH_NS * np.arange(N_BATCHES + 1)


def batch_index(t_rel_ns, edges_ns, guard_ns=0.0):
    """Assign each time to a batch number 1..6, or 0 when ambiguous/outside.

    An object is assigned to batch b if  edges[b-1] + guard <= t < edges[b] - guard.
    Times inside a guard band around any boundary, before batch 1 or after
    batch 6, get 0 and must be excluded from per-batch rate comparisons
    (they are *not* silently assigned to the nearest batch).

    guard_ns should be set to ~3x the per-object timing resolution relative to
    the spill reference (LRS: tens of ns; CRS anode-crossing t0: ~100-300 ns;
    MINERvA: ~few ns), plus the rise time of the batch-1 kicker waveform when
    cutting "carefully" around batch 1 as requested.
    """
    t = np.asarray(t_rel_ns, float)
    edges = np.asarray(edges_ns, float)
    idx = np.searchsorted(edges, t, side="right")          # 1..6 inside spill
    out = np.where((idx >= 1) & (idx <= N_BATCHES), idx, 0)
    if guard_ns > 0:
        lo = edges[np.clip(out - 1, 0, N_BATCHES)]
        hi = edges[np.clip(out, 0, N_BATCHES)]
        bad = (t < lo + guard_ns) | (t >= hi - guard_ns)
        out = np.where((out > 0) & ~bad, out, 0)
    return out.astype(int)


def beam_power_kw(pot_per_spill, cycle_s=1.2, proton_gev=120.0):
    """Instantaneous NuMI beam power [kW] for a given POT/spill and cycle time."""
    joules = pot_per_spill * proton_gev * 1.602176634e-10
    return joules / cycle_s / 1e3
