"""Per-batch, per-POT detector rates.

This module is deliberately data-format agnostic. A reader for a particular
file type (2x2 FLOW HDF5, MINERvA/Mx2 ROOT, ...) must produce, for each
recorded beam spill:

    spill_time_s   unix time of the spill's detector trigger [s]
    t_rel_ns       time of each counted object relative to the spill reference
                   (e.g. RWM leading edge digitised in the light system) [ns]
    weights        optional per-object weights (default 1)

and the caller provides the per-batch POT table from ifbeam.load_batches().

Timing quality is enforced here, not trusted: spills whose detector trigger is
not matched to an IFBeam $A9 record within `match_tol_s`, or whose reference
(RWM) time is missing/outlying, are dropped and counted in the returned audit.
This follows the analysis requirement that any spill with timing drift between
detector and accelerator clocks is removed rather than corrected.

Counting objects: what to count is a physics choice made by the reader:
  * 2x2 LRS: light flashes (sum-hit clusters) above threshold, in-spill.
  * 2x2 CRS: anode-crossing track segments (their earliest pixel time is the
    interaction time to ~100 ns), i.e. through-going rock muons.
  * Mx2: hits/tracks in the front veto planes (MINERvA timing ~ns).
"""
import numpy as np
import pandas as pd

from . import numi


def fit_batch_edges(t_rel_ns, t_first_guess_ns, search_ns=400.0, step_ns=2.0):
    """Calibrate where batch 1 starts in a detector's clock.

    Scans the start time t1 over +-search_ns around the guess and maximises the
    number of objects falling inside the nominal 6-batch window [t1, t1+SPILL_NS).
    Returns (t1_best, edges). With O(1e4) in-spill objects the window edges are
    found to a few tens of ns. Always inspect the time histogram as well
    (scripts/plot_batch_times.py) - this function only aligns, it does not
    validate that six batches are visible.
    """
    t = np.sort(np.asarray(t_rel_ns, float))
    grid = np.arange(t_first_guess_ns - search_ns, t_first_guess_ns + search_ns + step_ns, step_ns)
    inside = (np.searchsorted(t, grid + numi.SPILL_NS) - np.searchsorted(t, grid))
    # choose the centre of the plateau of maximal containment
    best = grid[inside == inside.max()]
    t1 = 0.5 * (best.min() + best.max())
    return t1, numi.nominal_batch_edges_ns(t1)


def match_spills(det_spill_times_s, beam_spill_times_s, match_tol_s=0.5):
    """Index of the nearest IFBeam spill for each detector spill, or -1."""
    beam = np.asarray(beam_spill_times_s, float)
    order = np.argsort(beam)
    bs = beam[order]
    det = np.asarray(det_spill_times_s, float)
    j = np.clip(np.searchsorted(bs, det), 1, len(bs) - 1)
    left, right = bs[j - 1], bs[j]
    pick = np.where(np.abs(det - left) <= np.abs(det - right), j - 1, j)
    ok = np.abs(det - bs[pick]) <= match_tol_s
    return np.where(ok, order[pick], -1)


def per_batch_counts(spill_ids, t_rel_ns, edges_ns, guard_ns, weights=None):
    """Counts per (spill, batch) using numi.batch_index with guard bands."""
    b = numi.batch_index(t_rel_ns, edges_ns, guard_ns)
    w = np.ones_like(b, float) if weights is None else np.asarray(weights, float)
    df = pd.DataFrame({"spill": spill_ids, "batch": b, "w": w})
    df = df[df["batch"] > 0]
    return df.groupby(["spill", "batch"])["w"].sum().rename("n").reset_index()


def guard_acceptance(edges_ns, guard_ns, batch):
    """Fraction of a batch's duration surviving the guard band (for POT scaling)."""
    width = edges_ns[batch] - edges_ns[batch - 1]
    return max(width - 2 * guard_ns, 0.0) / width


def rates_by_group(counts, batches, group_col, edges_ns, guard_ns):
    """Sum counts and POT by an arbitrary grouping of (spill, batch) rows.

    counts   output of per_batch_counts (spill ids must be IFBeam spill ids)
    batches  ifbeam.load_batches() rows restricted to *detector-live, matched,
             goodbeam* spills (the caller decides the selection; spills with no
             detector record must not be included, or the rate is biased low)
    group_col column of `batches` defining the comparison, e.g. 'is_batch1' or
             a position-quantile label ('x_quartile').
    POT per (spill,batch) is scaled by the guard-band acceptance so that the
    rate is per *effective* POT.
    """
    b = batches.copy()
    b["acc"] = [guard_acceptance(edges_ns, guard_ns, int(k)) for k in b["batch"]]
    b["pot_eff"] = b["pot_batch"] * b["acc"]
    m = b.merge(counts, on=["spill", "batch"], how="left").fillna({"n": 0.0})
    g = m.groupby(group_col).agg(n=("n", "sum"), pot=("pot_eff", "sum"), rows=("n", "size"))
    g["rate_per_1e13"] = g["n"] / g["pot"] * 1e13
    g["rate_err"] = np.sqrt(g["n"]) / g["pot"] * 1e13
    return g
