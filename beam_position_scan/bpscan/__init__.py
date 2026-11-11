"""bpscan: NuMI proton-beam position on target vs. 2x2 / Mx2 detector rates.

Package layout (see ../README.md for the full guide):

    numi.py           NuMI spill structure (RF buckets, Booster batches) and
                      BPM geometry constants shared with ../../beam_quality.
    ifbeam.py         Load the long-format CSV written by beam_quality/get_data.cpp
                      into tidy per-spill, per-batch tables (positions at the
                      target z, per-batch intensity fractions, per-batch POT).
    target_model.py   Geometric model: fraction of a Gaussian proton beam that
                      intercepts a finite-width target fin, and the curvature of
                      the resulting rate response around the optimum.
    stats.py          Poisson sensitivity: two-sample rate-ratio tests,
                      required counts, and Monte Carlo scan planning
                      (how well a scan locates an optimum / an edge).
    detector_rates.py Per-batch, per-POT counting of detector objects whose
                      times are known relative to the spill (LRS flashes,
                      anode-crossing CRS tracks, Mx2 hits). Data-format
                      agnostic: callers supply (spill_time, t_rel_ns) arrays.

Nothing here invents detector data. Functions that need detector inputs take
arrays and refuse to run (raise) when required quantities are missing.
"""
__all__ = ["numi", "ifbeam", "target_model", "stats", "detector_rates"]
