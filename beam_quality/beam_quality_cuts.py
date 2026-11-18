"""Shared loader for the goodbeam cut values that get_data.cpp applied.

get_data.cpp records the exact NOvA cut values it evaluated in
<stem>_cuts.json (and as cut_* TParameter objects in <stem>.root), so
every plot and report in this directory draws the cuts that were
actually applied instead of a hard-coded copy that can drift from the
evaluation. (Before this record existed, the scripts hard-coded a
symmetric +-2 mm box and 1.58 mm widths — NOvA's periods 1-10 values,
not the periods 13+ values that apply to all 2x2 running.)

Outputs written before the cut record existed fall back to NOvA's
current values (periods 13+), with a warning. For July 2024 (2x2 Run 1)
the old and current values give identical verdicts on every spill.
"""
import json
import os
import warnings

# NOvA IFDBSpillInfo.fcl standard_ifdbspillinfo values for NOvA periods
# 13+ (all 2x2 running), converted to absolute protons, kA, mm and s —
# identical to novaCutsForEpoch("P13+") in get_data.cpp, whose header
# comment restates where each value comes from.
NOVA_P13PLUS = dict(
    epoch="P13+",
    minPOT=2.00e12,                     # POT > 2e12 protons
    minHornI=-202.0, maxHornI=-196.4,   # kA, same window for FHC and RHC
    hornI0HC=1.0,                       # kA, |hornI| below this = horn off
    minPosX=-2.00, maxPosX=2.00,        # mm at target
    minPosY=-3.50, maxPosY=0.50,        # mm at target (NOvA docdb-57045)
    minWidthX=0.57, maxWidthX=1.88,     # mm
    minWidthY=0.57, maxWidthY=1.88,     # mm (NOvA docdb-57045)
    maxDeltaT=0.5,                      # s, trigger-to-spill timing (cut 6)
)


def load_cuts(base, stem):
    """Return the cut dictionary recorded for <base>/<stem>, or NOvA P13+."""
    path = os.path.join(base, f"{stem}_cuts.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    warnings.warn(f"{path} not found (this output predates the cut record); "
                  "assuming NOvA's current periods-13+ cut values")
    return dict(NOVA_P13PLUS)


def position_box(cuts):
    """(x0, y0, width, height) of the cuts 3-4 position window at the target [mm]."""
    return (cuts["minPosX"], cuts["minPosY"],
            cuts["maxPosX"] - cuts["minPosX"], cuts["maxPosY"] - cuts["minPosY"])


def position_label(cuts):
    """Short in-plot label for the position window. Kept short so it never runs
    under a legend; the box edges themselves show the window on the axes."""
    return f"NOvA goodbeam window ({cuts['epoch']})"


def position_window_text(cuts):
    """Full numeric description of the position window, for titles/captions."""
    return (f"x {cuts['minPosX']:+.2f}..{cuts['maxPosX']:+.2f} mm, "
            f"y {cuts['minPosY']:+.2f}..{cuts['maxPosY']:+.2f} mm")


def window_limits(values, lo_cut, hi_cut, pad=0.06):
    """Axis limits spanning both the data and the cut window, never clipping either."""
    lo = min(min(values), lo_cut)
    hi = max(max(values), hi_cut)
    p = (hi - lo) * pad
    return lo - p, hi + p
