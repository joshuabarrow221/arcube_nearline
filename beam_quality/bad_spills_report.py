"""Human-readable bad-spills report from a get_data quality CSV.

Reads get_data_<window>_quality.csv (written by get_data.cpp alongside
the 'quality' TTree) and writes get_data_<window>_bad_spills.txt next
to it: one block per spill failing the NOvA goodbeam cuts 1-5, with the
failing cuts spelled out and the likely reason derived from which
quantities are missing vs. out of range. The cut values quoted are the
ones get_data.cpp recorded in get_data_<window>_cuts.json (see
beam_quality_cuts.py); the "Beam-quality cuts" comment block in
get_data.cpp documents the catalogue and NOvA's provenance for each cut.

usage: python3 bad_spills_report.py [get_data_<t0>_to_<t1>] ...
"""
import os
import sys
import pandas as pd
from datetime import datetime, timedelta, timezone

import beam_quality_cuts as bqc

BASE = os.environ.get("BEAMQ_DIR", os.getcwd())
CDT = timezone(timedelta(hours=-5))


def cut_catalogue(c):
    """Human-readable requirement for each pass flag, from the recorded cuts."""
    return {  # key -> (name, requirement); key k corresponds to column pass_<k>
        "pot":    ("POT", f"> {c['minPOT']:.2e} protons"),
        "horn":   ("horn current", f"in ({c['minHornI']:g}, {c['maxHornI']:g}) kA"),
        "posx":   ("x position", f"in ({c['minPosX']:+.2f}, {c['maxPosX']:+.2f}) mm"),
        "posy":   ("y position", f"in ({c['minPosY']:+.2f}, {c['maxPosY']:+.2f}) mm"),
        "widthx": ("x width", f"in ({c['minWidthX']:.2f}, {c['maxWidthX']:.2f}) mm"),
        "widthy": ("y width", f"in ({c['minWidthY']:.2f}, {c['maxWidthY']:.2f}) mm"),
    }


def fmt(v, unit=""):
    return "MISSING" if pd.isna(v) else f"{v:.4g}{unit}"


def reason(row, cuts):
    """Most likely cause of a failure, judged from what is missing vs. out of range."""
    if pd.isna(row.spillpot):
        return ("no toroid reading in IFBeam at this spill (E:TRTGTD and "
                "E:TR101D both absent) -> POT unknown, POT cut fails")
    if row.spillpot < 1e11 and pd.isna(row.posx) and pd.isna(row.widthx):
        return ("near-zero POT 'ghost' spill (<1e11 protons); BPM and "
                "multiwire report no usable data")
    bad = [k for k in cuts if row[f"pass_{k}"] == 0]
    return "out-of-range: " + ", ".join(cuts[k][0] for k in bad)


def report(stem):
    q = pd.read_csv(f"{BASE}/{stem}_quality.csv")
    recorded = bqc.load_cuts(BASE, stem)
    cuts = cut_catalogue(recorded)
    bad = q[q.goodbeam == 0]
    out = f"{BASE}/{stem}_bad_spills.txt"
    with open(out, "w") as f:
        f.write(f"Bad-spill report for {stem}\n")
        f.write("goodbeam = NOvA cuts 1-5 (POT, horn current, x/y position, x/y width);\n")
        f.write("the trigger-timing cut 6 is not evaluable without a detector trigger stream.\n")
        f.write(f"Cut values (NOvA epoch {recorded['epoch']}): "
                + "; ".join(f"{n} {r}" for n, r in cuts.values()) + "\n")
        if "isRHC" in q:
            f.write(f"Horn polarity from E:HRNDIR: {int(q.isRHC.sum())} of {len(q)} spills RHC "
                    "(E:HRNDIR > 0); the horn-current window applies to both polarities.\n")
        f.write(f"\n{len(bad)} of {len(q)} spills fail goodbeam. "
                f"POT lost: {bad.spillpot.fillna(0).sum():.4g} of {q.spillpot.fillna(0).sum():.4g} "
                f"({100*bad.spillpot.fillna(0).sum()/q.spillpot.fillna(0).sum():.4f}%)\n\n")
        for _, r in bad.iterrows():
            t = datetime.fromtimestamp(r.time_s, CDT)
            fails = [cuts[k][0] for k in cuts if r[f"pass_{k}"] == 0]
            f.write(f"spill @ {r.time_s:.3f}  ({t:%Y-%m-%d %H:%M:%S.%f} CDT"[:60] + ")\n")
            f.write(f"  POT {fmt(r.spillpot)}  hornI {fmt(r.hornI,' kA')}  "
                    f"pos ({fmt(r.posx)}, {fmt(r.posy)}) mm  "
                    f"width ({fmt(r.widthx)}, {fmt(r.widthy)}) mm\n")
            f.write(f"  fails: {', '.join(fails)}\n")
            f.write(f"  reason: {reason(r, cuts)}\n\n")
    print(f"wrote {out} ({len(bad)} bad spills)")


if __name__ == "__main__":
    stems = sys.argv[1:] or [
        "get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500",
        "get_data_2024-07-09T000000-0500_to_2024-07-09T020500-0500",
    ]
    for s in stems:
        report(s)
