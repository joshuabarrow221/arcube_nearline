"""Human-readable bad-spills report from a get_data quality CSV.

Reads get_data_<window>_quality.csv (written by get_data.cpp alongside
the 'quality' TTree) and writes get_data_<window>_bad_spills.txt next
to it: one block per spill failing the NOvA goodbeam cuts 1-5, with the
failing cuts spelled out and the likely reason derived from which
quantities are missing vs. out of range. See the "Beam-quality cuts"
comment block in get_data.cpp for the cut catalogue.

usage: python3 bad_spills_report.py [get_data_<t0>_to_<t1>] ...
"""
import sys
import numpy as np
import pandas as pd
from datetime import datetime, timedelta, timezone

BASE = os.environ.get("BEAMQ_DIR", os.getcwd())
CDT = timezone(timedelta(hours=-5))
CUTS = {  # NOvA standard_ifdbspillinfo values (see get_data.cpp)
    "pot":  ("POT", "> 2.00e12 protons"),
    "horn": ("horn current", "in (-202, -196.4) kA"),
    "posx": ("x position", "in (-2, +2) mm"),
    "posy": ("y position", "in (-2, +2) mm"),
    "widthx": ("x width", "in (0.57, 1.58) mm"),
    "widthy": ("y width", "in (0.57, 1.58) mm"),
}

def fmt(v, unit=""):
    return "MISSING" if pd.isna(v) else f"{v:.4g}{unit}"

def reason(row):
    if pd.isna(row.spillpot):
        return ("no toroid reading in IFBeam at this spill (E:TRTGTD and "
                "E:TR101D both absent) -> POT unknown, POT cut fails")
    if row.spillpot < 1e11 and pd.isna(row.posx) and pd.isna(row.widthx):
        return ("near-zero POT 'ghost' spill (<1e11 protons); BPM and "
                "multiwire report no usable data")
    bad = [k for k in CUTS if row[f"pass_{'pot' if k=='pot' else k}"] == 0]
    return "out-of-range: " + ", ".join(CUTS[k][0] for k in bad)

def report(stem):
    q = pd.read_csv(f"{BASE}/{stem}_quality.csv")
    bad = q[q.goodbeam == 0]
    out = f"{BASE}/{stem}_bad_spills.txt"
    with open(out, "w") as f:
        f.write(f"Bad-spill report for {stem}\n")
        f.write("goodbeam = NOvA cuts 1-5 (POT, horn current, x/y position, x/y width);\n")
        f.write("the trigger-timing cut 6 is not evaluable without a detector trigger stream.\n")
        f.write(f"Cut values: " + "; ".join(f"{n} {r}" for n, r in CUTS.values()) + "\n")
        f.write(f"\n{len(bad)} of {len(q)} spills fail goodbeam. "
                f"POT lost: {bad.spillpot.fillna(0).sum():.4g} of {q.spillpot.fillna(0).sum():.4g} "
                f"({100*bad.spillpot.fillna(0).sum()/q.spillpot.fillna(0).sum():.4f}%)\n\n")
        for _, r in bad.iterrows():
            t = datetime.fromtimestamp(r.time_s, CDT)
            fails = [CUTS[k][0] for k in CUTS
                     if r[f"pass_{'pot' if k == 'pot' else k}"] == 0]
            f.write(f"spill @ {r.time_s:.3f}  ({t:%Y-%m-%d %H:%M:%S.%f} CDT"[:60] + ")\n")
            f.write(f"  POT {fmt(r.spillpot)}  hornI {fmt(r.hornI,' kA')}  "
                    f"pos ({fmt(r.posx)}, {fmt(r.posy)}) mm  "
                    f"width ({fmt(r.widthx)}, {fmt(r.widthy)}) mm\n")
            f.write(f"  fails: {', '.join(fails)}\n")
            f.write(f"  reason: {reason(r)}\n\n")
    print(f"wrote {out} ({len(bad)} bad spills)")

if __name__ == "__main__":
    stems = sys.argv[1:] or [
        "get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500",
        "get_data_2024-07-09T000000-0500_to_2024-07-09T020500-0500",
    ]
    for s in stems:
        report(s)
