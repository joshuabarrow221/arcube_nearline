#!/usr/bin/env python3
"""Collect every number quoted in the report into report/numbers.tex.

The LaTeX source never hard-codes a result: it uses macros such as
\\BoneDx defined here from results/*/*.json. Rerun after any analysis rerun:
    python3 scripts/make_report_numbers.py && (cd report && pdflatex ... )
Missing inputs produce a visible '??' rather than a stale number.
"""
import json
from pathlib import Path

R = Path(__file__).resolve().parents[1] / "results"
out = []


def mac(name, value, fmt="{:.3g}"):
    try:
        s = fmt.format(value)
    except Exception:
        s = "??"
    out.append(f"\\newcommand{{\\{name}}}{{{s}}}")


def load(p):
    try:
        return json.loads((R / p).read_text())
    except Exception:
        return {}


ib = load("ifbeam/ifbeam_batches.json")
b1 = ib.get("batch1_minus_rest", {})
mac("NspillGood", ib.get("n_spills_good"), "{:,}")
mac("NspillTot", ib.get("n_spills_total"), "{:,}")
mac("PotGood", ib.get("pot_good", 0) / 1e17, "{:.2f}")
mac("BoneDx", b1.get("dx_mean"), "{:+.3f}")
mac("BoneDxRms", b1.get("dx_rms"), "{:.3f}")
mac("BoneDxSem", b1.get("dx_sem"), "{:.4f}")
mac("BoneDy", b1.get("dy_mean"), "{:+.3f}")
mac("BoneDyRms", b1.get("dy_rms"), "{:.3f}")
mac("SigX", ib.get("widths", {}).get("sigx_median"), "{:.2f}")
mac("SigY", ib.get("widths", {}).get("sigy_median"), "{:.2f}")
sp = ib.get("spill_pos", {})
mac("PosX", sp.get("posx_mean"), "{:+.2f}"); mac("PosXrms", sp.get("posx_rms"), "{:.3f}")
mac("PosY", sp.get("posy_mean"), "{:+.2f}"); mac("PosYrms", sp.get("posy_rms"), "{:.2f}")
mac("HornI", ib.get("hornI", {}).get("median"), "{:.1f}")
pb = {r["batch"]: r for r in ib.get("per_batch", [])}
for k, w in zip(range(1, 7), "one two three four five six".split()):
    mac(f"Bx{w}", pb.get(k, {}).get("x_mean"), "{:+.3f}")
    mac(f"Bfrac{w}", 100 * pb.get(k, {}).get("frac_mean", float("nan")), "{:.2f}")
lev = ib.get("lever_arms", {})
mac("LevBoneXhalf", lev.get("B_batch1_x_halves", {}).get("dx"), "{:+.3f}")
mac("LevRestXhalf", lev.get("B_batches2to6_x_halves", {}).get("dx"), "{:+.3f}")
mac("LevRestYhalf", lev.get("B_batches2to6_y_halves", {}).get("dy"), "{:+.3f}")
mac("LevD", lev.get("D_batch1_left_vs_2to6_right", {}).get("dx"), "{:+.3f}")
pr = ib.get("predicted_fractional_rate_difference", {}).get("A_batch1_vs_2to6", {})
mac("PredAzero", 100 * pr.get("fin_centre_at_+0.0mm", float("nan")), "{:+.2f}")
mac("PredAmone", 100 * pr.get("fin_centre_at_-1.0mm", float("nan")), "{:+.2f}")
mac("PredAone", 100 * pr.get("fin_centre_at_+1.0mm", float("nan")), "{:+.2f}")

se = load("sensitivity/sensitivity.json")
mac("Curv", se.get("curvature_per_mm2"), "{:.4f}")
for k, n in (("1.0", "Rone"), ("2.0", "Rtwo"), ("3.0", "Rthree"), ("4.5", "Rfourhalf")):
    mac(n, se.get("R_at", {}).get(k), "{:.3f}")
md = se.get("min_delta_3sigma_at", {})
mac("MinDeltaOneEfive", 100 * md.get("1e+05", float("nan")), "{:.1f}")
mac("MinDeltaOneEsix", 100 * md.get("1e+06", float("nan")), "{:.2f}")
w = se.get("scan_wide_template", {})
for n, name in (("50", "ScanFifty"), ("200", "ScanTwoH"), ("500", "ScanFiveH"), ("2000", "ScanTwoK")):
    mac(name, w.get(n, {}).get("rms"), "{:.2f}")
nar = se.get("scan_narrow_parabola", {})
mac("NarrowTwoKfail", 100 * nar.get("2000", {}).get("frac_fail", float("nan")), "{:.0f}")
mac("NarrowTwoKrms", nar.get("2000", {}).get("rms"), "{:.2f}")
wt = se.get("scan_wide_wrong_template", {})
mac("WrongBias", wt.get("2000", {}).get("bias"), "{:+.2f}")
tb = se.get("time_budget_100kW", {})
mac("CountsQuarter", tb.get("0.25", {}).get("counts_per_point"), "{:.0f}")
mac("CountsHalf", tb.get("0.5", {}).get("counts_per_point"), "{:.0f}")
mac("CountsTenth", tb.get("0.1", {}).get("counts_per_point"), "{:.0f}")

br = load("batch_rates/batch_rates.json")
au = br.get("audit", {})
mac("BRspills", au.get("spills_used"), "{:,}")
mac("BRtrig", au.get("beam_triggers"), "{:,}")
mac("BRunmatched", au.get("unmatched_triggers"), "{:,}")
mac("BRdropT", au.get("dropped_timing_resid"), "{:,}")
mac("BRdropGB", au.get("dropped_not_goodbeam"), "{:,}")
mac("BRdropDup", au.get("dropped_duplicate_match"), "{:,}")
mac("BRoffset", au.get("trigger_minus_A9_offset_ms"), "{:.1f}")
mac("BRresid", au.get("trigger_resid_rms_ms"), "{:.1f}")
mac("BRtone", au.get("batch1_start_ns"), "{:.0f}")
mac("BRflashes", au.get("flashes_selected"), "{:,}")
r = br.get("batch1_over_rest", {})
mac("BRratio", r.get("ratio"), "{:.3f}")
mac("BRerr", r.get("err"), "{:.3f}")
mac("BRpval", r.get("pvalue"), "{:.2f}")
mac("BRmindelta", 100 * r.get("min_detectable_3sigma", float("nan")), "{:.0f}")
mac("BRnone", r.get("n1"), "{:,.0f}")
mac("BRnrest", r.get("n2"), "{:,.0f}")
fr = load("batch_rates/flash_rate_summary.json")
FMT = {"FRrunOneMinDelta": "{:.1f}", "FRpowerKw": "{:.0f}", "FRperHourScan": "{:.0f}",
       "FRrockPerHourScan": "{:.0f}", "FRperSpill": "{:.2f}", "FRrockPerSpill": "{:.2f}"}
for k, v in fr.items():
    if k == "FRrunOneTotal":
        e = int(__import__("math").floor(__import__("math").log10(v)))
        out.append(f"\\newcommand{{\\{k}}}{{{v / 10**e:.1f}\\times10^{{{e}}}}}")
    else:
        mac(k, v, FMT.get(k, "{:.3g}"))
(Path(__file__).resolve().parents[1] / "report" / "numbers.tex").write_text("\n".join(out) + "\n")
print("\n".join(out))


def _fmt_time(h):
    """Human time per scan point."""
    if not h == h:
        return "??"
    m = h * 60
    if m < 1:
        return "$<$1\\,min"
    if m < 90:
        return f"{m:.0f}\\,min"
    return f"{h:.1f}\\,h"


rates = {"All": fr.get("FRperHourScan", float("nan")), "Rock": fr.get("FRrockPerHourScan", float("nan")), "Nu": 50.0}
extra = []
for prec, pname in (("0.5", "Half"), ("0.25", "Quarter"), ("0.1", "Tenth")):
    n = tb.get(prec, {}).get("counts_per_point", float("nan"))
    for rname, r in rates.items():
        extra.append(f"\\newcommand{{\\T{pname}{rname}}}{{{_fmt_time(n / r)}}}")
with open(Path(__file__).resolve().parents[1] / "report" / "numbers.tex", "a") as fh:
    fh.write("\n".join(extra) + "\n")
print("\n".join(extra))
