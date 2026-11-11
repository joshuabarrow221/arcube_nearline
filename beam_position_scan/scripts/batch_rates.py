#!/usr/bin/env python3
"""Per-batch, per-POT 2x2 light-flash rates: batch 1 vs batches 2-6 (and per batch).

Inputs
  --flashes   *.flashes.npz from extract_lrs_flashes.py
  --ifbeam    get_data.cpp stems covering the same time span (per-batch POT,
              positions, goodbeam)
Pipeline (each step writes its counts to the audit in the output JSON)
  1. Beam triggers (trig_type==1) are matched to IFBeam $A9 spills using the
     light DAQ unix time (utime_ms). The trigger-minus-$A9 offset is fitted
     (median) and every spill whose residual exceeds --max-resid-ms is
     DROPPED (timing-drift rule), as are unmatched triggers and spills
     failing goodbeam.
  2. The batch-1 start in the light clock is fitted from the flash-time
     distribution (detector_rates.fit_batch_edges); six nominal 1581.8 ns
     batches follow.
  3. Flashes are assigned to batches with guard bands (--guard-ns, plus
     --guard-b1-extra-ns around batch 1, as requested: "cut around the batch 1
     timing window carefully").
  4. Rates = flashes / (per-batch POT x guard acceptance). Comparisons:
     batch 1 vs 2-6, each batch separately, and pile-up diagnostics: rate in
     batch k split by whether an earlier flash occurred in the same spill.

Flash selection (--min-nhits, --min-ntpc) is the "neutrino-equivalent" proxy:
  min-ntpc >= 2 favours through-going (rock-muon-like) activity.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import ifbeam, numi, stats, detector_rates as dr  # noqa: E402


def load_flashes(paths):
    trig, fl = [], []
    off = 0
    for k, pth in enumerate(paths):
        d = np.load(pth)
        b = d["trig_type"] == 1
        idx = np.where(b)[0]
        t = pd.DataFrame({"trig_row": idx + off, "utime_s": d["trig_utime_ms"][b] / 1e3, "file": k})
        trig.append(t)
        f = pd.DataFrame({c[6:]: d[c] for c in d.files if c.startswith("flash_")})
        f["trig_row"] = f["event"] + off
        fl.append(f)
        off += len(d["trig_type"])
    return pd.concat(trig, ignore_index=True), pd.concat(fl, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--flashes", nargs="+", required=True)
    p.add_argument("--ifbeam", nargs="+", required=True)
    p.add_argument("--out", default="results/batch_rates")
    p.add_argument("--max-resid-ms", type=float, default=10.0)
    p.add_argument("--guard-ns", type=float, default=50.0)
    p.add_argument("--guard-b1-extra-ns", type=float, default=50.0)
    p.add_argument("--min-nhits", type=int, default=3)
    p.add_argument("--min-ntpc", type=int, default=1)
    p.add_argument("--t1-guess-ns", type=float, default=1250.0)
    p.add_argument("--run1-pot", type=float, default=2.2e19,
                   help="Run 1 POT with CRS+LRS live (RunsDB mx2x2runs_v0.2: 2.2e19; LRS-only 3.0e19)")
    p.add_argument("--scan-power-kw", type=float, default=100.0)
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    audit = {}

    trig, fl = load_flashes(a.flashes)
    audit["beam_triggers"] = len(trig)
    frames = []
    for k, s in enumerate(a.ifbeam):
        b = ifbeam.load_batches(s); b["spill"] += 10_000_000 * k; frames.append(b)
    B = pd.concat(frames, ignore_index=True)
    spills = B.groupby("spill").first().reset_index()[["spill", "time_s", "goodbeam", "spillpot"]]

    # 1. spill matching + timing-drift veto -------------------------------
    j = dr.match_spills(trig.utime_s.values, spills.time_s.values, match_tol_s=0.5)
    trig["spill"] = np.where(j >= 0, spills.spill.values[np.maximum(j, 0)], -1)
    trig["resid_s"] = np.where(j >= 0, trig.utime_s.values - spills.time_s.values[np.maximum(j, 0)], np.nan)
    audit["unmatched_triggers"] = int((j < 0).sum())
    off = float(np.nanmedian(trig.resid_s))
    trig["resid_ms"] = (trig.resid_s - off) * 1e3
    bad_t = np.abs(trig.resid_ms) > a.max_resid_ms
    audit["trigger_minus_A9_offset_ms"] = off * 1e3
    audit["trigger_resid_rms_ms"] = float(np.nanstd(trig.resid_ms[~bad_t]))
    audit["dropped_timing_resid"] = int((bad_t & (j >= 0)).sum())
    dup = trig.spill.duplicated(keep=False) & (trig.spill >= 0)
    audit["dropped_duplicate_match"] = int(dup.sum())
    keep = (j >= 0) & ~bad_t & ~dup
    trig = trig[keep].merge(spills, on="spill")
    audit["dropped_not_goodbeam"] = int((trig.goodbeam != 1).sum())
    trig = trig[trig.goodbeam == 1]
    audit["spills_used"] = int(len(trig))

    # 2. flash selection and batch-edge calibration -------------------------
    f = fl[fl.trig_row.isin(trig.trig_row)].merge(trig[["trig_row", "spill"]], on="trig_row")
    f = f[(f.nhits >= a.min_nhits) & (f.ntpc >= a.min_ntpc)]
    t1, edges = dr.fit_batch_edges(f.t_rel.values, a.t1_guess_ns, search_ns=600)
    audit["batch1_start_ns"] = float(t1)
    audit["flashes_selected"] = int(len(f))
    # out-of-spill control region (before batch 1) -> flat background estimate
    pre = (f.t_rel > t1 - 1000) & (f.t_rel < t1 - 100)
    audit["background_flashes_per_ns"] = float(pre.sum() / 900.0)

    # 3. batch assignment with guards ---------------------------------------
    bi = numi.batch_index(f.t_rel.values, edges, a.guard_ns)
    extra = (np.abs(f.t_rel.values - edges[0]) < a.guard_ns + a.guard_b1_extra_ns) | \
            (np.abs(f.t_rel.values - edges[1]) < a.guard_ns + a.guard_b1_extra_ns)
    bi = np.where(extra, 0, bi)
    f["batch"] = bi
    counts = f[f.batch > 0].groupby(["spill", "batch"]).size().rename("n").reset_index()

    # 4. per-POT rates ----------------------------------------------------
    Bu = B[B.spill.isin(trig.spill) & B.pot_batch.notna()].copy()
    widths = np.diff(edges)
    acc = {b: max(widths[b - 1] - 2 * a.guard_ns - (2 * a.guard_b1_extra_ns if b == 1 else 0)
                  - (a.guard_b1_extra_ns if b == 2 else 0), 0) / widths[b - 1] for b in range(1, 7)}
    Bu["pot_eff"] = Bu.pot_batch * Bu.batch.map(acc)
    m = Bu.merge(counts, on=["spill", "batch"], how="left").fillna({"n": 0})
    per = m.groupby("batch").agg(n=("n", "sum"), pot=("pot_eff", "sum"), x=("x_tgt", "mean"), y=("y_tgt", "mean"))
    bkg = audit["background_flashes_per_ns"] / max(len(trig), 1)        # per spill per ns
    per["bkg"] = [bkg * (widths[b - 1] * acc[b]) * len(trig) for b in per.index]
    per["rate_e13"] = (per.n - per.bkg) / per.pot * 1e13
    per["err_e13"] = np.sqrt(per.n) / per.pot * 1e13
    res = dict(audit=audit, args=vars(a), per_batch=per.reset_index().to_dict(orient="records"))
    n1, p1 = per.loc[1, "n"] - per.loc[1, "bkg"], per.loc[1, "pot"]
    n2, p2 = (per.loc[2:6, "n"] - per.loc[2:6, "bkg"]).sum(), per.loc[2:6, "pot"].sum()
    r, e = stats.ratio_and_error(n1, p1, n2, p2)
    pval, ci = stats.exact_ratio_test(round(n1), p1, round(n2), p2)
    res["batch1_over_rest"] = dict(ratio=r, err=e, pvalue=pval, ci68=ci, n1=n1, n2=n2,
                                   min_detectable_3sigma=stats.min_detectable_delta(n1 + n2))

    # pile-up diagnostic: does an earlier flash in the same spill suppress later ones?
    f_sorted = f[f.batch > 0].sort_values(["spill", "t_rel"])
    f_sorted["n_before"] = f_sorted.groupby("spill").cumcount()
    res["pileup_fraction_with_earlier_flash_by_batch"] = (
        f_sorted.groupby("batch")["n_before"].apply(lambda s: float((s > 0).mean())).to_dict())
    # Masking test: for batch k>=2, compare the mean number of flashes in batch k
    # in spills WITH vs WITHOUT any selected flash in batches < k. For independent
    # Poisson arrivals (and the ~1% spill-intensity spread) the two are equal;
    # a deficit after an earlier flash measures pile-up / hit-finding loss.
    cnt = counts.pivot_table(index="spill", columns="batch", values="n", fill_value=0)
    cnt = cnt.reindex(index=trig.spill.unique(), columns=range(1, 7), fill_value=0)
    mask_rows = {}
    for k in range(2, 7):
        earlier = cnt.loc[:, 1:k - 1].sum(axis=1) > 0
        a1, a0 = cnt.loc[earlier, k], cnt.loc[~earlier, k]
        r = (a1.mean() / a0.mean()) if a0.mean() > 0 else float("nan")
        err = r * np.sqrt(1 / max(a1.sum(), 1) + 1 / max(a0.sum(), 1))
        mask_rows[k] = dict(n_with=int(earlier.sum()), n_without=int((~earlier).sum()),
                            ratio_with_over_without=float(r), err=float(err))
    res["masking_test"] = mask_rows
    # absolute rates and projections ---------------------------------------
    inwin = (f.t_rel >= edges[0]) & (f.t_rel < edges[-1])
    spot = trig.spillpot.mean()
    rock = inwin & (f.ntpc >= 4)
    cyc_spills_per_h = 3600.0 / 60.46 * 46          # 46 NuMI spills per ~60.46 s supercycle
    pot_per_spill_scan = spot * a.scan_power_kw / numi.beam_power_kw(spot)
    pot_h = pot_per_spill_scan * cyc_spills_per_h
    frs = dict(FRperSpill=float(inwin.sum() / len(trig)), FRperEthirteen=float(inwin.sum() / (spot * len(trig)) * 1e13),
               FRrockPerSpill=float(rock.sum() / len(trig)), FRspillPot=float(spot / 1e13),
               FRpowerKw=float(numi.beam_power_kw(spot)),
               FRperHourScan=float(inwin.sum() / (spot * len(trig)) * pot_h),
               FRrockPerHourScan=float(rock.sum() / (spot * len(trig)) * pot_h),
               FRrunOneTotal=float(inwin.sum() / (spot * len(trig)) * a.run1_pot),
               FRrunOneMinDelta=float(100 * stats.min_detectable_delta(inwin.sum() / (spot * len(trig)) * a.run1_pot)))
    (out / "flash_rate_summary.json").write_text(json.dumps(frs, indent=2))
    res["flash_rates"] = frs
    (out / "batch_rates.json").write_text(json.dumps(res, indent=2, default=float))

    # figures
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.hist(f.t_rel, bins=np.arange(t1 - 2000, t1 + numi.SPILL_NS + 3000, 32), histtype="step", color="#2a6fb0")
    for e_ in edges:
        ax.axvline(e_, color="#d0622b", lw=0.8)
    ax.set_xlabel("flash time relative to light beam trigger [ns]"); ax.set_ylabel("flashes / 32 ns")
    ax.set_title(f"{len(trig)} matched goodbeam spills; red: fitted batch edges", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "flash_times.pdf"); fig.savefig(out / "flash_times.png", dpi=150); plt.close(fig)
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    ax.errorbar(per.x, per.rate_e13, per.err_e13, fmt="o", color="#2a6fb0")
    for b_, row in per.iterrows():
        ax.annotate(str(b_), (row.x, row.rate_e13), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlabel("batch x at target [mm]"); ax.set_ylabel("flashes per 1e13 POT")
    fig.tight_layout(); fig.savefig(out / "rate_vs_batch_x.pdf"); plt.close(fig)
    print(json.dumps({k: res[k] for k in ("audit", "batch1_over_rest")}, indent=1, default=float))
    print(per)


if __name__ == "__main__":
    main()
