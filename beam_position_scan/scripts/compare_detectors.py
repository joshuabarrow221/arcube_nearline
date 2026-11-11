#!/usr/bin/env python3
"""Do 2x2 light flashes and MINERvA (Mx2) track-like slices track each other?

The hypothesis (J. Barrow): Mx2 track-like rates should be highly correlated
with 2x2 neutrino activity, so the bigger detector can stand in for the
smaller one in a beam-position scan. What this script can and cannot test:

  * Both rates scale with POT, so their *integrated* rates trivially track
    each other. That is not evidence for anything.
  * Spill-by-spill, each detector sees O(1) objects; Poisson noise swamps any
    common beam fluctuation. The spill-level correlation coefficient is
    reported with its Poisson-only expectation so nobody over-reads it.
  * Measured: the spill-level correlation is ~0.4, thousands of times the
    POT-only expectation, because many objects are the SAME particle seen by
    both detectors (see coincidence_2x2_mx2 in the JSON and the report).
  * The useful test is STABILITY OF THE RATIO: in blocks of --block spills
    (~1 min), R = N_mx2 / N_2x2. If the two respond identically to whatever
    the beam did (intensity, position jitter, horn), R is constant and the
    block chi2/ndf is ~1. Excess scatter means they respond differently or
    one has a varying efficiency (dead time, DAQ). With Run 1's tiny position
    jitter this cannot test *position* response; only a scan can.

Inputs: 2x2 *.flashes.npz, Mx2 *.mx2slices.npz, IFBeam stems. Per detector the
in-spill window is found with detector_rates.fit_batch_edges ('auto'); the
selection for each is a pandas query (--sel-2x2, --sel-mx2).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import ifbeam, numi, detector_rates as dr  # noqa: E402


def per_spill_counts(paths, spills, sel, max_resid_ms=10.0):
    """Return Series spill -> in-spill object count (only cleanly matched spills)."""
    trig, fl, off = [], [], 0
    for pth in paths:
        d = np.load(pth)
        b = d["trig_type"] == 1
        idx = np.where(b)[0]
        trig.append(pd.DataFrame({"row": idx + off, "utime_s": d["trig_utime_ms"][b] / 1e3}))
        f = pd.DataFrame({c[6:]: d[c] for c in d.files if c.startswith("flash_")})
        f["row"] = f["event"] + off
        fl.append(f)
        off += len(d["trig_type"])
    trig, fl = pd.concat(trig, ignore_index=True), pd.concat(fl, ignore_index=True)
    j = dr.match_spills(trig.utime_s.values, spills.time_s.values, 0.5)
    trig["spill"] = np.where(j >= 0, spills.spill.values[np.maximum(j, 0)], -1)
    res = trig.utime_s.values - spills.time_s.values[np.maximum(j, 0)]
    res = np.where(j >= 0, res, np.nan)
    res = (res - np.nanmedian(res)) * 1e3
    ok = (j >= 0) & (np.abs(res) <= max_resid_ms)
    trig = trig[ok]
    trig = trig[~trig.spill.duplicated(keep=False)]
    f = fl[fl.row.isin(trig.row)].query(sel)
    t1, _ = dr.fit_batch_edges(f.t_rel.values, float(np.median(f.t_rel)) - numi.SPILL_NS / 2,
                               search_ns=float(np.ptp(f.t_rel)), step_ns=8.0)
    t1, edges = dr.fit_batch_edges(f.t_rel.values, t1, search_ns=600)
    f = f[(f.t_rel >= edges[0]) & (f.t_rel < edges[-1])]
    n = f.merge(trig[["row", "spill"]], on="row").groupby("spill").size()
    return n.reindex(trig.spill.values, fill_value=0), float(t1)


def _objects(paths, spills, sel, t1):
    tr, fl, off = [], [], 0
    for pth in paths:
        d = np.load(pth)
        b = d["trig_type"] == 1
        tr.append(pd.DataFrame({"row": np.where(b)[0] + off, "u": d["trig_utime_ms"][b] / 1e3}))
        f = pd.DataFrame({c[6:]: d[c] for c in d.files if c.startswith("flash_")})
        f["row"] = f["event"] + off
        fl.append(f)
        off += len(d["trig_type"])
    tr, fl = pd.concat(tr), pd.concat(fl).query(sel)
    j = dr.match_spills(tr.u.values, spills.time_s.values, 0.5)
    tr["spill"] = np.where(j >= 0, spills.spill.values[np.maximum(j, 0)], -1)
    tr = tr[(tr.spill >= 0) & ~tr.spill.duplicated(keep=False)]
    f = fl.merge(tr[["row", "spill"]], on="row")
    f["t"] = f.t_rel - t1
    return f[(f.t > -200) & (f.t < numi.SPILL_NS + 200)]


def coincidences(a, spills, t1a, t1b, window_ns=60.0):
    """Same-spill 2x2-flash / Mx2-object time differences (batch-1 starts aligned).

    Returns the peak position, the fraction of 2x2 flashes with an Mx2 object
    within +-window_ns of the peak, the reverse fraction, and the accidental
    level estimated from the flat |dt| in 300-900 ns sidebands.
    """
    A = _objects(a.flashes, spills, a.sel_2x2, t1a)
    M = _objects(a.mx2, spills, a.sel_mx2, t1b)
    pr = A[["spill", "t"]].reset_index(names="ia").merge(
        M[["spill", "t"]].reset_index(names="im"), on="spill", suffixes=("_a", "_m"))
    dt = (pr.t_a - pr.t_m).values
    pk = float(np.median(dt[np.abs(dt) < 150]))
    inwin = np.abs(dt - pk) < window_ns
    side = (np.abs(dt - pk) > 300) & (np.abs(dt - pk) < 900)
    acc = side.sum() / 1200.0 * 2 * window_ns
    na = pr[inwin].drop_duplicates(["spill", "t_a"]).shape[0]
    nm = pr[inwin].drop_duplicates(["spill", "t_m"]).shape[0]
    return dict(dt_peak_ns=pk, window_ns=window_ns, pairs_in_window=int(inwin.sum()),
                accidental_pairs_est=float(acc), n_2x2=int(len(A)), n_mx2=int(len(M)),
                frac_2x2_with_mx2=na / max(len(A), 1), frac_mx2_with_2x2=nm / max(len(M), 1),
                dt_rms_ns=float(np.std(dt[inwin])))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--flashes", nargs="+", required=True)
    p.add_argument("--mx2", nargs="+", required=True)
    p.add_argument("--ifbeam", nargs="+", required=True)
    p.add_argument("--sel-2x2", default="nhits >= 3")
    p.add_argument("--sel-mx2", default="ntrk >= 1")
    p.add_argument("--block", type=int, default=50, help="spills per block (~1 min at 46 spills/min)")
    p.add_argument("--out", default="results/compare_detectors")
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    frames = []
    for k, s in enumerate(a.ifbeam):
        b = ifbeam.load_batches(s); b["spill"] += 10_000_000 * k; frames.append(b)
    sp = pd.concat(frames).groupby("spill").first().reset_index()[["spill", "time_s", "goodbeam", "spillpot"]]
    na, t1a = per_spill_counts(a.flashes, sp, a.sel_2x2)
    nb, t1b = per_spill_counts(a.mx2, sp, a.sel_mx2)
    both = pd.DataFrame({"n2x2": na, "nmx2": nb}).dropna().astype(int).reset_index(names="spill")
    both = both.merge(sp, on="spill").query("goodbeam == 1").sort_values("time_s").reset_index(drop=True)
    pot = both.spillpot.values
    res = dict(args=vars(a), spills_common=int(len(both)), t1_2x2_ns=t1a, t1_mx2_ns=t1b,
               per_e13_2x2=float(both.n2x2.sum() / pot.sum() * 1e13),
               per_e13_mx2=float(both.nmx2.sum() / pot.sum() * 1e13),
               mx2_over_2x2=float(both.nmx2.sum() / max(both.n2x2.sum(), 1)))
    # spill-level correlation and its Poisson-only expectation (common POT only)
    r = float(np.corrcoef(both.n2x2, both.nmx2)[0, 1])
    la, lb = both.n2x2.mean(), both.nmx2.mean()
    vpot = np.var(pot / pot.mean())
    r_exp = la * lb * vpot / np.sqrt((la + la**2 * vpot) * (lb + lb**2 * vpot))
    res["spill_corr"] = r
    res["spill_corr_expected_from_pot_only"] = float(r_exp)
    res["spill_corr_err"] = float(1 / np.sqrt(len(both)))
    # block ratio stability
    blk = np.arange(len(both)) // a.block
    g = both.groupby(blk).agg(n2x2=("n2x2", "sum"), nmx2=("nmx2", "sum"), pot=("spillpot", "sum"),
                              t=("time_s", "mean"), n=("spill", "size"))
    g = g[(g.n == a.block) & (g.n2x2 > 0)]
    R = res["mx2_over_2x2"]
    # Poisson: N_mx2 | (N_mx2+N_2x2) ~ Binomial(f = R/(1+R)) per block
    # The two counts are NOT independent: a through-going (rock) muon makes a
    # 2x2 flash and an Mx2 track in the same spill (measured: >50% of 2x2
    # flashes have an Mx2 track within ~60 ns). So the block variance of
    # N_mx2 - R*N_2x2 uses the empirical per-spill covariance, not Poisson.
    C = np.cov(both.nmx2, both.n2x2)
    var_spill = C[0, 0] + R**2 * C[1, 1] - 2 * R * C[0, 1]
    res["per_spill_cov"] = C.tolist()
    chi2 = float((((g.nmx2 - R * g.n2x2) ** 2) / (g.n * var_spill)).sum())
    res["blocks"] = int(len(g))
    res["block_ratio_chi2"] = chi2
    res["block_ratio_ndf"] = int(len(g) - 1)
    from scipy import stats as st
    res["block_ratio_pvalue"] = float(st.chi2.sf(chi2, len(g) - 1))
    # per-POT rate of each vs time (drift check) as chi2 vs constant
    for c in ("n2x2", "nmx2"):
        rate = g[c].sum() / g.pot.sum()
        res[f"{c}_block_rate_chi2_ndf"] = float((((g[c] - rate * g.pot) ** 2) / (rate * g.pot)).sum() / (len(g) - 1))
    res["coincidence_2x2_mx2"] = coincidences(a, sp, t1a, t1b)
    (out / "compare_detectors.json").write_text(json.dumps(res, indent=2, default=float))
    g.to_csv(out / "blocks.csv")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(6.4, 5), sharex=True)
    tt = (g.t - g.t.min()) / 3600
    for c, col, lab in (("n2x2", "#2a6fb0", "2x2 flashes"), ("nmx2", "#d0622b", "Mx2 track slices")):
        ax[0].errorbar(tt, g[c] / g.pot * 1e13, np.sqrt(g[c]) / g.pot * 1e13, fmt=".", color=col, label=lab)
    ax[0].set_ylabel("per 1e13 POT"); ax[0].legend(fontsize=7)
    rb = g.nmx2 / g.n2x2
    ax[1].errorbar(tt, rb, rb * np.sqrt(1 / g.nmx2.clip(1) + 1 / g.n2x2), fmt=".", color="k")
    ax[1].axhline(R, color="#d0622b", lw=0.8)
    ax[1].set_ylabel("Mx2 / 2x2"); ax[1].set_xlabel("hours since first block (gaps removed)")
    ax[0].set_title(f"{a.block}-spill blocks; ratio chi2/ndf = {chi2:.0f}/{len(g)-1}", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "rates_vs_time.pdf"); fig.savefig(out / "rates_vs_time.png", dpi=150)
    print(json.dumps({k: v for k, v in res.items() if k != "args"}, indent=1, default=float))


if __name__ == "__main__":
    main()
