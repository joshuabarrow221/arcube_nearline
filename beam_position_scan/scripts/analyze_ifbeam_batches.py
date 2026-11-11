#!/usr/bin/env python3
"""Step 1 of the proof of concept: is batch 1 really somewhere else on the target?

Reads one or more beam_quality/get_data.cpp output stems, builds per-batch
positions at the target z, and reports:
  * per-batch mean/RMS position and intensity share (goodbeam spills only)
  * per-spill (batch 1) - (batches 2-6, intensity weighted) offset: mean, RMS,
    and its stability over time (10-minute blocks)
  * the "lever arm" available to each proposed sub-division:
        A  batch 1 vs batches 2-6
        B  halves: each population split at its own median x (or y)
        C  quarters: x/y quadrants about each population's median
    i.e. the mean position difference between the sub-samples that would be
    compared, plus the POT share of each
  * horn-current sign and E:HRNDIR summary (polarity sanity check)
  * the predicted fractional rate difference for each comparison from the
    interception model, for an assumed fin width and the measured width

Outputs into --out: ifbeam_batches.json, batch_positions.pdf,
batch1_offset_hist.pdf, batch1_offset_vs_time.pdf

Usage:
  python3 scripts/analyze_ifbeam_batches.py STEM [STEM ...] --out results/ifbeam
  (STEM = path without .csv, e.g. .../get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500)
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
from matplotlib.patches import Ellipse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import ifbeam, target_model as tm  # noqa: E402

COL = ["#d0622b", "#2a6fb0", "#3a9a5b", "#8a4fb0", "#b08a2a", "#2aa0a0"]


def ellipse(ax, x, y, color, nsig=1.0, **kw):
    c = np.cov(np.vstack([x, y]))
    vals, vecs = np.linalg.eigh(c)
    ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
    ax.add_patch(Ellipse((x.mean(), y.mean()), 2 * nsig * np.sqrt(vals[1]), 2 * nsig * np.sqrt(vals[0]),
                         angle=ang, fill=False, color=color, lw=1.5, **kw))


def lever(df, mask_a, mask_b):
    """Mean x,y difference (a - b), POT-weighted, and POT shares."""
    pa, pb = df.loc[mask_a, "pot_batch"], df.loc[mask_b, "pot_batch"]
    xa = np.average(df.loc[mask_a, "x_tgt"], weights=pa); xb = np.average(df.loc[mask_b, "x_tgt"], weights=pb)
    ya = np.average(df.loc[mask_a, "y_tgt"], weights=pa); yb = np.average(df.loc[mask_b, "y_tgt"], weights=pb)
    tot = pa.sum() + pb.sum()
    return dict(dx=float(xa - xb), dy=float(ya - yb), xa=float(xa), xb=float(xb), ya=float(ya), yb=float(yb),
                pot_share_a=float(pa.sum() / tot), pot_share_b=float(pb.sum() / tot))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("stems", nargs="+")
    p.add_argument("--out", default="results/ifbeam")
    p.add_argument("--width", type=float, default=7.4)
    p.add_argument("--eps-miss", type=float, default=0.2)
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    frames = []
    for k, s in enumerate(a.stems):
        b = ifbeam.load_batches(s)
        b["window"] = Path(s).name
        b["spill"] = b["spill"] + 10_000_000 * k        # unique across windows
        frames.append(b)
    b = pd.concat(frames, ignore_index=True)
    good = (b["goodbeam"] == 1) & b["x_tgt"].notna() & b["y_tgt"].notna() & (b["pot_batch"] > 0)
    g = b[good].copy()
    res = dict(windows=[Path(s).name for s in a.stems],
               n_spills_total=int(b["spill"].nunique()), n_spills_good=int(g["spill"].nunique()),
               pot_good=float(g.groupby("spill")["spillpot"].first().sum()))

    summ = ifbeam.batch_summary(g)
    res["per_batch"] = summ.reset_index().to_dict(orient="records")
    d = ifbeam.batch1_minus_rest(g)
    res["batch1_minus_rest"] = dict(dx_mean=float(d.dx.mean()), dx_rms=float(d.dx.std()),
                                    dx_sem=float(d.dx.std() / np.sqrt(len(d))),
                                    dy_mean=float(d.dy.mean()), dy_rms=float(d.dy.std()),
                                    dy_sem=float(d.dy.std() / np.sqrt(len(d))), n=int(len(d)))
    # spill-level widths from the multiwire (quality CSV) and polarity
    sp = g.groupby("spill").first()
    res["widths"] = dict(sigx_median=float(sp.widthx.median()), sigy_median=float(sp.widthy.median()))
    res["spill_pos"] = dict(posx_mean=float(sp.posx.mean()), posx_rms=float(sp.posx.std()),
                            posy_mean=float(sp.posy.mean()), posy_rms=float(sp.posy.std()))
    res["hornI"] = dict(median=float(b.groupby("spill")["hornI"].first().median()))

    # time stability of the batch-1 offset (10-minute blocks)
    t = g[g.batch == 1].set_index("spill")["time_s"].reindex(d.index)
    blk = ((t - t.min()) // 600).astype(int)
    stab = d.groupby(blk).agg(dx=("dx", "mean"), dy=("dy", "mean"), n=("dx", "size"))
    stab["t_min"] = t.groupby(blk).min().values
    res["batch1_offset_blocks"] = stab.reset_index(drop=True).to_dict(orient="records")

    # lever arms ------------------------------------------------------------
    lev = {}
    is1 = g.batch == 1
    lev["A_batch1_vs_2to6"] = lever(g, is1, ~is1)
    for name, sel in (("batch1", is1), ("batches2to6", ~is1)):
        sub = g[sel]
        mx, my = sub.x_tgt.median(), sub.y_tgt.median()
        lev[f"B_{name}_x_halves"] = lever(g, sel & (g.x_tgt < mx), sel & (g.x_tgt >= mx))
        lev[f"B_{name}_y_halves"] = lever(g, sel & (g.y_tgt < my), sel & (g.y_tgt >= my))
        q = {}
        for qx in (0, 1):
            for qy in (0, 1):
                m = sel & ((g.x_tgt >= mx) == bool(qx)) & ((g.y_tgt >= my) == bool(qy))
                q[f"x{'+' if qx else '-'}y{'+' if qy else '-'}"] = dict(
                    x=float(np.average(g.x_tgt[m], weights=g.pot_batch[m])),
                    y=float(np.average(g.y_tgt[m], weights=g.pot_batch[m])),
                    pot_share=float(g.pot_batch[m].sum() / g.pot_batch.sum()))
        lev[f"C_{name}_quadrants"] = q
    # 'extreme' halves: batch1 left half vs batches2-6 right half (max lever arm)
    lev["D_batch1_left_vs_2to6_right"] = lever(
        g, is1 & (g.x_tgt < g[is1].x_tgt.median()), (~is1) & (g.x_tgt >= g[~is1].x_tgt.median()))
    res["lever_arms"] = lev

    # model predictions for each lever arm -------------------------------------
    sig = res["widths"]["sigx_median"]
    pred = {}
    for k, v in lev.items():
        if "dx" not in v:
            continue
        pred[k] = {f"fin_centre_at_{c:+.1f}mm": float(tm.two_position_difference(
            v["xa"] - c, v["xb"] - c, sig, a.width, eps_miss=a.eps_miss)) for c in (0.0, -0.5, -1.0, 0.5, 1.0)}
    res["predicted_fractional_rate_difference"] = pred
    (out / "ifbeam_batches.json").write_text(json.dumps(res, indent=2, default=float))

    # figures ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for bt in range(1, 7):
        s = g[g.batch == bt]
        n = min(len(s), 3000)
        s = s.sample(n, random_state=1)
        ax.scatter(s.x_tgt, s.y_tgt, s=2, color=COL[bt - 1], alpha=0.25, lw=0)
        ellipse(ax, s.x_tgt.values, s.y_tgt.values, COL[bt - 1])
        ax.plot([], [], "o", color=COL[bt - 1], label=f"batch {bt}")
    ax.axvline(0, color="#999", lw=0.5); ax.axhline(0, color="#999", lw=0.5)
    ax.set_xlabel("x at target z [mm]"); ax.set_ylabel("y at target z [mm]")
    ax.legend(frameon=False, fontsize=7, markerscale=1.2)
    ax.set_title("per-batch BPM positions extrapolated to target (goodbeam)", fontsize=8)
    fig.tight_layout(); fig.savefig(out / "batch_positions.pdf"); fig.savefig(out / "batch_positions.png", dpi=150); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(7, 3))
    for ax_, col, lab in ((axs[0], "dx", "x"), (axs[1], "dy", "y")):
        ax_.hist(d[col], bins=120, color=COL[1], histtype="stepfilled", alpha=0.7)
        ax_.axvline(d[col].mean(), color=COL[0])
        ax_.set_xlabel(f"batch1 - batches2-6, {lab} at target [mm]")
    axs[0].set_ylabel("spills")
    fig.tight_layout(); fig.savefig(out / "batch1_offset_hist.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot((t - t.min()) / 60, d.dx, ",", color=COL[1], alpha=0.5)
    ax.plot((stab.t_min - t.min()) / 60 + 5, stab.dx, "o-", color=COL[0], ms=3)
    ax.set_xlabel("minutes since first spill (windows concatenated)"); ax.set_ylabel("batch1 - rest, x [mm]")
    fig.tight_layout(); fig.savefig(out / "batch1_offset_vs_time.pdf"); plt.close(fig)
    print(json.dumps({k: res[k] for k in ("n_spills_good", "batch1_minus_rest", "widths", "spill_pos", "hornI")}, indent=1))


if __name__ == "__main__":
    main()
