#!/usr/bin/env python3
"""Statistical feasibility of using 2x2 / Mx2 rates to locate the beam optimum.

Produces (in --out, default ../results/sensitivity):
  intercept_response.pdf   relative rate vs beam offset (interception model)
  batch_split_reach.pdf    smallest batch-1 vs batches-2-6 difference visible
                           at 3 sigma vs total counts, with model predictions
  scan_precision.pdf       toy-MC precision on the optimum vs counts/point for
                           (a) a plateau scan (+-1 mm, parabola) and
                           (b) an edge-reaching scan (+-4.5 mm, template fit)
  sensitivity.json         every number quoted in the report

All detector rates are INPUT ASSUMPTIONS (command-line flags) until replaced by
measured Run 1 rates; the JSON records which values were used.

Example:
  python3 scripts/sensitivity_study.py --sigma 1.43 --width 7.4 \
      --rate-at-100kw 50 150 500 2000
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import target_model as tm, stats as st  # noqa: E402

INK, MUTED, C1, C2, C3, C4 = "#222222", "#777777", "#2a6fb0", "#d0622b", "#3a9a5b", "#8a4fb0"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sigma", type=float, default=1.43, help="beam RMS at target [mm] (Run 1 multiwire ~1.43)")
    p.add_argument("--width", type=float, default=7.4, help="fin full width [mm] (NOvA-era ME fin: 7.4)")
    p.add_argument("--eps-miss", type=float, default=0.2, help="relative nu yield of protons missing the fin")
    p.add_argument("--batch-offset", type=float, default=0.30, help="batch1 - batches2-6 offset at target [mm]")
    p.add_argument("--rate-at-100kw", type=float, nargs="+", default=[50, 150, 500, 2000],
                   help="counts/hour at 100 kW with the beam at the optimum (assumed)")
    p.add_argument("--scan-points", type=float, nargs="+", default=list(np.arange(-4.5, 4.51, 1.5)))
    p.add_argument("--toys", type=int, default=1000)
    p.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "results" / "sensitivity"))
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    R = lambda x: tm.relative_rate(x, a.sigma, a.width, eps_miss=a.eps_miss)
    res = dict(inputs=vars(a).copy())

    # 1. response curve ---------------------------------------------------
    x = np.linspace(-7, 7, 701)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    for s_, c in ((1.0, C3), (a.sigma, C1), (2.0, C2)):
        ax.plot(x, tm.relative_rate(x, s_, a.width, eps_miss=a.eps_miss), color=c, lw=2,
                label=fr"$\sigma$ = {s_:.2f} mm")
    ax.axvspan(-a.width / 2, a.width / 2, color="#000000", alpha=0.05, lw=0)
    ax.set_xlabel("beam centroid relative to fin centre [mm]")
    ax.set_ylabel("relative rate (interception model)")
    ax.set_ylim(0, 1.05); ax.legend(frameon=False, fontsize=8)
    ax.set_title(f"fin width {a.width} mm, miss-yield {a.eps_miss}", fontsize=9, color=MUTED)
    fig.tight_layout(); fig.savefig(out / "intercept_response.pdf"); plt.close(fig)
    res["curvature_per_mm2"] = float(tm.curvature(a.sigma, a.width))
    res["R_at"] = {f"{v:.1f}": float(R(v)) for v in (0.5, 1, 2, 3, 3.7, 4.5)}

    # 2. batch split reach ------------------------------------------------
    d = a.batch_offset
    centres = np.array([0.0, 0.5, 1.0, 1.5])
    pred = {f"{c:.1f}": float(abs(tm.two_position_difference(c - d / 2, c + d / 2, a.sigma, a.width,
                                                              eps_miss=a.eps_miss))) for c in centres}
    res["batch_split_predicted_delta"] = pred
    N = np.logspace(3, 8, 200)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.loglog(N, [st.min_detectable_delta(n) for n in N], color=C1, lw=2, label=r"3$\sigma$ reach, stat only")
    ax.loglog(N, [st.min_detectable_delta(n, pot_syst_frac=0.003) for n in N], color=C2, lw=2,
              ls="--", label="+ 0.3% per-batch POT-share syst.")
    for (c, v), col in zip(pred.items(), (C3, C4, MUTED, INK)):
        if v > 0:
            ax.axhline(v, color=col, lw=1, ls=":")
            ax.text(1.2e3, v * 1.08, f"model, mean offset {c} mm from fin centre", fontsize=7, color=col)
    ax.axvspan(3e4, 3e5, color=C1, alpha=0.08, lw=0)
    ax.text(3.3e4, 0.3, "Run 1 counts\n(order of mag.)", fontsize=7, color=C1)
    ax.set_xlabel("total in-spill counts (all six batches)")
    ax.set_ylabel(r"min. |$\Delta$rate/rate| at 3$\sigma$")
    ax.set_ylim(1e-4, 1); ax.legend(frameon=False, fontsize=8, loc="upper right")
    fig.tight_layout(); fig.savefig(out / "batch_split_reach.pdf"); plt.close(fig)
    res["N_needed_3sigma"] = {k: float(st.required_total_counts(v)) if v > 0 else None for k, v in pred.items()}
    res["min_delta_3sigma_at"] = {f"{n:.0e}": float(st.min_detectable_delta(n)) for n in (3e4, 1e5, 3e5, 1e6)}

    # 3. scan precision ---------------------------------------------------
    counts = np.array([10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4])
    wide, narrow = [], []
    for n in counts:
        wide.append(st.simulate_scan(a.scan_points, n, R, 0.3, n_toys=a.toys, fit="template", rng=7))
        narrow.append(st.simulate_scan([-1.0, 0.0, 1.0], n, R, 0.3, n_toys=a.toys, fit="parabola", rng=7))
    # robustness: fit with a WRONG template (sigma, eps_miss off) - bias check
    Rwrong = lambda x: tm.relative_rate(x, a.sigma * 1.3, a.width, eps_miss=0.0)
    wrong = [st.simulate_scan(a.scan_points, n, R, 0.3, n_toys=a.toys // 2, fit="template", rng=9)
             for n in (200, 2000)]
    res["scan_points"] = list(map(float, a.scan_points))
    res["scan_wide_template"] = {str(int(n)): w for n, w in zip(counts, wide)}
    res["scan_narrow_parabola"] = {str(int(n)): w for n, w in zip(counts, narrow)}
    # mis-modelled template: generate with R, fit with Rwrong
    rng = np.random.default_rng(3); xs = np.asarray(a.scan_points); biases = {}
    for n in (200, 2000):
        est = []
        for _ in range(300):
            k = rng.poisson(n * R(xs - 0.3))
            grid = np.arange(-3.7, 4.3, 0.01)
            sh = Rwrong(xs[None, :] - grid[:, None]); A = k.sum() / sh.sum(1)
            lam = np.maximum(A[:, None] * sh, 1e-12)
            est.append(grid[np.argmax((k * np.log(lam) - lam).sum(1))])
        biases[str(n)] = dict(bias=float(np.mean(est) - 0.3), rms=float(np.std(est)))
    res["scan_wide_wrong_template"] = biases

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.loglog(counts, [w["rms"] for w in wide], "o-", color=C1, lw=2,
              label=f"edge-reaching scan, {len(a.scan_points)} pts $\\pm${max(a.scan_points):.1f} mm, template fit")
    ax.loglog(counts, [w["rms"] for w in narrow], "s-", color=C2, lw=2,
              label=r"plateau scan, 3 pts $\pm$1 mm, parabola fit")
    for n, w in zip(counts, narrow):
        if w["frac_fail"] > 0.05:
            ax.annotate(f"{100*w['frac_fail']:.0f}% no max", (n, w["rms"]), fontsize=6, color=C2,
                        xytext=(2, 4), textcoords="offset points")
    ax.axhline(0.26, color=MUTED, ls=":", lw=1); ax.text(11, 0.28, "BPM-based extrapolation (2005 NuMI BBA)", fontsize=7, color=MUTED)
    ax.axhline(0.03, color=MUTED, ls=":", lw=1); ax.text(11, 0.033, "SEM-based extrapolation (2005 NuMI BBA)", fontsize=7, color=MUTED)
    for r, col in zip(a.rate_at_100kw, (C3, C4, INK, MUTED)):
        ax.axvline(r, color=col, lw=1, ls="--")
        ax.text(r * 1.05, 2.0, f"{r:g}/h", fontsize=7, color=col, rotation=90, va="top")
    ax.set_xlabel("counts per scan point with beam at optimum (= rate x hours at that point)")
    ax.set_ylabel("RMS of fitted optimum [mm]")
    ax.set_ylim(0.01, 3); ax.legend(frameon=False, fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(out / "scan_precision.pdf"); plt.close(fig)

    # time budget table: hours per point for 0.5 / 0.25 mm with each assumed rate
    target_rms = (0.5, 0.25, 0.1)
    rms = np.array([w["rms"] for w in wide])
    need = {}
    for tr in target_rms:
        n_need = float(np.exp(np.interp(np.log(tr), np.log(rms[::-1]), np.log(counts[::-1]))))
        need[str(tr)] = dict(counts_per_point=n_need,
                             hours_per_point={f"{r:g}": n_need / r for r in a.rate_at_100kw})
    res["time_budget_100kW"] = need
    (out / "sensitivity.json").write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
