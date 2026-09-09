"""Animated GIF illustrating the batch-1 vs batches-2-6 lobe structure,
now with NOvA goodbeam quality flags, POT accounting, and delivery
efficiency.

Each 10-spill group produces TWO consecutive frames, mimicking the
arrival order within the NuMI spill train:
  frame A: the group's batch-1 readings appear (orange) — the head of
           each spill, which rides the extraction-kicker rise and lands
           ~0.3 mm to the left;
  frame B: the group's batches 2-6 follow (blue).
Points from spills that PASS all NOvA goodbeam cuts (1-5, evaluated by
get_data.cpp into the *_quality.csv) wear the full orange/blue hue;
points from FAILING spills wear a pale tint of the same hue. Every
point shown then joins a light-gray accumulation, building the full
two-lobe profile over the hour.

The main panel carries running "POT passed" / "POT missed" integrators
(spills with no usable toroid reading contribute 0 POT but are counted
as missed spills), and a third right-hand panel tracks the POT-weighted
beam-delivery efficiency: one dot per 10-spill group plus the
cumulative running efficiency. Axes are never clipped: they span the
union of the data and the NOvA goodbeam position box (|x|,|y| < 2 mm),
drawn in red. Rendered at 3000x2000.

The full catalogue of NOvA's goodbeam criteria, and exactly which of
them this pipeline applies, illustrates, or omits, is documented in the
"Beam-quality cuts" comment block of get_data.cpp.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from PIL import Image
from datetime import datetime, timedelta, timezone
import os, shutil, sys

BASE = os.environ.get("BEAMQ_DIR", os.getcwd())
# usage: python3 beam_batches_gif.py [get_data_<t0>_to_<t1>]
STEM = sys.argv[1] if len(sys.argv) > 1 else \
    "get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500"
CSV = f"{BASE}/{STEM}.csv"
QCSV = f"{BASE}/{STEM}_quality.csv"
OUT = f"{BASE}/beam_position_batches_{STEM.replace('get_data_', '')}.gif"
SPILLS_PER_GROUP = 10
CUT = 2.0                  # NOvA goodbeam position cut: |x|,|y| < 2 mm
CDT = timezone(timedelta(hours=-5))

# palette (dataviz reference, light mode, committed single look)
BLUE = "#2a78d6"      # categorical slot 1: batches 2-6, goodbeam
BLUE_F = "#a6c8ee"    # tint of slot 1: batches 2-6, fails goodbeam
ORANGE = "#eb6834"    # categorical slot 2: batch 1, goodbeam
ORANGE_F = "#f4b490"  # tint of slot 2: batch 1, fails goodbeam
RED = "#e34948"       # status "serious" — reserved for the cut boundary
INK = "#1a1a19"       # text primary
INK2 = "#6b6a60"      # text secondary
CTX = "#c9c8c0"       # accumulated context
GRID = "#e8e7e0"
SURF = "#ffffff"

print("loading CSVs...", flush=True)
df = pd.read_csv(CSV)
q = pd.read_csv(QCSV).sort_values("time_s").reset_index(drop=True)
q["pot"] = q["spillpot"].fillna(0.0)     # unknown POT counts as 0
q["good"] = q["goodbeam"].astype(bool)

def batches(device, col):
    d = df[(df.device == device) & (df["index"] >= 1) & (df.value != 0.0)]
    return d.rename(columns={"time_s": "t", "value": col, "index": "idx"})[["t", "idx", col]]

H = batches("E:HPTGT[]", "x")
V = batches("E:VPTGT[]", "y")
pieces = []
for i in range(1, 7):
    h = H[H.idx == i][["t", "x"]].sort_values("t")
    v = V[V.idx == i][["t", "y"]].sort_values("t")
    p = pd.merge_asof(h, v, on="t", tolerance=0.5, direction="nearest").dropna()
    p["idx"] = i
    pieces.append(p)
P = pd.concat(pieces).sort_values("t", kind="stable").reset_index(drop=True)
spill_times = np.sort(P["t"].unique())
P["spill"] = np.searchsorted(spill_times, P["t"])
# per-spill goodbeam flag from the quality CSV (times agree within ms)
spill_good = pd.merge_asof(
    pd.DataFrame({"t": spill_times}), q.rename(columns={"time_s": "t"})[["t", "good"]],
    on="t", tolerance=0.5, direction="nearest")["good"].fillna(False).values
print(f"spills: {len(spill_times)}, batch points: {len(P)}, "
      f"goodbeam: {spill_good.sum()}/{len(spill_times)} visualized spills", flush=True)

n_groups = len(spill_times) // SPILLS_PER_GROUP
t0 = spill_times[0]
x = P["x"].values
y = P["y"].values
tmin = (P["t"].values - t0) / 60.0
spill = P["spill"].values
is_b1 = (P["idx"] == 1).values
pt_good = spill_good[spill]

# ---- POT accounting and efficiency per group (over ALL quality spills,
# including the few with no BPM data that never appear as points) ----
group_tend = spill_times[(np.arange(n_groups) + 1) * SPILLS_PER_GROUP - 1]
qt = q["time_s"].values
qpot = q["pot"].values
qgood = q["good"].values
cum_pass = np.zeros(n_groups); cum_miss = np.zeros(n_groups)
cum_ngood = np.zeros(n_groups, int); cum_nfail = np.zeros(n_groups, int)
grp_eff = np.full(n_groups, np.nan)
prev_t = -np.inf
for g in range(n_groups):
    upto = qt <= group_tend[g] + 0.3
    ingrp = upto & (qt > prev_t + 0.3)   # prev_t starts at -inf: whole prefix
    cum_pass[g] = qpot[upto & qgood].sum()
    cum_miss[g] = qpot[upto & ~qgood].sum()
    cum_ngood[g] = (upto & qgood).sum()
    cum_nfail[g] = (upto & ~qgood).sum()
    gpot, ggood = qpot[ingrp], qgood[ingrp]
    if gpot.sum() > 0:
        grp_eff[g] = 100.0 * gpot[ggood].sum() / gpot.sum()
    prev_t = group_tend[g]
cum_eff = 100.0 * cum_pass / np.where(cum_pass + cum_miss > 0, cum_pass + cum_miss, np.nan)

def lims(c, pad=0.06):
    lo, hi = min(np.min(c), -CUT), max(np.max(c), CUT)
    p = (hi - lo) * pad
    return lo - p, hi + p
xlim, ylim = lims(x), lims(y)
tlim = (tmin[0] - 1, tmin[-1] + 1)
gmin = (group_tend - t0) / 60.0

fig = plt.figure(figsize=(15, 10), dpi=200)
fig.patch.set_facecolor(SURF)
gs = GridSpec(3, 2, width_ratios=[1.35, 1], height_ratios=[1, 1, 0.8],
              hspace=0.42, wspace=0.22, left=0.06, right=0.97, top=0.86, bottom=0.07)
axM = fig.add_subplot(gs[:, 0])
axX = fig.add_subplot(gs[0, 1])
axY = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, 1])

for ax in (axM, axX, axY, axE):
    ax.set_facecolor(SURF)
    ax.grid(color=GRID, linewidth=0.8)
    ax.tick_params(colors=INK2, labelsize=11)
    for sp in ax.spines.values():
        sp.set_color(GRID)

win0 = datetime.fromtimestamp(spill_times[0], CDT)
win1 = datetime.fromtimestamp(spill_times[-1], CDT)
fig.suptitle("NuMI spill-train arrival order and goodbeam quality at the target  (IFBeam $A9 spills)",
             fontsize=19, color=INK, x=0.06, ha="left", y=0.965)
fig.text(0.06, 0.915,
         f"{win0:%Y-%m-%d  %H:%M} → {win1:%H:%M} CDT  ·  batch 1 arrives first, batches 2-6 follow  ·  "
         "pale points fail NOvA goodbeam cuts 1-5  ·  shown points accumulate in gray",
         fontsize=12, color=INK2)

# --- main panel ---
axM.set_xlabel("Horizontal position at target [mm]", fontsize=13, color=INK)
axM.set_ylabel("Vertical position at target [mm]", fontsize=13, color=INK)
axM.set_xlim(*xlim); axM.set_ylim(*ylim)
axM.add_patch(Rectangle((-CUT, -CUT), 2 * CUT, 2 * CUT, fill=False,
                        edgecolor=RED, linewidth=2.2, linestyle=(0, (6, 3)), zorder=2))
axM.text(-CUT + 0.06, -CUT + 0.06, "NOvA goodbeam position cut (±2 mm)",
         color=RED, fontsize=11, va="bottom")
acc_sc = axM.scatter([], [], s=4, color=CTX, linewidths=0, zorder=1)
b1p_sc = axM.scatter([], [], s=46, color=ORANGE, edgecolors=SURF, linewidths=0.8, zorder=5)
b1f_sc = axM.scatter([], [], s=46, color=ORANGE_F, edgecolors=SURF, linewidths=0.8, zorder=5)
b26p_sc = axM.scatter([], [], s=42, color=BLUE, edgecolors=SURF, linewidths=0.8, zorder=4)
b26f_sc = axM.scatter([], [], s=42, color=BLUE_F, edgecolors=SURF, linewidths=0.8, zorder=4)
TXTBOX = dict(facecolor=SURF, edgecolor="none", pad=1.5)
stamp = axM.text(0.03, 0.985, "", transform=axM.transAxes, fontsize=13,
                 color=INK, va="top", family="monospace", bbox=TXTBOX)
phase_txt = axM.text(0.03, 0.955, "", transform=axM.transAxes, fontsize=12,
                     color=INK2, va="top", family="monospace", bbox=TXTBOX)
pot_txt = axM.text(0.03, 0.925, "", transform=axM.transAxes, fontsize=12,
                   color=INK, va="top", family="monospace", bbox=TXTBOX)
axM.legend(handles=[
    Line2D([], [], marker="o", ls="", color=ORANGE, ms=9, label="Batch 1, goodbeam"),
    Line2D([], [], marker="o", ls="", color=ORANGE_F, ms=9, label="Batch 1, fails cuts"),
    Line2D([], [], marker="o", ls="", color=BLUE, ms=9, label="Batches 2-6, goodbeam"),
    Line2D([], [], marker="o", ls="", color=BLUE_F, ms=9, label="Batches 2-6, fails cuts"),
    Line2D([], [], marker="o", ls="", color=CTX, ms=7, label="Accumulated profile")],
    loc="lower right", frameon=True, facecolor=SURF, edgecolor="none",
    framealpha=1.0, fontsize=11, labelcolor=INK)

# --- time-series panels: every batch reading, no averaging ---
for ax, lab in ((axX, "Horizontal"), (axY, "Vertical")):
    ax.set_xlim(*tlim)
    ax.axhline(+CUT, color=RED, lw=1.4, linestyle=(0, (6, 3)))
    ax.axhline(-CUT, color=RED, lw=1.4, linestyle=(0, (6, 3)))
    ax.set_title(f"{lab} position vs time (per batch)", fontsize=13, color=INK, loc="left")
    ax.set_ylabel("[mm]", fontsize=12, color=INK)
axX.set_ylim(*lims(x))
axY.set_ylim(*lims(y))
accX = axX.scatter([], [], s=1.5, color=CTX, linewidths=0)
accY = axY.scatter([], [], s=1.5, color=CTX, linewidths=0)
curBX = axX.scatter([], [], s=6, color=BLUE, linewidths=0)
curBY = axY.scatter([], [], s=6, color=BLUE, linewidths=0)
curOX = axX.scatter([], [], s=6, color=ORANGE, linewidths=0)
curOY = axY.scatter([], [], s=6, color=ORANGE, linewidths=0)
curX = axX.axvline(0, color=INK2, lw=1.2, alpha=0.7)
curY = axY.axvline(0, color=INK2, lw=1.2, alpha=0.7)

# --- efficiency panel ---
axE.set_xlim(*tlim)
axE.set_ylim(0, 108)
axE.set_title("POT-weighted delivery efficiency (goodbeam)",
              fontsize=13, color=INK, loc="left")
axE.set_ylabel("[%]", fontsize=12, color=INK)
axE.set_xlabel(f"Minutes after {win0:%H:%M} CDT", fontsize=12, color=INK)
eff_dots = axE.scatter([], [], s=14, color=BLUE, linewidths=0, zorder=3)
eff_line, = axE.plot([], [], color=INK2, lw=1.6, zorder=2)
curE = axE.axvline(0, color=INK2, lw=1.2, alpha=0.7)
axE.legend(handles=[
    Line2D([], [], marker="o", ls="", color=BLUE, ms=7, label="per 10 spills"),
    Line2D([], [], color=INK2, lw=1.6, label="cumulative")],
    loc="lower right", frameon=True, facecolor=SURF, edgecolor="none",
    framealpha=1.0, fontsize=10, labelcolor=INK)

EMPTY = np.empty((0, 2))
def offs(mask, cx, cy):
    return np.c_[cx[mask], cy[mask]] if mask.any() else EMPTY

# Frames are spooled to disk as PNGs and reloaded lazily at assembly
# time, so long windows (1000+ frames) don't exhaust RAM.
FRAMES_DIR = f"{BASE}/_frames_{STEM.replace('get_data_', '')}"
shutil.rmtree(FRAMES_DIR, ignore_errors=True)
os.makedirs(FRAMES_DIR)
frame_files = []
for g in range(n_groups):
    a, b = g * SPILLS_PER_GROUP, (g + 1) * SPILLS_PER_GROUP
    past = spill < a
    grp = (spill >= a) & (spill < b)
    g_b1, g_b26 = grp & is_b1, grp & ~is_b1
    ts0 = datetime.fromtimestamp(spill_times[a], CDT).strftime("%H:%M:%S")
    ts1 = datetime.fromtimestamp(spill_times[b - 1], CDT).strftime("%H:%M:%S")
    tcur = (spill_times[b - 1] - t0) / 60.0

    acc_sc.set_offsets(offs(past, x, y))
    accX.set_offsets(offs(past, tmin, x))
    accY.set_offsets(offs(past, tmin, y))
    curX.set_xdata([tcur]); curY.set_xdata([tcur]); curE.set_xdata([tcur])
    stamp.set_text(f"{ts0}-{ts1} CDT   spills {a}-{b - 1}")
    ge = f"{grp_eff[g]:5.1f}%" if np.isfinite(grp_eff[g]) else "  n/a"
    pot_txt.set_text(
        f"POT passed {cum_pass[g]:.4g} ({cum_ngood[g]} spills)   "
        f"POT missed {cum_miss[g]:.3g} ({cum_nfail[g]} spills)\n"
        f"efficiency: this group {ge}   cumulative {cum_eff[g]:6.2f}%")
    keep = np.isfinite(grp_eff[: g + 1])
    eff_dots.set_offsets(np.c_[gmin[: g + 1][keep], grp_eff[: g + 1][keep]])
    eff_line.set_data(gmin[: g + 1], cum_eff[: g + 1])

    for phase in ("A", "B"):
        if phase == "A":     # heads of the spill trains arrive first
            b1p_sc.set_offsets(offs(g_b1 & pt_good, x, y))
            b1f_sc.set_offsets(offs(g_b1 & ~pt_good, x, y))
            curOX.set_offsets(offs(g_b1, tmin, x))
            curOY.set_offsets(offs(g_b1, tmin, y))
            b26p_sc.set_offsets(EMPTY); b26f_sc.set_offsets(EMPTY)
            curBX.set_offsets(EMPTY); curBY.set_offsets(EMPTY)
            phase_txt.set_text("batch 1 arriving ...")
        else:                # ... followed by the rest of each train
            b26p_sc.set_offsets(offs(g_b26 & pt_good, x, y))
            b26f_sc.set_offsets(offs(g_b26 & ~pt_good, x, y))
            curBX.set_offsets(offs(g_b26, tmin, x))
            curBY.set_offsets(offs(g_b26, tmin, y))
            phase_txt.set_text("batches 2-6 arriving")
        fig.canvas.draw()
        img = Image.frombuffer("RGBA", fig.canvas.get_width_height(),
                               fig.canvas.buffer_rgba()).convert("RGB")
        fname = f"{FRAMES_DIR}/f{len(frame_files):05d}.png"
        img.quantize(colors=128, method=Image.MEDIANCUT).save(fname)
        frame_files.append(fname)
    if g % 25 == 0:
        print(f"group {g}/{n_groups}", flush=True)

print("assembling GIF...", flush=True)
first = Image.open(frame_files[0])
rest = (Image.open(f) for f in frame_files[1:])
first.save(OUT, save_all=True, append_images=rest,
           duration=70, loop=0, optimize=True)
size = first.size
shutil.rmtree(FRAMES_DIR)
print(f"done: {OUT}  {os.path.getsize(OUT)/1e6:.1f} MB, "
      f"{len(frame_files)} frames, {size[0]}x{size[1]}", flush=True)
