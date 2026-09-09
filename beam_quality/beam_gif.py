"""Animated GIF of NuMI beam position at the target vs time.

Reads the tidy CSV written by get_data.cpp and plots the INDIVIDUAL
batch-by-batch BPM readings (array indices 1-6, zeros = bad readings
excluded) — no averaging at all, so the per-batch "lobe" structure is
visible. Each batch's 121- and TGT-station readings are linearly
extrapolated to the target z (NOvA's BpmProjection geometry), so the
plotted coordinate is the one the goodbeam position cut is defined in.
Each frame shows the 60 batch points of 10 consecutive spills. Axes are
never clipped: they span the union of the data and the NOvA "goodbeam"
position box cut (IFDBSpillInfo.fcl: posx and posy within +-2 mm at
target), which is drawn on the plot. Rendered at 3000x2000.

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
from PIL import Image
from datetime import datetime, timedelta, timezone
import os, sys

BASE = os.environ.get("BEAMQ_DIR", os.getcwd())
# usage: python3 beam_gif.py [get_data_<t0>_to_<t1>]
STEM = sys.argv[1] if len(sys.argv) > 1 else \
    "get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500"
CSV = f"{BASE}/{STEM}.csv"
OUT = f"{BASE}/beam_position_spills_{STEM.replace('get_data_', '')}.gif"
SPILLS_PER_FRAME = 10
TRAIL_SPILLS = 50          # how many previous spills stay visible, fading
CUT = 2.0                  # NOvA goodbeam position cut: |x|,|y| < 2 mm
CDT = timezone(timedelta(hours=-5))

# palette (dataviz reference, light mode, committed single look)
BLUE = "#2a78d6"      # series slot 1
RED = "#e34948"       # status "serious" — reserved for the cut boundary
INK = "#1a1a19"       # text primary
INK2 = "#6b6a60"      # text secondary
CTX = "#c9c8c0"       # context marks
GRID = "#e8e7e0"
SURF = "#ffffff"

print("loading CSV...", flush=True)
df = pd.read_csv(CSV)

def batches(device, col):
    d = df[(df.device == device) & (df["index"] >= 1) & (df.value != 0.0)]
    return d.rename(columns={"time_s": "t", "value": col, "index": "idx"})[["t", "idx", col]]

# Surveyed BPM station z-positions [feet] and the NOvA linear
# extrapolation of each batch to the target (z = 0): the same constants
# and formula as computeBpmPosition() in get_data.cpp (ported from
# NOvA's IFDBSpillInfo extrapolate_position/BpmProjection), so the
# plotted coordinates match the ones the goodbeam position cut acts on.
Z_HP121, Z_VP121 = -68.04458, -66.99283
Z_HPTGT, Z_VPTGT = -31.25508, -30.16533
FX = (0.0 - Z_HP121) / (Z_HPTGT - Z_HP121)
FY = (0.0 - Z_VP121) / (Z_VPTGT - Z_VP121)

HT = batches("E:HPTGT[]", "xt"); HU = batches("E:HP121[]", "xu")
VT = batches("E:VPTGT[]", "yt"); VU = batches("E:VP121[]", "yu")
# pair the four devices per batch index; same-spill rows differ by ~ms
pieces = []
for i in range(1, 7):
    p = HT[HT.idx == i][["t", "xt"]].sort_values("t")
    for d, col in ((HU, "xu"), (VT, "yt"), (VU, "yu")):
        p = pd.merge_asof(p, d[d.idx == i][["t", col]].sort_values("t"),
                          on="t", tolerance=0.5, direction="nearest")
    p = p.dropna()
    p["x"] = p.xu + (p.xt - p.xu) * FX   # extrapolated to target z
    p["y"] = p.yu + (p.yt - p.yu) * FY
    pieces.append(p)
P = pd.concat(pieces).sort_values("t", kind="stable").reset_index(drop=True)
spill_times = np.sort(P["t"].unique())
P["spill"] = np.searchsorted(spill_times, P["t"])
print(f"spills: {len(spill_times)}, batch points: {len(P)}", flush=True)

n_frames = len(spill_times) // SPILLS_PER_FRAME
t0 = spill_times[0]
x = P["x"].values
y = P["y"].values
tmin = (P["t"].values - t0) / 60.0
spill = P["spill"].values

# never clip: axes cover the data AND the goodbeam box
def lims(c, pad=0.06):
    lo, hi = min(np.min(c), -CUT), max(np.max(c), CUT)
    p = (hi - lo) * pad
    return lo - p, hi + p
xlim, ylim = lims(x), lims(y)
tlim = (tmin[0] - 1, tmin[-1] + 1)

fig = plt.figure(figsize=(15, 10), dpi=200)
fig.patch.set_facecolor(SURF)
gs = GridSpec(2, 2, width_ratios=[1.35, 1], hspace=0.32, wspace=0.22,
              left=0.06, right=0.97, top=0.86, bottom=0.08)
axM = fig.add_subplot(gs[:, 0])
axX = fig.add_subplot(gs[0, 1])
axY = fig.add_subplot(gs[1, 1])

for ax in (axM, axX, axY):
    ax.set_facecolor(SURF)
    ax.grid(color=GRID, linewidth=0.8)
    ax.tick_params(colors=INK2, labelsize=11)
    for sp in ax.spines.values():
        sp.set_color(GRID)

win0 = datetime.fromtimestamp(spill_times[0], CDT)
win1 = datetime.fromtimestamp(spill_times[-1], CDT)
fig.suptitle("NuMI beam position at target  (E:HPTGT / E:VPTGT, IFBeam $A9 spills)",
             fontsize=19, color=INK, x=0.06, ha="left", y=0.965)
fig.text(0.06, 0.915,
         f"{win0:%Y-%m-%d  %H:%M} → {win1:%H:%M} CDT  ·  batch positions extrapolated to target z, "
         "10 spills per frame  ·  red box: NOvA goodbeam cut (|x|,|y| < 2 mm)",
         fontsize=12, color=INK2)

# --- main panel: beam spot ---
axM.set_xlabel("Horizontal position at target [mm]", fontsize=13, color=INK)
axM.set_ylabel("Vertical position at target [mm]", fontsize=13, color=INK)
axM.set_xlim(*xlim); axM.set_ylim(*ylim)
axM.scatter(x, y, s=4, color=CTX, linewidths=0, zorder=1)   # full-hour context
axM.add_patch(Rectangle((-CUT, -CUT), 2 * CUT, 2 * CUT, fill=False,
                        edgecolor=RED, linewidth=2.2, linestyle=(0, (6, 3)), zorder=2))
axM.text(-CUT + 0.06, -CUT + 0.06, "NOvA goodbeam position cut (±2 mm)",
         color=RED, fontsize=11, va="bottom")
trail_sc = axM.scatter([], [], s=16, linewidths=0, zorder=3)
cur_sc = axM.scatter([], [], s=42, color=BLUE, edgecolors=SURF,
                     linewidths=0.8, zorder=5)
stamp = axM.text(0.03, 0.97, "", transform=axM.transAxes, fontsize=13,
                 color=INK, va="top", family="monospace")

# --- time-series panels: every batch reading, no averaging ---
for ax, c, lab in ((axX, x, "Horizontal"), (axY, y, "Vertical")):
    ax.set_xlim(*tlim)
    ax.set_ylim(*lims(c))
    ax.scatter(tmin, c, s=1.5, color=CTX, linewidths=0)
    ax.axhline(+CUT, color=RED, lw=1.4, linestyle=(0, (6, 3)))
    ax.axhline(-CUT, color=RED, lw=1.4, linestyle=(0, (6, 3)))
    ax.set_title(f"{lab} position vs time (per batch)", fontsize=13, color=INK, loc="left")
    ax.set_ylabel("[mm]", fontsize=12, color=INK)
axY.set_xlabel(f"Minutes after {win0:%H:%M} CDT", fontsize=12, color=INK)
progX = axX.scatter([], [], s=1.5, color=BLUE, linewidths=0)
progY = axY.scatter([], [], s=1.5, color=BLUE, linewidths=0)
curX = axX.axvline(0, color=INK2, lw=1.2, alpha=0.7)
curY = axY.axvline(0, color=INK2, lw=1.2, alpha=0.7)

blue_rgb = matplotlib.colors.to_rgb(BLUE)

frames = []
for i in range(n_frames):
    a, b = i * SPILLS_PER_FRAME, (i + 1) * SPILLS_PER_FRAME
    curm = (spill >= a) & (spill < b)
    trailm = (spill >= max(0, a - TRAIL_SPILLS)) & (spill < a)
    # fading trail of the previous TRAIL_SPILLS spills' batch points
    trail_sc.set_offsets(np.c_[x[trailm], y[trailm]])
    age = a - spill[trailm]                       # 1 (newest) .. TRAIL_SPILLS
    alphas = 0.45 * (1.0 - age / (TRAIL_SPILLS + 1))
    trail_sc.set_facecolors([(*blue_rgb, al) for al in alphas])
    # the 10 current spills' individual batch points, full strength
    cur_sc.set_offsets(np.c_[x[curm], y[curm]])
    ta, tb = spill_times[a], spill_times[b - 1]
    ts0 = datetime.fromtimestamp(ta, CDT).strftime("%H:%M:%S")
    ts1 = datetime.fromtimestamp(tb, CDT).strftime("%H:%M:%S")
    stamp.set_text(f"{ts0}-{ts1} CDT   spills {a}-{b - 1}")
    donem = spill < b
    progX.set_offsets(np.c_[tmin[donem], x[donem]])
    progY.set_offsets(np.c_[tmin[donem], y[donem]])
    curX.set_xdata([(tb - t0) / 60.0]); curY.set_xdata([(tb - t0) / 60.0])

    fig.canvas.draw()
    img = Image.frombuffer("RGBA", fig.canvas.get_width_height(),
                           fig.canvas.buffer_rgba()).convert("RGB")
    frames.append(img.quantize(colors=128, method=Image.MEDIANCUT))
    if i % 25 == 0:
        print(f"frame {i}/{n_frames}", flush=True)

print("assembling GIF...", flush=True)
frames[0].save(OUT, save_all=True, append_images=frames[1:],
               duration=80, loop=0, optimize=True)
print(f"done: {OUT}  {os.path.getsize(OUT)/1e6:.1f} MB, "
      f"{len(frames)} frames, {frames[0].size[0]}x{frames[0].size[1]}", flush=True)
