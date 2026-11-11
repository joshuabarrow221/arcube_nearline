"""Turn beam_quality/get_data.cpp outputs into per-spill, per-batch tables.

Inputs (written by get_data.cpp for a window <stem>):
    <stem>.csv          long format: device,time_s,time_iso,index,value
                        array devices (E:HPTGT[] ...) have index 0..6,
                        element 0 = BPM auto-tune average (discarded),
                        elements 1..6 = Booster batches 1..6.
    <stem>_quality.csv  one row per spill on the reference timeline with
                        spillpot, hornI, posx, posy, widthx, widthy, pass_*,
                        goodbeam (NOvA goodbeam cuts 1-5).

Output of load_batches(): a pandas DataFrame with one row per (spill, batch):
    spill        integer spill id (row number of the quality CSV)
    time_s       quality-timeline unix time of the spill [s]
    batch        1..6
    x_tgt, y_tgt per-batch position extrapolated to the target z [mm]
    x121, y121, xtgtbpm, ytgtbpm   raw BPM readings [mm]
    ix, iy       per-batch BPM intensities (arbitrary units; used as weights)
    frac         per-batch share of the spill intensity (from ix)
    pot_batch    spillpot * frac  [protons]
    goodbeam ... copied from the quality CSV

Matching: every device is matched to the quality timeline with
pandas.merge_asof within +-0.5 s (same tolerance as get_data.cpp). Readings of
exactly 0 for a batch are treated as missing (the acquisition's "bad batch"
convention, as in computeBpmPosition()).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from . import numi

ARRAY_DEVICES = {
    "E:HP121[]": "x121", "E:VP121[]": "y121",
    "E:HPTGT[]": "xtgtbpm", "E:VPTGT[]": "ytgtbpm",
    "E:HITGT[]": "ix", "E:VITGT[]": "iy",
}
MATCH_TOL_S = 0.5


def _read_long(csv_path, devices):
    """Read only the requested devices from the (large) long-format CSV."""
    chunks = []
    for chunk in pd.read_csv(csv_path, chunksize=2_000_000,
                             usecols=["device", "time_s", "index", "value"]):
        chunks.append(chunk[chunk["device"].isin(devices)])
    return pd.concat(chunks, ignore_index=True)


def load_quality(stem):
    q = pd.read_csv(f"{stem}_quality.csv")
    q.insert(0, "spill", np.arange(len(q)))
    return q


def load_batches(stem):
    """Return the per-(spill, batch) DataFrame described in the module docstring."""
    stem = str(stem)
    q = load_quality(stem)
    raw = _read_long(f"{stem}.csv", list(ARRAY_DEVICES))
    raw = raw[raw["index"].between(1, numi.N_BATCHES)]
    ref = q[["spill", "time_s"]].sort_values("time_s")

    per_dev = []
    for dev, col in ARRAY_DEVICES.items():
        d = raw[raw["device"] == dev].pivot_table(index="time_s", columns="index",
                                                  values="value", aggfunc="first")
        d = d.reset_index().sort_values("time_s").rename(columns={"time_s": "t_dev"})
        m = pd.merge_asof(ref, d, left_on="time_s", right_on="t_dev",
                          direction="nearest", tolerance=MATCH_TOL_S)
        long = m.melt(id_vars=["spill", "time_s", "t_dev"],
                      value_vars=list(range(1, numi.N_BATCHES + 1)),
                      var_name="batch", value_name=col)
        per_dev.append(long[["spill", "batch", col]].set_index(["spill", "batch"]))
    b = pd.concat(per_dev, axis=1).reset_index()
    b["batch"] = b["batch"].astype(int)
    for col in ARRAY_DEVICES.values():
        b.loc[b[col] == 0.0, col] = np.nan              # bad-batch convention

    b["x_tgt"] = numi.extrapolate(b["x121"], numi.Z_HP121, b["xtgtbpm"], numi.Z_HPTGT)
    b["y_tgt"] = numi.extrapolate(b["y121"], numi.Z_VP121, b["ytgtbpm"], numi.Z_VPTGT)
    # Per-batch share of the spill: horizontal BPM intensity (clipped at 0),
    # normalised over the batches that have a valid reading in that spill.
    ixc = b["ix"].clip(lower=0)
    b["frac"] = ixc / ixc.groupby(b["spill"]).transform("sum")
    b = b.merge(q, on="spill", how="left", suffixes=("", "_q"))
    b["pot_batch"] = b["spillpot"] * b["frac"]
    return b.sort_values(["spill", "batch"]).reset_index(drop=True)


def batch_summary(b, select=None):
    """Mean, RMS and standard error of per-batch positions and intensity shares."""
    d = b if select is None else b[select]
    g = d.groupby("batch")
    out = pd.DataFrame({
        "n_spills": g["x_tgt"].count(),
        "x_mean": g["x_tgt"].mean(), "x_rms": g["x_tgt"].std(),
        "y_mean": g["y_tgt"].mean(), "y_rms": g["y_tgt"].std(),
        "frac_mean": g["frac"].mean(), "frac_rms": g["frac"].std(),
    })
    out["x_sem"] = out["x_rms"] / np.sqrt(out["n_spills"])
    out["y_sem"] = out["y_rms"] / np.sqrt(out["n_spills"])
    return out


def batch1_minus_rest(b, select=None):
    """Per-spill (batch 1) - (intensity-weighted mean of batches 2-6) offsets.

    Doing the subtraction spill by spill removes common-mode spill-to-spill
    steering jitter, so the result isolates the batch-1 trajectory difference.
    Returns a DataFrame indexed by spill with dx, dy [mm].
    """
    d = b if select is None else b[select]
    one = d[d["batch"] == 1].set_index("spill")[["x_tgt", "y_tgt"]]
    rest = d[d["batch"] > 1].copy()
    w = rest["frac"].fillna(0)
    rest["wx"], rest["wy"] = rest["x_tgt"] * w, rest["y_tgt"] * w
    rest["w"] = w.where(rest["x_tgt"].notna(), 0)
    s = rest.groupby("spill")[["wx", "wy", "w"]].sum()
    mean_rest = pd.DataFrame({"x_tgt": s["wx"] / s["w"], "y_tgt": s["wy"] / s["w"]})
    out = (one - mean_rest).rename(columns={"x_tgt": "dx", "y_tgt": "dy"})
    return out.dropna()


def find_stems(directory):
    """All get_data stems (files with a matching _quality.csv) in a directory."""
    return sorted(str(p)[: -len("_quality.csv")] for p in Path(directory).glob("*_quality.csv"))
