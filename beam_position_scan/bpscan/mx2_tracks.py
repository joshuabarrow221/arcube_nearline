"""MINERvA (Mx2) DST reader: per-gate time slices with reconstructed tracks.

Why this exists
  The 2x2 light flashes give ~1.5 objects per spill at ~800 kW. MINERvA, which
  surrounds the 2x2 in the Mx2 configuration, sees the same NuMI flux over a
  much larger area and mass (plus rock muons), so it should give many more
  timed, track-like objects per spill. More objects per POT means a smaller
  Poisson error on the batch-1 / batches-2-6 ratio and on each scan point.

What a DST holds (checked on TS1_00050017_0005_numib_v09_2407090337_DST.root)
  TTree "minerva", one entry per beam gate (1.2 s apart), 219 branches.
  Used here:
    ev_gps_time_sec, ev_gps_time_usec     gate GPS time (UTC unix, 1 us)
    ev_trigger_type, ev_gate              trigger type (16 in numib files), gate number
    n_slices                              number of time slices in the gate
    clus_id_time, clus_id_time_slice,     inner-detector clusters: time [ns, gate clock],
    clus_id_pe, clus_id_z                 time-slice index, light, z position [mm]
    n_tracks, trk_time_slice, trk_hits,   reconstructed tracks and the slice they belong
    trk_nodes, trk_theta                  to; hits and nodes on the track; polar angle
  Not used (heavy, or not needed for timing): trk_node_* [][300], blobs, MC.

The unit of counting is a TIME SLICE that contains at least one reconstructed
track ("track-like slice"), not a track. One interaction (or one rock muon
with its delta rays) can give several tracks; counting slices keeps the
objects close to independent Poisson events, which the error bars assume.
Slice index 0 is MINERvA's "unsliced" bucket and is skipped [likely - check
against the MINERvA slicer docs if a reviewer asks].

Slice time = light-weighted mean of its cluster times. t_min (earliest
cluster) is stored as well; the batch-edge fit tells which aligns better.

Output schema deliberately mirrors the 2x2 flash .npz (trig_* per gate,
flash_* per object) so scripts/batch_rates.py can run on either; select Mx2
objects there with --select (a pandas query on the flash_* columns).
"""
import numpy as np
import pandas as pd

BRANCHES = [
    "ev_gps_time_sec", "ev_gps_time_usec", "ev_trigger_type", "ev_gate", "n_slices",
    "clus_id_time", "clus_id_time_slice", "clus_id_pe", "clus_id_z",
    "n_tracks", "trk_time_slice", "trk_hits", "trk_nodes", "trk_theta",
]


def _flat(arr, gate0):
    """Jagged awkward array -> (flat numpy values, gate index per value)."""
    import awkward as ak
    counts = ak.to_numpy(ak.num(arr))
    gate = np.repeat(np.arange(len(counts)) + gate0, counts)
    return ak.to_numpy(ak.flatten(arr)), gate


def slices_from_chunk(ch, gate0=0):
    """Build one row per (gate, time slice >= 1) from a chunk of DST entries.

    ch: dict-like of awkward arrays for BRANCHES (from uproot iterate).
    gate0: row offset of this chunk within the file.
    Returns (gates DataFrame, slices DataFrame).
    """
    import awkward as ak
    gates = pd.DataFrame({
        "row": np.arange(len(ch["ev_gate"])) + gate0,
        "utime_ms": ak.to_numpy(ch["ev_gps_time_sec"]).astype(np.float64) * 1e3
                    + ak.to_numpy(ch["ev_gps_time_usec"]) / 1e3,
        "trigger_type": ak.to_numpy(ch["ev_trigger_type"]),
        "gate": ak.to_numpy(ch["ev_gate"]),
        "n_slices": ak.to_numpy(ch["n_slices"]),
        "n_tracks": ak.to_numpy(ch["n_tracks"]),
    })
    t, g = _flat(ch["clus_id_time"], gate0)
    c = pd.DataFrame({"event": g, "slice": _flat(ch["clus_id_time_slice"], gate0)[0], "t": t,
                      "pe": _flat(ch["clus_id_pe"], gate0)[0], "z": _flat(ch["clus_id_z"], gate0)[0]})
    c = c[c.slice >= 1]
    c["tpe"] = c.t * c.pe
    s = c.groupby(["event", "slice"]).agg(nclus=("t", "size"), pe=("pe", "sum"), tpe=("tpe", "sum"),
                                          t_min=("t", "min"), t_med=("t", "median"),
                                          zmin=("z", "min"), zmax=("z", "max")).reset_index()
    s["t_rel"] = np.where(s.pe > 0, s.tpe / s.pe.where(s.pe > 0, 1), s.t_med)
    s = s.drop(columns="tpe")
    ts, tg = _flat(ch["trk_time_slice"], gate0)
    tr = pd.DataFrame({"event": tg, "slice": ts, "hits": _flat(ch["trk_hits"], gate0)[0],
                       "nodes": _flat(ch["trk_nodes"], gate0)[0], "theta": _flat(ch["trk_theta"], gate0)[0]})
    tr = tr.groupby(["event", "slice"]).agg(ntrk=("hits", "size"), trk_hits_max=("hits", "max"),
                                            trk_nodes_max=("nodes", "max"),
                                            theta_min=("theta", "min")).reset_index()
    s = s.merge(tr, on=["event", "slice"], how="left")
    for col in ("ntrk", "trk_hits_max", "trk_nodes_max"):
        s[col] = s[col].fillna(0).astype(np.int32)
    s["theta_min"] = s["theta_min"].fillna(np.nan)
    return gates, s


def read_dst(path_or_url, step_size="40 MB", max_entries=None):
    """Read a DST (local path or https URL) and return (gates, slices) DataFrames."""
    import uproot
    gl, sl = [], []
    with uproot.open(path_or_url) as f:
        tree = f["minerva"]
        n = tree.num_entries if max_entries is None else min(tree.num_entries, max_entries)
        start = 0
        for ch in tree.iterate(BRANCHES, step_size=step_size, entry_stop=n, library="ak"):
            g, s = slices_from_chunk(ch, start)
            gl.append(g); sl.append(s)
            start += len(g)
    return pd.concat(gl, ignore_index=True), pd.concat(sl, ignore_index=True)
