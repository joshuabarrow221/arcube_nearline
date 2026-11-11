"""Read the light-system timing needed for per-batch rates from 2x2 FLOW files.

Why the light system: a NuMI spill is six ~1.58 us batches. Charge-readout
(CRS) timestamps of an interaction are smeared by the unknown drift time
(0-~200 us) unless a t0 is supplied, so the CRS alone cannot say which batch an
interaction came from. The light readout (LRS) samples at 62.5 MHz (16 ns)
in a ~16 us window opened by the beam trigger, so every light pulse has a time
relative to the trigger good to ~10-20 ns -- ~100x finer than a batch.

What this module reads (2x2 Run 1 FLOW, reflow v11, flow v1.5.0 layout):
  light/events/data        per trigger: trig_type (1 = beam, 0 = other; see
                           note), utime_ms[8], tai_ns[8] (per ADC)
  light/sum_hits/data      per summed-channel pulse: tpc, det, busy_ns
                           (= peak time relative to the trigger [ns], i.e.
                           sample_idx*16), ns_spline (sub-sample correction,
                           ns), sum, max
  light/events/ref/light/sum_hits/ref   event <-> hit association (preferred)

trig_type == 1 matched the IFBeam $A9 spill count in the file checked
(454 light beam triggers vs 455 IFBeam spills in 10 min) but has not been
checked against LRS documentation. Callers should still match each beam
trigger to an IFBeam spill (detector_rates.match_spills) and drop unmatched
ones.

NOTE on the `ns` field of sum_hits: in v11 it is the *event* timestamp stored
as float64 unix ns (precision ~256 ns at 1.7e18), not the hit time. Never use it
for in-spill timing; use busy_ns + ns_spline.

Remote reading: HTTP byte ranges through fsspec. FLOW datasets use tiny lzf
chunks (light/events: 16 rows/chunk), so the default 256 kB block size of the
lifetime reader over-fetches ~50x. Use block_size ~16-32 kB here, or better,
run on NERSC with local paths (no network at all).
"""
from contextlib import contextmanager

import numpy as np
import h5py

EVENTS = "light/events/data"
TNS = "light/t_ns/data"
HITS = "light/sum_hits/data"
HIT_REF = "light/events/ref/light/sum_hits/ref"
HIT_FIELDS = ("tpc", "det", "busy_ns", "ns_spline", "sum", "max")


@contextmanager
def open_any(path, block_size=32768):
    """Open a local or HTTP(S) FLOW file read-only (HTTP: byte ranges)."""
    if not str(path).startswith(("http://", "https://")):
        with h5py.File(path, "r") as f:
            yield f
        return
    import fsspec
    fs = fsspec.filesystem("http", skip_instance_cache=True, client_kwargs={"trust_env": True})
    with fs.open(path, block_size=block_size, cache_type="blockcache",
                 headers={"Accept-Encoding": "identity"}) as stream:
        with h5py.File(stream, "r") as f:
            yield f


def read_light(f):
    """Return dict of numpy arrays: events (trig_type, utime_ms, tai_ns, t_ns) and hits.

    hits['event'] is the light-event row each hit belongs to (from the ref
    dataset when present, otherwise by matching the hit's `ns` field to the
    event's t_ns, which is exact in v11 because `ns` IS the event time).
    """
    ev = f[EVENTS].fields(["trig_type", "utime_ms", "tai_ns"])[:]
    t_ns = f[TNS]["t_ns"][:, 0] if TNS in f else ev["tai_ns"][:, 0].astype(float)
    hits = f[HITS].fields(list(HIT_FIELDS) + ["ns"])[:]
    if HIT_REF in f:
        ref = f[HIT_REF][:]
        hit_event = np.full(len(hits), -1, np.int64)
        hit_event[ref[:, 1]] = ref[:, 0]
        how = "ref"
    else:
        order = np.argsort(t_ns)
        j = np.clip(np.searchsorted(t_ns[order], hits["ns"]), 1, len(order) - 1)
        left, right = order[j - 1], order[j]
        hit_event = np.where(np.abs(hits["ns"] - t_ns[left]) <= np.abs(hits["ns"] - t_ns[right]), left, right)
        ok = np.abs(hits["ns"] - t_ns[hit_event]) < 1e3
        hit_event = np.where(ok, hit_event, -1)
        how = "time-match"
    return dict(trig_type=ev["trig_type"], utime_ms=ev["utime_ms"][:, 0].astype(np.int64),
                tai_ns=ev["tai_ns"][:, 0].astype(np.int64), t_ns=t_ns,
                hit_event=hit_event, hit_t_rel=hits["busy_ns"] + hits["ns_spline"],
                hit_tpc=hits["tpc"], hit_det=hits["det"], hit_sum=hits["sum"], hit_max=hits["max"],
                association=how)


def cluster_flashes(event, t_rel, tpc, amp, gap_ns=80.0):
    """Group hits of the same trigger into flashes (one per physical interaction).

    Hits are merged across all TPCs when consecutive peak times differ by less
    than gap_ns (default 80 ns = 5 samples). A through-going rock muon lights
    several TPCs simultaneously and becomes ONE flash. Returns a structured
    array: event, t_rel (earliest hit), nhits, ntpc, tpc_mask, amp_sum.
    """
    order = np.lexsort((t_rel, event))
    e, t, k, a = event[order], t_rel[order], tpc[order].astype(np.int64), amp[order]
    new = np.ones(len(e), bool)
    new[1:] = (e[1:] != e[:-1]) | (np.diff(t) > gap_ns)
    gid = np.cumsum(new) - 1
    n = gid[-1] + 1 if len(gid) else 0
    out = np.zeros(n, dtype=[("event", "i8"), ("t_rel", "f8"), ("nhits", "i4"), ("ntpc", "i4"),
                             ("tpc_mask", "i4"), ("amp_sum", "f8")])
    if n == 0:
        return out
    out["event"] = e[new]
    out["t_rel"] = t[new]
    out["nhits"] = np.bincount(gid)
    mask = np.zeros(n, np.int64)
    np.bitwise_or.at(mask, gid, 1 << k)
    out["tpc_mask"] = mask
    out["ntpc"] = [bin(m).count("1") for m in mask]
    out["amp_sum"] = np.bincount(gid, weights=a)
    return out
