#!/usr/bin/env python3
"""Extract per-gate MINERvA (Mx2) time slices with tracks into compact .npz.

One output per DST: <out>/<DST basename>.mx2slices.npz with
  trig_utime_ms, trig_type (=1 for every gate in a numib file), trig_gate,
  trig_n_slices, trig_n_tracks                    (one row per gate)
  flash_event, flash_slice, flash_t_rel, flash_t_min, flash_t_med, flash_nclus,
  flash_pe, flash_zmin, flash_zmax, flash_ntrk, flash_trk_hits_max,
  flash_trk_nodes_max, flash_theta_min            (one row per time slice >= 1)
  meta (json string)
The "flash_" prefix is kept so scripts/batch_rates.py reads 2x2 and Mx2 alike.
See bpscan/mx2_tracks.py for what each quantity means and why slices, not
tracks, are counted.

Inputs: local DST paths (NERSC /global/cfs/cdirs/dune/www/data/2x2/Mx2/dst/)
or https://portal.nersc.gov/project/dune/data/2x2/Mx2/dst/... URLs. Only the
listed branches are read. Existing outputs are skipped (resumable); a failure
writes <name>.failed.txt and the loop continues.

Files covering the July 2024 IFBeam windows used in the report:
  TS1_00050017_0005_numib_v09_2407090337_DST.root   (Jul 8 22:37 - Jul 9 01:02 CDT)
  TS1_00050017_0006_numib_v09_2407090615_DST.root   (Jul 9 01:15 - ~03:40 CDT)
  TS1_00050018_0017_numib_v09_2407120509_DST.root   (Jul 12 00:09 - ~02:35 CDT)
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import mx2_tracks as mx  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--max-entries", type=int, default=None, help="debug: read only the first N gates")
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for src in a.inputs:
        name = Path(str(src).split("?")[0]).name
        target = out / (name + ".mx2slices.npz")
        if target.exists():
            print("skip", name); continue
        t0 = time.time()
        try:
            g, s = mx.read_dst(src, max_entries=a.max_entries)
        except Exception as err:
            print("FAIL", name, type(err).__name__, err, flush=True)
            (out / (name + ".failed.txt")).write_text(f"{type(err).__name__}: {err}\n")
            continue
        meta = dict(source=str(src), n_gates=int(len(g)), n_slices=int(len(s)),
                    n_track_slices=int((s.ntrk > 0).sum()), trigger_types=sorted(map(int, set(g.trigger_type))),
                    seconds=round(time.time() - t0, 1))
        tmp = target.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, trig_utime_ms=g.utime_ms.values, trig_type=np.ones(len(g), np.int32),
                            trig_gate=g.gate.values, trig_n_slices=g.n_slices.values,
                            trig_n_tracks=g.n_tracks.values, trig_raw_type=g.trigger_type.values,
                            **{f"flash_{c}": s[c].values for c in s.columns}, meta=json.dumps(meta))
        tmp.rename(target)
        print("ok", name, json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
