#!/usr/bin/env python3
"""Extract per-beam-trigger light flashes from 2x2 FLOW files into compact .npz.

One output per input file: <out>/<input basename>.flashes.npz with
  trig_utime_ms, trig_tai_ns, trig_t_ns, trig_type      (all light triggers)
  hit_*                                                 (beam-trigger hits only)
  flash_*                                               (clustered, beam triggers only)
  meta (json string): source, size, association method, code version

Inputs can be local paths (NERSC: /global/cfs/cdirs/dune/www/data/2x2/...)
or https://portal.nersc.gov/project/dune/data/2x2/... URLs (byte-range
streaming; nothing else in the file is downloaded).

Examples
  # NERSC (fast, local filesystem):
  python3 extract_lrs_flashes.py --out flashes/ \
     /global/cfs/cdirs/dune/www/data/2x2/reflows/v11/flow/beam/july8_2024/nominal_hv/packet-0050017-2024_07_09_0*.FLOW.hdf5
  # laptop over HTTP:
  python3 extract_lrs_flashes.py --out flashes/ --block-size 32768 https://portal.nersc.gov/.../file.FLOW.hdf5
Existing outputs are skipped (resumable). Use -I with python when the inputs
directory is untrusted.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bpscan import flow_light as fl  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--block-size", type=int, default=32768)
    p.add_argument("--gap-ns", type=float, default=80.0)
    a = p.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for src in a.inputs:
        name = Path(str(src).split("?")[0]).name
        target = out / (name + ".flashes.npz")
        if target.exists():
            print("skip", name); continue
        t0 = time.time()
        try:
            with fl.open_any(src, a.block_size) as f:
                d = fl.read_light(f)
        except Exception as err:            # keep going; record the failure
            print("FAIL", name, type(err).__name__, err, flush=True)
            (out / (name + ".failed.txt")).write_text(f"{type(err).__name__}: {err}\n")
            continue
        beam_ev = np.where(d["trig_type"] == 1)[0]
        m = np.isin(d["hit_event"], beam_ev)
        fls = fl.cluster_flashes(d["hit_event"][m], d["hit_t_rel"][m], d["hit_tpc"][m], d["hit_sum"][m], a.gap_ns)
        meta = dict(source=str(src), association=d["association"], n_triggers=int(len(d["trig_type"])),
                    n_beam_triggers=int(len(beam_ev)), n_beam_hits=int(m.sum()), n_flashes=int(len(fls)),
                    gap_ns=a.gap_ns, seconds=round(time.time() - t0, 1))
        tmp = target.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, trig_utime_ms=d["utime_ms"], trig_tai_ns=d["tai_ns"], trig_t_ns=d["t_ns"],
                            trig_type=d["trig_type"], hit_event=d["hit_event"][m], hit_t_rel=d["hit_t_rel"][m],
                            hit_tpc=d["hit_tpc"][m], hit_det=d["hit_det"][m], hit_sum=d["hit_sum"][m],
                            hit_max=d["hit_max"][m], **{f"flash_{k}": fls[k] for k in fls.dtype.names},
                            meta=json.dumps(meta))
        tmp.rename(target)
        print("ok", name, json.dumps(meta), flush=True)


if __name__ == "__main__":
    main()
