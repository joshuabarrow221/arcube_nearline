# beam_position_scan — NuMI proton-beam position on target vs. 2x2 / Mx2 rates

Can the 2x2 ND-LAr demonstrator (and MINERvA+2x2, "Mx2") tell the NuMI
accelerator where the proton beam should sit on the target? This directory
holds the analysis code, the feasibility study, and the write-up
(`report/beam_position_scan.pdf`).

**Read the report first.** Its executive summary states what is established,
what is model-dependent, and what has not been run on detector data yet.

Builds on `../beam_quality/` (IFBeam fetcher + NOvA goodbeam cuts), started by
**Gianfranco Ingratta** with **Bruce Howard** (York U.) and extended by
J. L. Barrow (UMN).

## The idea in one paragraph

Each NuMI spill is six ~1.58 µs Booster batches. Batch 1 is extracted on the
rising edge of the kicker and lands ~0.3 mm (horizontal) away from batches 2-6
on the target, every spill. If the detector can tell which batch an
interaction came from (it can, with the light system or MINERvA timing, which
resolve ~10 ns against a 1.58 µs batch), then batch 1 and batches 2-6 are two
simultaneous "beam positions" sharing every detector and beam condition. A
per-POT rate difference between them measures the local slope of rate vs.
beam position, i.e. which side of the optimum the beam is on. The same
machinery generalises to a dedicated scan: count per-POT rates at each
setpoint and fit the known interception shape for the target-fin centre.

## Status (keep this table current)

| Step | State | Where |
|---|---|---|
| Spill/batch timing model, batch assignment with guard bands | done, unit-tested | `bpscan/numi.py` |
| IFBeam per-batch positions at target z, per-batch POT share | done, tested on synthetic data; run on July 9/12 2024 exports | `bpscan/ifbeam.py`, `scripts/analyze_ifbeam_batches.py` |
| Interception model, Poisson sensitivity, scan toy MC | done | `bpscan/target_model.py`, `bpscan/stats.py`, `scripts/sensitivity_study.py` |
| 2x2 FLOW light reader -> per-trigger flashes (.npz), local or HTTP byte ranges | done; run on 19 v11 files (run 50017 2024-07-09 23:54-02:05 CDT, run 50018 2024-07-12 00:00-01:00 CDT), outputs in `results/flashes/` | `bpscan/flow_light.py`, `scripts/extract_lrs_flashes.py` |
| Spill matching, timing-drift veto, duplicate-trigger veto, batch-edge fit, per-batch per-POT rates, masking test, follower-bias estimate | done on the 19 files (7,868 spills, 11,937 flashes): batch1/rest = 0.992 +- 0.026 (all), 1.000 +- 0.030 (>=2 TPCs). Null, as predicted. Delayed followers bias the ratio by -2.1% (all) / -0.3% (>=2 TPCs) | `bpscan/detector_rates.py`, `scripts/batch_rates.py`, `results/batch_rates*/` |
| Run 1 at scale (all v11 files, ~2.2e19 POT + matching IFBeam per-batch data) | **not run** - needs IFBeam per-batch exports for the whole run and the extractor at NERSC; expected 3-sigma reach ~1% | see "Picking this up" |
| Data defect | two v11 files (`packet-0050017-2024_07_09_01_14_40`, `packet-0050018-2024_07_12_00_00_19`) contain duplicated light beam triggers (687 and 456 triggers share a spill). Dropped; report to reflow maintainers | `batch_rates.json: duplicate_triggers_by_file` |
| Mx2 (MINERvA) reader | **not started**. DSTs exist: `portal.nersc.gov/project/dune/data/2x2/Mx2/dst/` (502 files, ~374 MB each; runs 50000-50018). Tree `minerva`, one entry per gate, per-hit `hit_time`, `hit_module`, `hit_plane`, `hit_pe` - enough to assign hits to batches | — |

## Layout

```
beam_position_scan/
  README.md                 this file
  bpscan/                   importable library (no detector-format code inside)
    numi.py                 RF bucket / batch structure, BPM z geometry, batch_index()
    ifbeam.py               get_data.cpp CSV -> per-(spill,batch) table, batch1-minus-rest
    target_model.py         Gaussian beam on finite fin: interception, curvature, R(x)
    stats.py                ratio tests, required counts, scan toy MC (parabola/template)
    detector_rates.py       spill matching, batch-edge calibration, per-batch per-POT rates
  scripts/
    analyze_ifbeam_batches.py   step 1: is batch 1 really elsewhere on target? lever arms
    sensitivity_study.py        step 3: Run 1 reach, ~100 kW scan precision, time budgets
  tests/test_bpscan.py      pytest suite (python3 -m pytest -q tests/)
  results/                  JSON + figures written by the scripts (committed for the report)
  report/                   LaTeX source + built PDF
```

## Quick start

```bash
pip install numpy scipy pandas matplotlib h5py pytest
python3 -m pytest -q tests/

# 1. Beam side (needs get_data.cpp outputs; see ../beam_quality/README.md)
python3 scripts/analyze_ifbeam_batches.py \
    /path/get_data_2024-07-09T000000-0500_to_2024-07-09T020500-0500 \
    /path/get_data_2024-07-12T000100-0500_to_2024-07-12T010100-0500 \
    --out results/ifbeam

# 2. Feasibility (rates are assumptions until measured; pass measured ones)
python3 scripts/sensitivity_study.py --sigma 1.43 --width 7.4 --rate-at-100kw 50 150 500 2000

# 3. Report
cd report && pdflatex beam_position_scan.tex && pdflatex beam_position_scan.tex
```

## Rules this code enforces (do not relax them silently)

1. **Rates are per POT, per batch.** Batch intensities differ spill to spill;
   POT per batch = spill POT × that batch's BPM-intensity share.
2. **Unmatched or mistimed spills are dropped, not corrected.** A detector
   spill must match an IFBeam $A9 record within the tolerance, and its
   reference (RWM) time must be present and consistent; otherwise the spill is
   removed and counted in the audit. The goodbeam selection is applied on the
   beam side.
3. **Guard bands around batch edges.** Objects within `guard_ns` of any batch
   boundary are not assigned to a batch; the POT is scaled by the surviving
   fraction of the batch duration. Use a wider guard around batch 1 (kicker
   rise) as requested.
4. **Batch edges are calibrated from data** (`fit_batch_edges`), never assumed
   from cable-length arithmetic. Look at the time histogram: six plateaus
   should be visible.
5. **No simulated numbers stand in for data.** Feasibility numbers are labelled
   as model/assumption in every output JSON.

## Picking this up (for the next person or agent)

* Data access from a Claude cloud container is blocked (portal.nersc.gov and
  dbdata3vm.fnal.gov return proxy 403). Run fetches on NERSC or a laptop.
* The FLOW streaming pattern (HTTP byte ranges, h5py over fsspec, strong
  ETag) lives on branch `dev/purity-overlay`, `actions/lifetime/flow_input.py`.
* Next concrete task: run the extractor over all Run 1 v11 beam FLOW files
  (at NERSC: `/global/cfs/cdirs/dune/www/data/2x2/reflows/v11/flow/beam/*/*/*.FLOW.hdf5`),
  fetch per-batch IFBeam data for the same span with `../beam_quality/get_data.cpp`
  in 1-hour windows, then
  `python3 scripts/batch_rates.py --flashes flashes/*.npz --ifbeam <stems> --out results/batch_rates`.
  Look at `masking_test` in the output before quoting the batch-1 ratio.
* Then: an Mx2 DST reader emitting (gate time, hit_time relative to the spill)
  for the front veto planes, fed to the same rate code.
* Note: in `../beam_quality/README.md` / `get_data.cpp` the horn current reads
  about -200 kA and is described as FHC, but Run 1 was antineutrino (RHC) running.
  The stripline sum does not carry polarity; do not use the horn cut to tell FHC from RHC.
* Open physics questions are listed at the end of the report.
