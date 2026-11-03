# Fit selected tracks over a recent interval

`actions/lifetime/recent_tracks.py` fits the pooled segments from a trailing
**six elapsed hours** by default. It does not average per-file lifetime fits or
require individual files to fit successfully first. This is a rolling query;
the configured monitor's existing calendar six-hour bins remain separate.

For a NERSC validation run, use Python 3.10+ and install this repository's
analysis requirements into the chosen environment:

```bash
python -m pip install -r actions/lifetime/requirements.txt
```

The default `quality-stored` mode needs only this repository and those
dependencies; it does **not** need a separate `ndlar_flow` checkout, H5Flow,
database access, or any local OneDrive artifacts. The optional `reselect` mode
also requires an installed compatible `ndlar_flow` and its dependencies.

From the repository root, in the lifetime Python environment:

```bash
python actions/lifetime/recent_tracks.py /data/run/FLOW \
  --output /data/lifetime_studies/recent
```

Inputs can be files, recursively searched directories, or quoted glob patterns.
Keep each invocation within an acquisition cohort with stable detector settings.
The default `--selection quality-stored` checks the original hit references for
each stored segment using the hit-quality filter, preserving accepted segment
charge and length. It requires full FLOW files with the original associations.
Missing references are errors; the command does not silently substitute an
unvetted selection. For FLOW without stored track selections, explicitly use
`--selection reselect` to run the local `ndlar_flow` selector and quality filter;
this needs that package and its dependencies in the Python environment.

An optional `--quality-config /path/to/quality.json` supplies the quality-filter
parameters. Defaults and their limitations are described in the
[track-quality guide](lifetime-track-quality.md) and source. Thresholds need validation against detector data;
successful execution is not proof that every accepted segment is noise-free.

For a reproducible historical interval, specify its exclusive end with an
explicit UTC offset:

```bash
python actions/lifetime/recent_tracks.py '/data/run/FLOW/*.FLOW.hdf5' \
  --hours 6 --end 2026-11-03T12:00:00Z \
  --output /data/lifetime_studies/replay
```

This includes acquisitions starting in `[06:00, 12:00)` UTC on November 3.
Default end time is the actual time the command starts, so old archived inputs
require `--end`; the window never shifts silently to the newest available file.
Filename `CST` and `CDT` suffixes retain their recorded offsets across DST.
Modification times are not acquisition times. For a selected-segment sidecar,
the recorded `source_file` or `source_flow_url` supplies the acquisition filename.

Membership is **whole-file by acquisition start**, not an exact wall-time cut
on individual segments. The common segment schema lacks an absolute event
timestamp; its `t` field is electron drift time. Included acquisitions can
extend past the end boundary, and acquisitions starting before the start
boundary are excluded. The JSON receipt states this scope explicitly.

## Pool already stored segment arrays

The notebook-style gather-and-fit option is explicit:

```bash
python actions/lifetime/recent_tracks.py '/data/run/FLOW/*.FLOW.hdf5' \
  --selection stored --dataset analysis/rock_muon_segments/data \
  --hours 6 --end 2026-11-03T12:00:00Z \
  --output /data/lifetime_studies/stored
```

`--selection stored` labels the noise quality as **unverified upstream
selection**. Numerical `valid_segments` filtering cannot identify hot pixels
or synchronous hardware firing from `dQ`, `dx`, `nhits` and `t` alone.
Without `--dataset`, the command uses both disjoint timing categories when
present, otherwise the legacy dataset. It never additionally includes the
legacy copy when split selections are available. A single available timing
category must be requested explicitly using `--dataset`.

Missing requested datasets are listed as skipped. Empty datasets are counted;
nonfinite rows and rows with nonpositive `dx`, `dQ` or `nhits` are excluded.
Repeated file-local segment IDs from different files are preserved. Repeated
paths, symlinks, hard links and recorded source identities are deduplicated;
exact duplicate nonempty segment payloads are rejected to catch renamed copies.

The pool requires matching recorded geometry, calibration and selector metadata
and matching array schemas. It refuses incompatible inputs rather than mixing
them into one estimate. Acknowledging `--allow-missing-metadata` permits older
incomplete metadata and assumes a 0.1-µs drift tick only when absent. A recorded
different tick period is always rejected. This flag cannot establish equal
unrecorded calibrations, selections or live hardware settings.

## Outputs and interpretation

Each invocation creates a new `recent-*` directory below `--output`:

- `summary.json`: interval, inputs and skips/errors, numerical counts, metadata,
  quality-selection mode, fit status and result;
- `segments.h5`: the valid pooled segments and parallel `member_index`,
  `source_dataset_index` and `source_row` datasets identifying each row in its
  recorded selection file;
- `selections/`: vetted or freshly reconstructed per-file sidecars, when used;
- `fit.png`: the diagnostic figure if the shared fitter accepts the estimate.

An input error prevents fitting a silently incomplete pool. No valid segments,
failed fits and rejected estimates produce a JSON receipt with null lifetime
and process exit status **2**; an accepted fit returns **0**. Argument/discovery
errors also return a nonzero status. Source FLOW files are never modified.

Input errors print their file and reason immediately, and the final summary
prints the blocking message and error count. Completed selection sidecars keep
unique input-index filenames even when a later compatibility check rejects a
file. A broad or ambiguous PPS interval is handled as a recorded io_group-level
segment veto, not an input-file error; independently evaluable groups and
events remain available for the fit.

During fitting, each distinct `pylandau` log message is displayed once at its
original severity. `summary.json` retains every distinct message and its full
occurrence count in `fit_warnings`, so repeated optimizer warnings do not flood
the terminal. The temporary filter is removed even after a fit exception;
SciPy/Python warnings and other loggers are unchanged.

The command calls the same `lifetime.fit_samples` and `langau_lifetime` path as
the existing track analysis, including its usable-bin, covariance and 20-ms
acceptance checks. Bins are exactly:

```python
time_bins = np.linspace(0, 1960, 20)
dqdx_bins = np.linspace(0, 250, 70)
nhits_bins = np.linspace(0, 30, 42)
```

`t` is in **0.1-µs ticks**, so the drift range above is **0–196 µs**, not
0–1960 µs. Charge, length and hit-density inputs retain their existing FLOW
units and conventions. Fit errors are conditional statistical uncertainties;
they exclude selection, calibration and timing systematics. A six-hour pool
also assumes lifetime and acquisition conditions are sufficiently stable during
that interval. Inspect the fit diagnostics and selection losses before treating
the numerical lifetime as a validated argon-purity measurement.
