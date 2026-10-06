# Beam/cosmic electron-lifetime test

The normal nearline configuration keeps using
`rock_muon_selection_data.yaml`, whose output combines all selected events. To
test separate externally triggered and off-beam muon samples, configure the
worker environment before running `actions/flow_charge_data.sh`:

```bash
export ARCUBE_NEARLINE_MUON_WORKFLOW=yamls/proto_nd_flow/workflows/rock_muon_selection_beam_cosmic_data.yaml
```

This workflow must be available in the installed `ndlar_flow` checkout. It
writes these segment datasets:

- `analysis/beam_rock_muon_segments/data`: events with `n_ext_trigs > 0`
- `analysis/cosmic_muon_segments/data`: events with `n_ext_trigs == 0`

These names describe timing categories. An external trigger is not particle
truth, and an off-beam activity trigger is only cosmic-enriched.

`lifetime.py` auto-detects the split datasets and falls back to the legacy
`analysis/rock_muon_segments/data` dataset. New JSON rows carry `sample`, input
file, dataset, segment count, and fit status fields. The time-series plotter
also accepts the historical three-field JSON rows and labels them `mixed`.
Successful fits also receive a diagnostic `fit_quality` flag. It identifies
nonphysical, poorly constrained, and out-of-display-range results but does not
silently remove them; physics-approved alarm and exclusion criteria are still
needed before production use.

To attach each per-file fit to its selected-segment HDF5 group, set:

```bash
export ARCUBE_NEARLINE_WRITE_LIFETIME_METADATA=1
```

This metadata write is opt-in because it reopens an otherwise completed flow
file in read/write mode. Test it on copied flow files before enabling it in a
nearline worker.

## Unified purity monitor

`lifetime.py` is the central command for fitting files, importing slow-controls
snapshots, querying databases, updating history, and rendering the overlay.
`lifetime_timeseries.py` remains a compatibility wrapper. The default PNG is
exactly **3000 × 2000 pixels**; the companion is `<plot>.png.html`, retaining the
existing shifter URL convention. Connecting lines join successive valid points across missing/rejected windows
as visual guides; no measurements are interpolated or filled. Each source's latest measurement date appears below the PNG.

Install the standalone analysis dependencies (Python 3.10 or newer):

```bash
python3 -m venv ../.venv-purity
../.venv-purity/bin/python -m pip install -r actions/lifetime/requirements.txt
```

No FireWorks, MongoDB, MPI, ROOT, or full FLOW installation is needed to read
already-produced FLOW files or make plots. Running `watcher.py` itself still
requires the normal nearline environment; creating missing track selections
requires the `ndlar_flow`/`h5flow` environment documented in that repository.

### One FLOW file and its native packet companion

```bash
../.venv-purity/bin/python actions/lifetime/lifetime.py \
  --input_file /data/packet-0070002-2026_09_29_09_17_12_CDT.FLOW.hdf5 \
  --packet-file /data/packet-0070002-2026_09_29_09_17_12_CDT.h5 \
  --output_file_plot /output/file.lifetime.png \
  --output_file_json /output/json_elifetime.json \
  --output-timeseries /output/elifetime_time_series.png \
  --slow-controls /snapshots/slowcontrols-20260929.json
```

Omit unavailable inputs. `--slow-controls` can be repeated; pass all archived
snapshots that overlap the aggregation windows, not only the latest partial
window. `actions/lifetime.sh` does this automatically for its snapshot directory.
An existing history can be replotted with only `--output_file_json` and
`--output-timeseries`. `--timestamp` supplies an offset-aware timestamp for
nonstandard input names. No source FLOW data are changed unless the explicit
`--write-hdf5-metadata` option is used.

Track samples are fit once per FLOW file. When both split datasets exist, a
third **all selected MIP candidates** fit pools their disjoint segments; it is
not an average of fitted lifetimes. Legacy `rock_muon_segments` supplies the
all-candidate fit only. The existing selection uses long, straight,
through-going candidates; it is not an independent, truth-level identification
of every MIP or of rock versus cosmic origin. The waveform/event timing and
selection definitions remain those in `ndlar_flow`. `Q` in the inspected
`CalibHitBuilder` is not lifetime-corrected (`E` is).

Ahmad's reference fitter is run with `--detector-wide`: charge objects from all
IO groups are pooled before each drift-bin Landau-Gaussian fit. Its physical
channel grouping still occurs before pooling. The central command uses
`ADC - pedestal` and 28-tick successive-hit spacing, matching his README example.
Both choices, and the automatic EXT stream choice, remain **provisional**.
`--ext-io`, `--successive-ticks`, and `--mpv-bootstrap` are explicit controls.
The output retains fit configuration, counts, coverage and diagnostics. No
packet lifetime error bar is invented: the reference fitter supplies bootstrap
MPV errors, but does not estimate an uncertainty on the final lifetime.

The uploaded `packet_lifetime_final` helper was absent. `packet_pedestal.py`
replaces that dependency by fitting the static identity
`Q_raw = slope * (ADC - pedestal)` per physical channel, using explicit
FLOW hit-to-packet references. This avoids assuming a common gain or calibration
epoch. Constant-ADC channels and inconsistent/nonlinear channels cannot provide
an identifiable pedestal and are omitted; inspect the resulting coverage.
The FLOW file must correspond to the native packet file. An explicit
`--pedestal` selects the legacy panel JSON instead; as Ahmad's PDF notes, a
panel mean need not represent the static electronics baseline.

DeMario's existing histogram-convolution observable and drift-bin choices are
preserved. Drift `t` is in 0.1 µs ticks; the existing midpoint division by 20
includes both midpoint averaging and tick conversion. Empty histogram support,
failed slice fits and invalid MPV uncertainties are now failures rather than
initial guesses accepted as measurements. These checks establish software
behavior, not physics approval of the estimator or alarm thresholds.

### Slow-controls data and units

The source configuration contains **no passwords**. Set connection URLs in the
named environment variables using the approved credentials on the export host.
Copy `sources.example.json` to `sources.local.json` and verify its host/database
and timestamp conventions before running live exports.

| Measurement | Source established from repository / supplied Grafana panel | Stored unit |
|---|---|---|
| PRM lifetime | PostgreSQL `prm_table.prm_lifetime` (BlobCraft2x2) | seconds |
| H₂O | PostgreSQL `query_cryo_float(start, end, 1893)` | ppb |
| O₂, ppb channel | `query_cryo_float(start, end, 1890)` | ppb |
| O₂, ppm channel | `query_cryo_float(start, end, 1874)` | ppm |
| N₂ | `query_cryo_float(start, end, 1871)` | ppm |

The example shows separate estimates for both O₂ channels and archives N₂ with
`use_for_lifetime: false`. Either O₂ channel can be excluded with that same flag.
Showing both is not proof that both analyzers are valid during every historical
period. Never combine or silently switch the two O₂ channels. Their units come from the
supplied Grafana field settings and should be checked against calibration and
returned values. Eva's script contains no N₂ conversion coefficient.

BlobCraft2x2 also documents the older O₂ Influx mapping:
`cryo_readonly / AE-1015A / magnitude`, accessed through a DAQ05 tunnel.
The exporter supports `kind: influx`, with `connection`, `connection_env`,
`database`, `measurement`, `field`, `quantity` and `unit` settings. The current
Grafana gas channels can instead use PostgreSQL directly, avoiding that extra
Influx dependency for gases. PRM Influx can be configured if its actual
measurement/field/unit are established; no PRM Influx name is guessed.
`kind: ignition` also supports monthly `table_prefix`, `tag_id` archives, while
`kind: ignition_function` uses the supplied Grafana query interface.

PRM observations are converted to µs and averaged arithmetically in
America/Chicago calendar days. Bars show sample SEM (absent for a single
observation), not instrument uncertainty or day-to-day spread. Gas
concentrations are converted to ppb, averaged in local 00–06, 06–12, 12–18,
18–24 windows, and then converted using Eva's formula:

```
tau_us = 1000 / (mean_O2_ppb / 0.299 + mean_H2O_ppb / 17)
O2_only_tau_us = 299 / mean_O2_ppb
```

The two estimates are distinctly labeled. Missing water is not treated as
zero. Zero concentrations do not yield a fabricated finite lifetime. Gas
measurement/conversion uncertainties are not supplied and no bars are drawn.
Eva's original code explicitly questions the conversion; gas-phase readings
must not be interpreted as a validated liquid-phase equivalence. No linear
projections from the original illustrative script are used for monitoring.

Snapshot schema (one raw observation, not a daily average):

```json
{"schema_version": 1, "measurements": [
  {"timestamp": "2026-09-29T08:15:00-05:00", "quantity": "prm_lifetime",
   "value": 0.0015, "unit": "s", "source": "PSQL PRM"}
]}
```

The numbers above illustrate the schema, not an actual measurement.
`quantity` can be `prm_lifetime`, `o2`, `h2o`; other archived quantities must have
`use_for_lifetime: false`. Optional `valid: false` excludes an invalid reading.
Timestamps must include UTC offsets. Export configuration `naive_timezone`
explicitly interprets database timestamps lacking an offset; the example uses
UTC following BlobCraft2x2, but this needs verification for the live database.
Exact overlapping observations are deduplicated; conflicting values at the
same source/quantity/time fail rather than being silently averaged.

### DAQ05 → NERSC operation

Run the export on a host with database access (e.g. DAQ05 or through an approved
read-only tunnel). PostgreSQL runs in a read-only transaction with a 60 s
statement timeout; range limits are inclusive start/exclusive end. Example:

```bash
python3 actions/lifetime/lifetime.py \
  --export-slow-controls /exports/purity/slowcontrols-20260929.json \
  --source-config actions/lifetime/sources.local.json \
  --start 2026-09-29T00:00:00-05:00 --end 2026-09-30T00:00:00-05:00
```

Use immutable timestamped filenames when exporting successive intervals.
Exports are atomic and a failed source query does not overwrite a good
snapshot. A successful empty query is distinguishable from a query failure in
the snapshot's source counts. Archive raw samples to allow averaging changes
without querying the control systems again.

On NERSC, configure an approved source path and existing SSH route:

```bash
export ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE=user@host:/exports/purity/
export ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR=/global/cfs/cdirs/dune/www/data/2x2/nearline_run3/slow_controls
bash extra/data_mgmt/sync_purity.sh
```

This performs ordinary rsync temporary-file/rename transfers without deletions,
then refreshes the plot even if no FLOW arrives. Schedule it using the site's
existing nearline scheduling mechanism; no remote cron or service is installed
by this change. An alternative is a separate watcher for immutable snapshots:

```bash
./watcher.py actions/lifetime.sh --ext json --path "$ARCUBE_NEARLINE_PURITY_SNAPSHOT_DIR"
./watcher.py actions/lifetime.sh --path /global/cfs/cdirs/dune/www/data/2x2/nearline_run3/flowed_charge
```

The watcher deduplicates by file path, so overwriting one snapshot filename
will **not** enqueue another task. Periodic `sync_purity.sh` refreshes do not
have this restriction. Use `--dry-run` to inspect watcher discovery first.

`actions/lifetime.sh` finds native packets through the existing mirrored
`flowed_charge` → `packet` paths, loads all slow-control snapshots, and calls
only the central command. `ARCUBE_NEARLINE_PYTHON` can select the analysis venv;
`ARCUBE_NEARLINE_PACKET_EXT_IO` and `ARCUBE_NEARLINE_PACKET_PEDESTAL` configure
packet inputs. Total FLOW-file size no longer prevents a small segment read.
The original paths remain `lifetime/jsons/json_elifetime.json` and
`lifetime/plots/elifetime_time_series.png[.html]`.

History updates are locked and atomic. Retries replace a result for the same
source file/sample rather than appending duplicate points. Failed fits remain
as explicit status rows with null lifetimes. PNG/HTML rendering is serialized
with history updates to avoid a slower worker overwriting a newer plot.
Input snapshots and FLOW files should be complete before watcher discovery.

Run regression tests with:

```bash
../.venv-purity/bin/python -m pip install pytest
../.venv-purity/bin/python -m pytest tests/test_lifetime.py -q
```

### Plateau and calibration quality gates

Before a gas lifetime is accepted, each contributing analyzer needs at least
three hours of observations within the six-hour bin, without a sampling gap
over 30 minutes. A flat run lasting at least three hours rejects every
intersecting bin, including the run's initial observations. There is no
interpolation or forward filling. Short initial/current bins remain unavailable
until sufficient history exists. Rejection reasons and input means/counts are
retained under `analyzer_quality` in the history JSON.

`--gas-quality-config gas_quality.example.json` overrides these defaults and
accepts per-source overrides under `sources`, keyed by the exact source name.
`flat_tolerance_ppb` is a peak-to-peak tolerance after unit conversion. Its
conservative default of zero detects exactly repeated digital values; an
analyzer's calibrated resolution is needed to choose a tolerance that also
rejects noisy saturation plateaus. `minimum_ppb` and `maximum_ppb` can enforce
known calibration limits. There is no assumed upper limit by default.
A stable reading is not itself proof of instrument failure; these are explicit
monitoring exclusion rules requested for the overlay, not calibration findings.
If H₂O fails but O₂ passes, the combined result is withheld while the separately
labeled O₂-only result may still be shown.

### Files with no selected muons

`--select-muons` runs the installed `ndlar_flow` `RockMuonSelection` methods
against read-only event→hit references and the geometry stored in the input.
It writes a small `.segments.h5` sidecar, not a copy or modification of the
multi-GB FLOW file. All events are processed. The output stores event/track/
segment counts and a SHA-256 of the selector implementation. This requires
`ndlar_flow`, `h5flow`, and scikit-learn in the Python environment. Set
`ARCUBE_NEARLINE_SELECT_MUONS=1` to enable this fallback in the watcher action.
An empty sample is a recorded failure, not zero lifetime.

The internal sample key `beam` is retained for compatibility, but the displayed
label is **externally triggered**: in a cosmic commissioning run, an external
trigger need not be a beam spill. Physical rock/cosmic separation needs an
independently verified beam timing/trigger convention. The combined fit and
externally triggered fit will coincide if every accepted candidate belongs to
that category.

### Live source verified during local validation

DAQ05's accessible cryogenics reader configuration uses
`ifdb11.fnal.gov:5461`, database `argoncube2x2_prd`; the older `ifdb08` name did
not resolve. The supplied jump route was verified with read-only queries:
`ssh -J acdcs@acd-gw06.fnal.gov acdcs@acd-daq05.fnal.gov`.
For export on that host, a source may specify `credential_ini` and optional
`credential_section` instead of `url_env`. For the existing installation these
are `/home/acd/acdcs/2x2/SlowControls2x2/Cryogenics/config.ini` and `secrets`.
Credentials are consumed in memory, not included in the snapshot. This does not
install or run anything remotely by itself.

Month-wide queries hit the existing 60-second timeout during testing. The
exporter therefore queries `query_cryo_float` sequentially in at most six-hour
windows (`query_hours` is configurable), with no parallel database queries.
The source's inclusive end is filtered back to half-open windows locally.

### Local validation on 2026-10-05

The smallest available commissioning FLOW file,
`packet-0070002-2026_09_29_09_17_12_CDT.FLOW.hdf5` (6,724,434,264 bytes),
and its native packet file (787,334,869 bytes) were downloaded from the Run 3
NERSC portal. Processing all 3,851 events produced 15 selected externally
triggered tracks and 746 segments, with **1.632 ± 0.149 ms** fitted lifetime.
The combined selection gives the same result because no off-beam candidates
passed these cuts. The source FLOW file was opened read-only.

The detector-wide packet attempt recovered pedestals for 81.1% of selected
packets and formed 7,379 fit objects, but only one drift bin passed the
reference MPV quality requirements. It therefore produced **no accepted packet
lifetime**. The failed status and diagnostic files are retained; the selection
was not relaxed to manufacture a point.

Read-only DAQ05 queries exported 49,393 PRM/analyzer observations from
2026-09-28 onward. The extended local overlay includes 26 daily PRM means from September 5 onward, both
oxygen channels separately, the O2-only and O2+H2O conversions, and track
lifetimes independently recalculated from 22 FLOW files. Published numerical
lifetimes are retained only for corroboration, never as plotted input. Short boundary windows
were withheld. Accepted windows did not trigger the default exact-value
plateau check; noisy plateau detection still requires an instrument-appropriate
`flat_tolerance_ppb` and calibration limits.

The regression suite passed 30 tests, including plateau rejection, calibration
bounds, unit conversion, DST boundaries, duplicate handling, concurrent history
writers, packet pedestal references, and a PRM-only CLI refresh. The real-data
PNG was checked at exactly 3000×2000 pixels, and the watcher action was run
locally without muon data. The environment passed `pip check`. Artifacts and
inputs are under the sibling `purity_validation/` directory and are not Git
inputs. No NERSC scheduler, DAQ service, or GitHub push was changed by validation.

### Six-hour track summaries and connecting lines

The overlay now displays an arithmetic mean of accepted per-file track
lifetimes in each local 00–06, 06–12, 12–18, and 18–24 window, separately for
each selection. Each file is fitted directly from its FLOW segment arrays.
The six-hour summary weights those file estimates equally; it does not re-fit
segments pooled across files or weight by exposure.
The original per-file results remain under `lifetimes`; derived displayed
rows are stored under `plot_lifetimes` with window boundaries, member
provenance, accepted/rejected counts, and uncertainty components.

Track error bars are the larger of the propagated independent fit errors
`sqrt(sum(sigma_i**2))/n` and the between-file SEM, where available. This keeps
large observed file-to-file scatter visible; shared calibration systematics
and correlations are not estimated. A singleton retains its original fit
error. A failed-only window remains unavailable. These averages describe
files with successful fits, not guaranteed coverage of an entire six hours.

Lines match their markers. They are guides between observed estimates, not
measurements between files. All valid points are connected, including across
missing or rejected windows. Invalid measurements remain unplotted. O2-only estimates use dashed lines and crosses;
O2+H2O uses solid lines and squares. Tag 1874 is purple with open squares/plus
markers; tag 1890 is green with filled squares/x markers. The combined legend
explicitly names the additional H2O tag 1893. Both conversions remain
provisional and are not independent measurements of purity.


### September 29 beam check and published-file provenance

The 21 reference timestamps were matched one-to-one to distinct NERSC FLOW
files, whose segment arrays were read directly and independently fitted.
Seventeen are hot-pixel-hunt trials on October 2 (12:09:44–14:39:22 CDT);
their independently calculated six-hour mean is 1.489 ms, with 0.084 ms
between-file SEM. The October 5 inputs are explicitly named induced-noise
runs. These means combine successful file fits across the stated windows;
they are not evidence that detector settings were constant or that the
observed variation was entirely a purity change.

The newly processed September 29 FLOW file spans 09:17:11.606–09:18:07.808
CDT and has 1,221 events with external triggers. Its 15 accepted tracks
provide the 1.632 ± 0.149 ms result. The `beam` sample key means only
`n_ext_trigs > 0`, as documented by `RockMuonSelection`; it does not imply an
IFBeam match or physical rock-muon origin.

Read-only IFBeam queries returned HTTP 200 with zero NuMI `$A9` records for
all 15 real devices used by `beam_quality/get_data.cpp` throughout September
29 (local midnight to midnight). The existing macro was run unchanged. A
July 12, 2024 00:01–00:02 CDT control returned 46 spills, all passing its
cuts 1–5. September 29 therefore has **no IFBeam-confirmed beam** in this
check. With no spills, beam quality and detector-to-spill coincidence are
not evaluable; missing archive rows are not measured zero POT. The existing
horn cut is FHC-specific and requires polarity validation for future runs.

The raw queries, hashes, FLOW timing summary, macro logs and control outputs
are retained under sibling `purity_validation/ifbeam/`, including the
reproducible `check_september29.py` and `september29_beam_check.json`.
No detector or control-system state was changed. The revised overlay passed
30 regression tests and was visually checked at 3000×2000 pixels.


### Direct FLOW inputs and reference-only corroboration

Track history must identify a direct FLOW fit. Existing unprovenanced track
records are migrated into `reference_lifetimes`, excluded from all means and
plots. `--reference-history FILE` imports a published JSON only into that
reference section. `reference_comparisons` matches exact timestamps and the
corresponding selection, retains unmatched/ambiguous cases, and reports
per-file differences without replacing or adjusting the calculated result.
The original per-file fits remain in `lifetimes`; six-hour summaries remain
in `plot_lifetimes`.

`--input_file` accepts a local FLOW path or an HTTP(S) FLOW URL. Remote inputs
use conditional HTTP byte-range reads with a strong ETag and an in-memory
block cache. Only the selected segment datasets needed by the fit are read;
HDF5 lifetime attributes and published plots/JSON never supply fit values.
The source URL, ETag, full-file size, dataset shapes/dtypes and SHA-256 hashes
are retained with each result. `requested_range_bytes` counts requested
cache blocks (an upper bound on response bytes when the final block is clipped).
No remote FLOW files are saved by default. `--flow-cache DIRECTORY` optionally
retains exact segment arrays for subsequent offline fits.

```bash
../.venv-purity/bin/python actions/lifetime/lifetime.py \
  --input_file "https://portal.nersc.gov/project/dune/data/2x2/nearline_run3/flowed_charge/ColdCommissioning/20260928_Cosmic/packet-0070002-2026_09_29_07_10_26_CDT.FLOW.hdf5" \
  --output_file_plot /output/file.lifetime.png \
  --output_file_json /output/history.json \
  --output-timeseries /output/elifetime_time_series.png \
  --reference-history /references/published_json_elifetime.json
```

All 21 files behind the published reference were independently refitted via
streaming: all produced finite results, in about 145 seconds, with about
311 MB of requested ranges versus 427 GB of full files. The fits were made
from the existing FLOW segment observables, not by repeating hit-level track
finding in those 21 files. The separate downloaded September 29 09:17 file
had no segments, so its hit-level selection was run locally as documented
above. Remote files without selected segments require local full-file
selection; missing results are never filled from a reference.

The recalculated values are not uniformly identical to the reference: the
largest difference is September 29 07:10, 1.424 ms recalculated versus 1.829 ms
published (−22.2%). Repeating the historical `dx != 0` mask on the same arrays
gives 1.424 ms as well, so the added validity filtering does not explain that
discrepancy. The published environment and exact input generation are not
available here. Preserve and review the differences; do not treat the JSON
comparison as validation of the physics estimator or force agreement.


### September 5 extension and Eva comparison

An additional read-only export recovered 594,790 observations from September
5 to September 28, then combined them with the later snapshots. There are
26 days with positive recorded PRM lifetimes. Missing days are not filled;
lines nevertheless connect successive valid means as requested. The local
PNG and interactive companion cover September 5 through the latest available
October measurements at 3000×2000 pixels for the PNG.

The entire public Run 3 FLOW tree was inventoried: 330 completed files, the
first dated September 28 at 17:47:37 CDT. All 62 September 28 files were
checked for selected-segment datasets and none contained them. No earlier
track points are inferred from a PNG or JSON. The 70 `.hdf5.tmp` files in
the inventory are excluded as incomplete outputs. Earlier track monitoring
requires additional completed FLOW data or full hit-level selection.

`purity_validation/output/eva_corroboration.csv` compares the available
September 4–10 source-script values with our database-derived daily PRM means
and six-hour gas estimates. September 5 gives 66.546 µs versus Eva's 66.5 µs;
September 6 gives 155.699 µs versus 156 µs. September 10 differs (153.400 µs
daily mean versus the script's 215 µs point), because these are different
sampling/aggregation choices. The newer month-long PNG source data are not
available, so it is a visual comparison, not a numeric reconstruction.

The user supplied Brandon How's September 17 10:32 ELog entry: the DF-560
range changed from **0–10 to 0–1 ppm**, with Ignition already reflecting the
change. The event is annotated at 10:32 CDT; no extra factor-of-ten scaling
is applied to the archive. This corrects the `ppb` wording in the reference
image's annotation without guessing a different historical concentration.

Eva-style grey shading marks the approximate September 17–28 PRM signal-
concern interval visible in her image. It is reference context, not a data-
quality mask. The PRM table supplies peak/time diagnostics but no explicit
quality flag, and a validated waveform cut is not available. Positive
recorded lifetimes are therefore retained with daily SEM; especially within
that interval, they should not be assumed to be signal-validated measurements.
No numeric point was digitized from the image.

Use `--annotations FILE` for sourced timestamp/label entries and optional
`end` spans. Annotations are retained in history and survive later refreshes;
both PNG and HTML show them. The range-change event and approximate shading
are stored in `purity_validation/output/event_annotations.json`. The related
`september_extension_summary.json`, `flow_corroboration.csv`, and raw query
snapshots document the comparison and its limits.
