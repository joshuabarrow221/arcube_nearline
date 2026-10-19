# NERSC Codex handoff: 2×2 lifetime monitoring

Prepared for Josh Barrow, 19 October 2026. This is an engineering and scientific
handoff for a new Codex session on NERSC. It supplements the concise
[operator guide](lifetime-integration.md), [methods](lifetime-methods.md), and
[shifter README](../actions/lifetime/README.md); it is not the shifter landing page.

## 1. Objective and current deployment status

Integrate direct FLOW track fits, detector-wide native-packet fits, daily PRM
measurements, and six-hour gas-analyzer lifetime equivalents into one recent
monitoring plot. Preserve these served names:

- `elifetime_time_series.png` — 3000 × 2000 static plot;
- `elifetime_time_series.png.html` — self-contained interactive counterpart;
- `elifetime_diagnostics.html` — failures and source/fit details;
- `elifetime_status.json` — machine-readable estimates, membership and provenance.

Target page:
<https://portal.nersc.gov/project/dune/data/2x2/nearline_run3/lifetime/plots/elifetime_time_series.png.html>

The implementation has been exercised locally on actual source data and compact
captures. **The running NERSC checkout, worker environment, DAQ05 transfer schedule,
and served-directory behavior have not been inspected or changed by this session.**
Local correctness is not evidence that the production service has been upgraded.
Use a preview publication directory first. Inspect the existing writer and stop
it before promoting the replacement; two unrelated publishers must not compete.

Work is on `dev/purity-overlay` in
<https://github.com/joshuabarrow221/arcube_nearline>. Do not assume `main` contains
the monitor. Establish the current remote and checked-out commits at deployment.
Important implementation commits already prepared before this handoff are:

| Commit | Purpose |
|---|---|
| `61655dc` | Configured watcher, pooled monitoring, snapshot ingestion/publication |
| `174c976` | Correct PRM source wall clock; UTC shifter display |
| `17c74b7` | UTC detailed plots and gas midpoint verification |
| `a82ec79` | Include parent-track settings and schema in pool compatibility |

The handoff's own commit will follow these. The local refresh's `delivery.json`
records the exact final local/remote commit and package checksums after push.

## 2. Decisions that should survive the move

1. **Source detector lifetimes from FLOW/native data.** Published PNG/JSON are
   corroboration only. Never digitize their points or substitute their fitted
   lifetimes into our calculation.
2. **Keep PRM daily.** Use six-hour pooled detector observations in the simplified
   monitor and six-hour concentration means for gases. Windows were fixed before
   fitting; do not choose durations or cuts by agreement with another method.
3. **Keep packets in the monitoring plot.** Show finite, positive failed-check
   packet candidates as open markers with diagnostic reasons. Do not invent a
   value for undefined, nonpositive, singular or boundary-constrained fits.
4. The short track labels are **All tracks** and **Cosmic-enriched tracks**.
   The externally triggered subset remains diagnostic rather than another nearly
   redundant monitoring curve.
5. Use true **UTC display axes** everywhere. Averaging boundaries remain Chicago
   calendar windows; changing a display clock must not change membership.
6. Gas values belong at the **exact midpoint of their six-hour window**. PRM
   and pooled detector points use observed-time centroids in the simplified view.
7. Use matching colored points/lines, no translucent gray box. Lines connect
   available accepted values across gaps; they do not imply continuous exposure.
8. Keep failures and all prospectively chosen files in the evidence. Do not
   remove files merely because their resulting lifetime is high or low.
9. Keep source files read-only and retain array hashes, calibration provenance,
   file membership, failed status and uncertainty meaning.

## 3. Local evidence and transfer package

The previous workspace was:

```text
/mnt/c/Users/inuya/OneDrive/UMN-Postdoc/Research/DUNE/NDPrototype_2x2/
  arcube_nearline/                    implementation checkout
  ndlar_flow/                        local upstream context/selection checkout
  BlobCraft2x2/                      original slow-controls/database context
  .venv-purity/                      tested local Python environment
  nearline_lifetime_integration/     explicit code overlay, patch and manifest
  purity_validation/
    pooling_week_20261002_08/        compact, receipted FLOW selections
    refresh_20261009/                corrected PRM history and full fit replay
    refresh_20261014/                six newer FLOW/native pairs and PR253 review
    refresh_20261019/                latest controls, replay, plots and archive
/mnt/d/VMWare/2x2_work/               mounted external drive; older native inputs
```

These are provenance locations, not NERSC path defaults. The local virtualenv
and WSL lock setting must not be copied into a multi-host NERSC environment.

The October 19 `slowcontrols_20260905_20261019.zip` contains the corrected raw
observation history, a manifest, source/quality/conversion configurations, a
checksum-verification helper, and this handoff. Its `snapshots/` directory has
45 UTC-date partitions. A partition is a storage choice, **not** an averaging
window. Import all of them together so Chicago day boundaries and plateau runs
can cross file boundaries correctly.

The canonical archive contains **746,708 distinct observations** from the nine
reviewed input exports; **50,489 repeated observations** were deduplicated only
after checking values, units and analysis flags for agreement:

| Source | Raw observations | Unit retained in archive |
|---|---:|---|
| PRM | 265 | seconds |
| O₂ tag 1890 | 495,145 | ppb |
| O₂ tag 1874 | 67,196 | ppm |
| H₂O tag 1893 | 106,287 | ppb |
| N₂ tag 1871 | 77,815 | ppm; archived, not converted |

Seven nonpositive PRM readings remain in raw evidence and are rejected by the
existing analysis loader. N₂ also retains its source values. The archive does
not contain credentials, native detector HDF5 files, detector shards or fitted
published NERSC histories. It is sufficient to reconstruct PRM/gas summaries,
not historical detector fits.

Only corrected historical snapshots are admitted. Older local exports labeled
PRM wall time as UTC; **do not add those superseded exports to the new archive**.
Exact-observation deduplication cannot repair a wrong timestamp: it would count
the same physical reading twice at different instants. Input-file hashes and
observation counts are retained in `manifest.json`.

### Latest successful database read

The read-only query ended at **2026-10-19 21:16:01.622770 UTC**. PRM was re-queried
from September 5, and gas from October 12 with overlap against the previous
export. Credentials were used in place on DAQ05, never copied into the archive.

- Latest PRM record: **2026-10-16 14:17:18 UTC / 09:17:18 CDT**.
- Latest H₂O: **2026-10-19 21:15:43.697 UTC**.
- Latest O₂ 1890: **2026-10-19 21:15:06.796 UTC**.
- Latest O₂ 1874: **2026-10-19 21:15:56.864 UTC**.
- Latest N₂: **2026-10-19 21:15:33.833 UTC**.

The PRM query succeeded but returned no newer record in the queried table.
This does not establish that the instrument itself stopped taking data; check
the producer and database ingestion if the operator expects later measurements.
Fresh page-generation time must not be presented as fresh PRM measurement time.

The October 16 PRM mean is **247.713 ± 15.649 µs SEM** from three positive
readings. In the latest partial gas window (October 19, 17–23 UTC, observed
through about 21:15 UTC), raw concentration means are O₂ 1890 **145.892 ppb**,
O₂ 1874 **10.843 ppb** after ppm normalization, and H₂O **642.576 ppb**. The
unchanged conversion gives combined equivalents of **1.648 µs** and **9.957 µs**,
respectively. These pass current software checks, but no verified instrument
upper calibration limits or gas-to-liquid correction have been supplied. Do
not promote the numerical result to a validated LAr lifetime. Inspect the
updated gas concentration plots and source operating conditions before drawing
a physical conclusion; no thresholds were tuned to suppress these readings.
The fresh raw export also contains an O₂ 1874 spike to **31.03046 ppm** at
**October 15, 20:26:26.628 UTC**. Its contribution remains in the affected mean
under the unchanged policy. Its cause and valid instrument range are unverified;
the gas concentration companion uses an explicitly logarithmic 1874 panel so
the full spike and ordinary approximately 10-ppb readings remain visible.

## 4. Actual execution path and code map

```text
completed FLOW + native companions          received raw slow-control snapshots
                    \                      /
                     watcher.py discovery
                              |
                     FireWorks ScriptTask
                              |
                     actions/lifetime.sh
                              |
            actions/lifetime/lifetime.py --monitor-config ...
                              |
             monitor.py ingestion, fixed pools, publication
```

`watcher.py` schedules work; it does not run expensive fits in the polling loop.
The wrapper and complete `actions/lifetime/` package are required. Copying only
`watcher.py` or only `lifetime.py` is insufficient. Configured mode is opt-in;
omitting the monitor configuration retains the historical command-line path.

| File | Responsibility |
|---|---|
| `watcher.py` | Detect completed input generations; queue/revisit dependencies |
| `actions/lifetime.sh` | Forward input and monitor profile to central Python |
| `monitor_config.py` | Path resolution, paired filenames, discovery revisions |
| `monitor.py` | Serialized ingestion, active state, fitting and publication |
| `track_shards.py` | Read stored selections, partition, fingerprint and pool |
| `lifetime_funcs.py` | DeMario-derived slice and attenuation fitting |
| `packet_lifetime.py` | Detector-wide native charge-object extraction |
| `packet_pooling.py` | Compatible object pooling, bootstrap and quality gates |
| `calibration_pedestal.py`, `packet_pedestal.py` | Static/linked calibration checks |
| `purity_sources.py` | Read-only export, import/deduplication, means, conversion |
| `monitor_plot.py` | Simplified PNG/HTML and diagnostics |
| `lifetime.py` | Central CLI and separate detailed-study rendering interface |
| `lifetime_io.py` | Publication permissions and locking backend |
| `extra/data_mgmt/export_purity.sh` | Database-host snapshot producer |
| `extra/data_mgmt/sync_purity.sh` | Receiving-host rsync and monitor refresh |

FLOW/native pairing preserves the **complete relative hierarchy**, replacing
`packet-NAME.FLOW.hdf5` with `packet-NAME.h5`. Basename equality across different
runs is not identity. A late native companion triggers a new revision of the
same FLOW input. Its active contribution is replaced rather than counted twice.
Packet-first arrival waits for FLOW. A snapshot can update the plot with no
new FLOW input. No-input manual refresh uses existing state and snapshots;
it does not discover/backfill unprocessed FLOW files by itself.

## 5. NERSC startup procedure

### Establish a separate preview

Inspect existing repository/worker configuration and active publishers first.
Use an isolated checkout or reviewed branch update; preserve site-specific
uncommitted changes. Record the exact commit before starting a long backfill.
The commands below use placeholders that the NERSC operator must set to actual
paths; they are not claims that these directories already exist.

```bash
# In the reviewed checkout of dev/purity-overlay:
git status --short --branch
git log -5 --oneline

# Select/create the site's Python 3.10+ analysis environment first.
python3 -m pip install -r actions/lifetime/requirements.txt
# FireWorks is also needed for queueing; retain the site's configured version
# and LaunchPad credentials. Do not reset a populated FireWorks database.

export ARCUBE_LIFETIME_FLOW_ROOT=/actual/nearline_run3/flowed_charge
export ARCUBE_LIFETIME_PACKET_ROOT=/actual/nearline_run3/packet
export ARCUBE_LIFETIME_SNAPSHOT_ROOT=/actual/private/purity/snapshots
export ARCUBE_LIFETIME_STATE_ROOT=/actual/private/purity/state
export ARCUBE_LIFETIME_PREVIEW_ROOT=/actual/served/purity-preview/plots
export ARCUBE_NEARLINE_PYTHON=/actual/analysis/environment/bin/python
export ARCUBE_NEARLINE_LIFETIME_CONFIG="$(realpath config/lifetime.example.json)"
unset ARCUBE_NEARLINE_LOCAL_OUTPUT
```

Watcher and workers must resolve the same paths/environment variables. A
site-specific profile with absolute paths is often easier to audit. Set
`display_timezone` to `UTC`, `timezone` to `America/Chicago`, and choose
`display_days` (local simplified plots use 7; example profile defaults to 35).
Use separate state/publication trees, outside every watched tree. If needed,
set `calibration_root` to relocate the exact FLOW-named pedestal basename.
Do not substitute a convenient pedestal file merely because its name is similar.

### Import the corrected history and smoke-test

Extract the slow-controls ZIP, run its `verify_archive.py`, then copy only the
completed `snapshots/slowcontrols-*.json` into the configured flat snapshot
directory. Keep manifest/config/handoff files outside that watched directory.
Do not use `rsync --inplace` or `--delete` against the snapshot archive.

```bash
# Read-only scheduling preview; no MongoDB connection.
python3 watcher.py actions/lifetime.sh \
  --lifetime-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG" --dry-run --once

# Controls-only publication from the imported archive and any existing state.
bash actions/lifetime.sh --monitor-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG"

# Test one actual completed FLOW/native pair and its named calibration.
bash actions/lifetime.sh /actual/flowed_charge/run/packet-NAME.FLOW.hdf5 \
  --monitor-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG"
```

Check the preview's PNG, interactive HTML, diagnostic failures, source freshness,
UTC times and source membership. A passing process can still contain rejected
scientific fits. Then test the actual queue/worker path with the site's existing
LaunchPad configuration:

```bash
python3 watcher.py actions/lifetime.sh \
  --lifetime-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG" --min-file-age 60
# In the configured worker environment:
rlaunch rapidfire --nlaunches infinite
```

One watcher per configuration is sufficient. Source producers should close and
rename completed files into place. A minimum modification age is only a heuristic,
not proof that a producer has closed the file. Size worker RAM for complete
native timing arrays. Updates to a state namespace are serialized; many workers
do not parallelize the same growing pool. Use worker wall limits of at most six
hours, below the documented 24-hour shared-filesystem lock lease.

### Promote and recover

After the NERSC preview has passed, arrange the publication directory and stop
the old lifetime publisher. Preserve the public basenames so shifter URLs remain
valid. Verify symlink serving and group permissions under the real web server.
Each new generation is rendered under `.lifetime-generations/`, and
`.lifetime-current` switches atomically. Stable PNG/HTML/status aliases follow
that link. Initial replacement of legacy regular files is sequential, with
backups in `.previous/`; subsequent generation switches are atomic.

For rollback, stop current writers and switch to a retained generation or restore
the old files/publisher. Inspect failed tasks rather than repeatedly resubmitting
unchanged bad input. After repairing an extraction/IO failure, use the normal
FireWorks rerun command, e.g. `lpad rerun_fws -i <fw_id>`. Scientific fit failures
remain diagnostic records. If intentional reprocessing preserves file size and
mtime, change the profile's `reprocess_token`. Do not mutate calibrations in place.

## 6. DAQ05 and future automatic source updates

**The watcher does not SSH to DAQ05 or query a database.** Export, transfer and
watching are separate implemented steps whose schedules/access remain site work.
The local read-only route authorized by the user was:

```bash
ssh -J acdcs@acd-gw06.fnal.gov acdcs@acd-daq05.fnal.gov
```

Treat DAQ05 as a live acquisition host. Read-only SQL and a narrow Python process
were used; no services, DAQ configurations or database contents were changed.
Do not copy credentials into Git, JSON snapshots, CLI arguments or public logs.
The local exporter read the existing credential file in place on DAQ05; the
production wrapper instead resolves credentials through protected environment
variables referenced by `sources.example.json`.

Verified source contract:

| Quantity | Query handle | Clock/unit |
|---|---|---|
| PRM | `prm_table.timestamp`, `prm_lifetime` | timezone-free Chicago wall time; seconds |
| H₂O | `query_cryo_float(..., 1893)` | epoch milliseconds; ppb |
| O₂ low-range channel | `query_cryo_float(..., 1890)` | epoch milliseconds; ppb |
| O₂ second channel | `query_cryo_float(..., 1874)` | epoch milliseconds; ppm |
| N₂ | `query_cryo_float(..., 1871)` | epoch milliseconds; ppm; not converted |

Names/units are based on the operator's Grafana panel definitions and actual
query results. A legacy Influx adapter exists but was not the source of this
successful PRM export. Do not assume Influx must be configured for this table.

The PRM source clock is the operator's stated LabView/Windows convention, not an
independent audit of that program. Localize naive values in `America/Chicago`
then convert to UTC; convert query bounds into the same source wall clock before
dropping timezone information. Do not hard-code +5 hours: winter uses +6. The
production adapter rejects ambiguous/nonexistent DST wall times rather than
guessing. Ignition epoch instants receive **no** PRM-style correction.

On the database-reachable export host, configure `PURITY_PSQL_URL` and
`PURITY_IGNITION_URL` securely, then:

```bash
bash extra/data_mgmt/export_purity.sh /actual/sources.json \
  "$START_ISO" "$END_ISO" /actual/purity_exports
```

Use explicit timestamp offsets. Export overlapping intervals (e.g. a rolling
48-hour query) so late arrivals and cross-boundary plateau evidence are retained.
The exporter publishes only after every enabled query succeeds. A successful
zero-row query is distinct from a failed query or an instrument zero.

On the receiver, with the monitor environment loaded:

```bash
export ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE='user@host:/actual/purity_exports/'
# Configure RSYNC_RSH for the approved site route if needed.
bash extra/data_mgmt/sync_purity.sh
```

The receiving wrapper copies completed snapshots and refreshes the monitor.
Watcher discovery can also queue the same idempotent update. Agree on an actual
export/transfer schedule and execution host; neither script installs one. Keep
older snapshots and their manifests for audit. Future overlapping exports are
deduplicated; conflicting same-source/time readings fail rather than silently
choosing the latest value.

## 7. Scientific definitions and errors

### Tracks and external triggers

“All tracks” means **available upstream-selected through-going MIP candidates**,
not every reconstructed track or a pure particle-identity sample. The inclusive
fit uses their segment union. “Cosmic-enriched” uses the same selection restricted
to parent events with `n_ext_trigs == 0`; it is an overlapping subset, hence its
errors are correlated with All tracks. The legacy `beam` key means
`n_ext_trigs > 0`, not verified beam origin. External timing/trigger records can
occur without neutrino beam; no IFBeam/POT/quality cut is applied by this label.
Do not interpret an externally triggered September 29 point as evidence of beam.

The configured monitor consumes stored FLOW selections. Explicit associations
partition legacy segments when available; missing associations permit only the
stored combined selection. No downstream method can restore candidates excluded
by the upstream selector. `actions/flow_charge_data.sh` supports
`ARCUBE_NEARLINE_MUON_WORKFLOW`, but retains its existing workflow by default.
Installing this repository does not automatically install/change ndlar_flow.

Pool segments before fitting; do not average fitted lifetimes in the monitoring
view. Retain full-directory acquisition cohorts, including Step/trial distinctions,
and reconstruction/calibration/selection compatibility. Commit `a82ec79` adds
parent-track attributes and schema to the fingerprint because unified selection
parameters can live there rather than on segment outputs. Matching absent
metadata does not establish equal code versions or stable detector settings.

DeMario's existing estimator fits 19 drift slices over 0–1960 ticks (0.1 µs/tick).
Each needs at least five valid segments and supported fixed charge/hit-density
histograms. Peak-fit prerequisites include bins outside the final attenuation
subset. An early unsupported slice can therefore abort a fit even if later
drift bins look populated. The error is from SciPy `curve_fit` covariance, **not
MINUIT**: slice covariance includes its residual scaling, then attenuation uses
MPV errors with `absolute_sigma=True`; `sigma_tau = sqrt(C_tau,tau)`. It excludes
shared calibration, event/track correlations and selection/model systematics.

### Native packets

The monitor requires matched native `packet-*.h5`, FLOW and the exact static
pedestal JSON named by FLOW. It verifies linked charge/ADC relationships,
groups signed `ADC - pedestal` by physical channel using the inherited 28-tick
successive-packet rule, and fits charge objects relative to EXT (default IO6).
This is one detector-wide object fit, not an average of TPC lifetimes. Check
trigger wiring and calibration provenance on the actual site setup.

Prospective gates: at least 50 objects per 10-µs slice; positive Landau–Gaussian
MPV and slice chi-square/ndf <3; at least 80% successful bootstrap replicas and
at least 10 successes; at least eight accepted 20–180-µs fit slices spanning
100 µs; interior/nonsingular attenuation fit with chi-square/ndf ≤3 and inverse
lifetime resolved beyond one covariance sigma; pedestal coverage ≥95% in **each**
member file. These are monitoring checks in our implementation; do not describe
every threshold as a requirement established by Ahmad's presentation.

The current alpha bound is [-5,5] ms⁻¹, implying an accepted positive interior
lifetime **above 200 µs**. A boundary result is unavailable, not an accepted
200-µs point. The October 14 afternoon pool hit this bound and had poor model
agreement (diagnostic chi-square/ndf about 20). More pooling does not remove the
bound. A scientifically justified extension is separate estimator-development
work; do not relax it merely to match PRM or track values.

Use 100 fixed-seed object-bootstrap replicas per slice. MPV uncertainty is the
sample SD of successful replica MPVs, **not** SEM of those replicas. The weighted
attenuation Jacobian gives `pinv(J.T @ J)`, without reduced-chi-square scaling;
`tau_us = 1000/alpha` and `sigma_tau_us = 1000*sigma_alpha/alpha**2`. Do not divide
again by the number of files. Shared EXT/event/channel/calibration effects are
not covered by independent-object resampling. Failed candidates have errors only
when a meaningful covariance survives; an absent bar never means zero error.

### PRM and gas equivalents

Daily PRM is the arithmetic mean of distinct positive finite recorded lifetimes.
Do not average inverse lifetimes and invert back. Its error is sample SD / sqrt(N)
for N>1; a singleton has no estimated error. Individual PRM readings have no
provided waveform-fit errors, so do not attach the day's SEM to each reading.

For gas, evaluate quality over the entire retained time sequence before grouping.
Defaults require at least three hours' span, no gap over 30 minutes, nonnegative
values, and no run of exactly repeated readings lasting at least three hours.
Runs can cross export/window boundaries; gaps reset them. Tolerance defaults to
**zero**, so a nearly flat noisy/saturated signal can pass. No instrument-specific
upper calibration limit has been validated. Review resolution/scale/range and
sampling-path knowledge before tightening `flat_tolerance_ppb` or bounds; do not
claim that the current gate detects every calibration plateau.

Average concentrations over six hours, then convert once:

```text
tau_us = 1 / (k_O2 * mean_O2_ppb + k_H2O * mean_H2O_ppb)
k_O2  = 0.0037505851298700343 us^-1 ppb^-1
k_H2O = 0.000093             us^-1 ppb^-1
```

These are Eva's `eva_20261006` model at 212.125 V/cm reference field. The two O₂
channels stay separate. There is no verified zero-offset or gas-to-liquid phase
correction, nor an instrument/conversion uncertainty model; label the results
provisional equivalents. Missing or rejected H₂O is not zero H₂O. An O₂-only
diagnostic exists separately. N₂ is archived with `use_for_lifetime=false`.

The operator's ELog quotation says the DF-560 scale changed **0–10 to 0–1 ppm**
on September 17 at 10:32; an older image's “ppb” annotation conflicts with this.
Do not infer a retrospective calibration correction from that image. The two
O₂ channels also need not track one another: the October 14 comparison showed
1890's combined equivalent falling from about 131 to 28 µs while 1874 stayed
around 31 µs. A flat equivalent alone does not prove saturation or constant LAr
purity; inspect the raw concentration companion and channel provenance.

## 8. Three plot products are deliberately different

| Local product | Detector estimator | PRM/gas |
|---|---|---|
| `elifetime_averages.png` | Six-hour means of per-file track fits; 24-hour packet-object pools | Daily PRM; six-hour gas |
| `elifetime_averages_raw.png` | Same detailed estimates plus individual track/PRM readings; lines follow means | Same summaries |
| `elifetime_time_series.png` | Six-hour pooled segment/object fits; simplified recent window | Daily PRM; six-hour gas |

The first two are expert study views. The configured watcher generates the third
plus its diagnostics/status; it does not automatically regenerate both legacy
study views. Their histories, packet periods and track estimators are not
interchangeable. The local refresh scripts reproduce all three explicitly.

The simplified track/packet points use mean member-file start time; a nominal
six-hour pool may contain minutes of data. PRM uses mean observed reading time.
Detailed averages use nominal centers; their raw companion adjusts track/PRM
means to observed-time centroids but retains packet-window centers. Thus older
detailed MIP/packet points can be horizontally displaced even for the same
files. This display inconsistency was explained but not silently changed. It
is a future presentation choice, not a detector time delay. Gas centers remain
fixed in every view, as explicitly requested.

## 9. Upstream selection and unresolved comparison

The local review of <https://github.com/DUNE/ndlar_flow/pull/253> on October 14
used head `350979656b63fc35b0d90aaa5b19e41444329773`; it was open/unmerged then.
**Recheck its current state on NERSC.** The installed local ndlar_flow was not
automatically changed. The review covered x/z boundary correction, YAML selector
parameters, crossing tags, PCA components and speed changes.

All 18,812 events in six fixed October 14 files were reselected with that pinned
source for comparison. This tested clustering/selection/segmentation through a
segment-sidecar adapter, not the full H5Flow stage or tag serialization. Morning
successful lifetimes matched while some afternoon segment counts changed; this
limited agreement does not establish general physics equivalence. Main plots
continue to use actual stored primary-FLOW selections.

NERSC's published October 14 JSON is corroboration-only. For the 02:55 CDT
acquisition, our direct fit gives 184.227 µs while that reference gave 134.546 µs.
The baseline fitter and pinned-PR comparison also gave 184.227 locally. A 15:46:45
CDT acquisition fails our drift-support check while the reference lists 88.593 µs.
The reference lacks file hashes/selector versions. Production generation,
settings or runtime differences remain unresolved; do not force agreement by
discarding source files or substituting reference values.

## 10. Validation and next work

The October 19 local suite ran all 90 tests: 89 passed in the sandbox, and the
HTTP-range test initially failed because the sandbox forbids its local listening
socket. That specific test then passed with the required socket access. No
scientific assertion was weakened to obtain a pass. Shell syntax checks and
independent archive/mean/plot checks are recorded in the final refresh receipts.

The fresh replay uses the complete retained selection/object inventory, including
the 235 weekly primary-FLOW captures and the six October 14 primary/native pairs.
Historical compact arrays are refitted, not substituted published lifetime
numbers. This is not a survey of all NERSC files or new detector coverage after
October 14. The external drive was available during this replay. See the local
refresh README and validation records for exact successful/failing fit counts.

Before declaring NERSC operational:

- Check the installed commit, Python numerical environment and upstream FLOW
  selection/configuration; record them with the first validation run.
- Verify read-only inputs, matched native files, calibration accessibility,
  real queue/worker execution, late-packet revision handling and snapshot-only
  publication in a preview tree.
- Verify actual web permissions/symlink handling and rollback before promotion.
- Establish DAQ05 export and transfer scheduling without embedding credentials.
- Investigate the stale PRM source if newer operator measurements are expected.
- Retain the known packet low-lifetime bound, gas calibration/plateau limitations,
  conditional uncertainty definitions and unresolved NERSC comparison in any
  physics interpretation. Treat improvements as explicit reviewed changes.

Do not broaden a backfill or bulk-download detector data by default: previous
weekly inventories approached 824 GB. Prefer existing NERSC files or verified
compact captures. The operator explicitly preferred all downloaded files, with
any extra streaming bounded to about an hour. Update the user on meaningful
progress and report unavailable evidence directly.
