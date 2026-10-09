# Lifetime monitor: operator guide

For plot interpretation see the [lifetime README](../actions/lifetime/README.md).
[Methods and implementation](lifetime-methods.md) covers fitting details.

## Configuration

Install the Python 3.10+ analysis dependencies in the worker environment:

```bash
python3 -m pip install -r actions/lifetime/requirements.txt
```

Use [config/lifetime.example.json](../config/lifetime.example.json) as the deployment
profile. Paths may be absolute, relative to the config file, or environment
variables. Watcher and workers must resolve them identically.

| Setting | Purpose |
|---|---|
| `flow_root` | Completed `.FLOW.hdf5` files with stored track selections |
| `packet_root` | Native `packet-*.h5` files in the matching relative hierarchy |
| `snapshot_root` | Flat archive of immutable `slowcontrols-*.json` files |
| `state_root` | Private shards, receipts, fit caches and history |
| `publish_root` | Served plot, diagnostics and status files |
| `display_days` | Recent interval displayed; does not change pooling |
| `display_timezone` | Axis clock; example uses UTC, independently of Chicago averaging windows |
| `maximum_display_us` | Display cap; omitted values remain in diagnostics |
| `cohort_rules`, `diagnostic_patterns` | Reviewed acquisition grouping and diagnostic-only runs |
| `gas_quality_config`, `gas_conversion_config` | Analyzer validity and conversion policy |

The example uses environment variables for all five roots. Set these on both
watcher and workers, with actual site paths:

```bash
export ARCUBE_LIFETIME_FLOW_ROOT=/data/nearline/flowed_charge
export ARCUBE_LIFETIME_PACKET_ROOT=/data/nearline/packet
export ARCUBE_LIFETIME_SNAPSHOT_ROOT=/data/purity/snapshots
export ARCUBE_LIFETIME_STATE_ROOT=/data/purity/state
export ARCUBE_LIFETIME_PREVIEW_ROOT=/data/purity/preview/plots
export ARCUBE_NEARLINE_LIFETIME_CONFIG="$(realpath config/lifetime.example.json)"
# Optional: select the installed analysis interpreter explicitly.
export ARCUBE_NEARLINE_PYTHON=/path/to/analysis/environment/bin/python
```

State and publication directories must be disjoint and outside watched trees.
Use a separate preview directory before changing the served output. Choose
`display_days` for the desired monitoring view, such as seven recent days; the
example's default is 35. This setting limits display, not retained history.
PRM remains daily and the configured detector/gas windows remain six hours.

Packet extraction defaults to enabled, detector-wide, with EXT I/O group 6 and
100 MPV bootstrap trials. Verify trigger wiring and calibration at the site.
Optional `calibration_root` relocates the exact FLOW-named pedestal basename.
Source files and calibration products must be immutable once published.

## Watcher and workers

Load the normal nearline FireWorks configuration, then run:

```bash
python3 watcher.py actions/lifetime.sh \
  --lifetime-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG" --dry-run --once
python3 watcher.py actions/lifetime.sh \
  --lifetime-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG" --min-file-age 60
# In the configured worker environment:
rlaunch rapidfire --nlaunches infinite
```

Dry-run discovery requires only Python's standard library and does not contact
MongoDB. Queueing requires FireWorks. Run one watcher per configuration. Producers
should close files and rename them into place; minimum age alone cannot prove
completion. Do not use `rsync --inplace` into watched names.

The watcher checks present FLOW/native dependencies, queues a canonical FLOW
input, and revisits it when a late packet or changed revision appears. Packet-first
arrival waits for FLOW. Snapshot arrivals independently trigger publication.
Full relative paths distinguish acquisitions with identical basenames. A worker
serializes state updates and replaces each acquisition's active descriptor, so
retries do not duplicate its contribution. Unchanged pool fits are cached.

For a manual update or initial smoke test:

```bash
bash actions/lifetime.sh --monitor-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG"
bash actions/lifetime.sh /data/nearline/flowed_charge/run/packet-NAME.FLOW.hdf5 \
  --monitor-config "$ARCUBE_NEARLINE_LIFETIME_CONFIG"
```

Use an actual completed acquisition filename with the standard nearline timestamp.
A no-input refresh uses saved shards and snapshots; initial backfill requires
processing original FLOW/native inputs, not importing published lifetime values.

## Slow-control export and transfer

**The watcher does not connect to DAQ05.** Export and transfer are implemented
companion steps; their credentials, SSH route and schedules are site configuration.

| Step | Host | Action |
|---|---|---|
| Export | DAQ05 or another database-reachable host | `export_purity.sh` queries sources read-only and writes completed snapshots |
| Transfer | NERSC/receiving host | `sync_purity.sh` copies snapshots with `rsync` and refreshes the monitor |
| Watch | NERSC/receiving host | `watcher.py` detects received snapshots and queues the same idempotent update |

Configure [sources.example.json](../actions/lifetime/sources.example.json) with the
verified database schema and units. Its current mapping is PSQL PRM lifetimes in
seconds; Ignition H₂O tag 1893 in ppb, O₂ tags 1890 in ppb and 1874 in ppm, and N₂
tag 1871 in ppm. N₂ is archived but not converted to lifetime. A legacy Influx
adapter also exists and requires an explicit measurement, field and units.

The PRM `timestamp` column has no timezone. Interpret it as the acquisition
computer's Chicago wall time (`naive_timezone: America/Chicago`), as specified
by the operator, and normalize to UTC only after localization. Query boundaries
use that same wall clock. Ignition epoch timestamps remain UTC instants. Earlier
exports that labeled PRM wall time as UTC must be corrected or re-exported before
import; do not mix those mislabeled readings with corrected snapshots.

On the export host, supply protected `PURITY_PSQL_URL` and `PURITY_IGNITION_URL`
environment variables and explicit ISO timestamps with offsets:

```bash
bash extra/data_mgmt/export_purity.sh /path/to/sources.json \
  "$START_ISO" "$END_ISO" /path/to/purity_exports
```

On the receiver, using the configured monitor environment:

```bash
export ARCUBE_NEARLINE_PURITY_RSYNC_SOURCE='user@host:/path/to/purity_exports/'
# Set RSYNC_RSH if the approved site route requires an SSH jump host.
bash extra/data_mgmt/sync_purity.sh
```

Schedule exports and transfers independently; running the watcher does not install
these schedules. Overlapping exports, for example a rolling 48-hour query, cover
late readings and plateau lookback. Retain older snapshots; import deduplicates
identical readings and rejects conflicting same-source/time values. Exports
publish only after every enabled query succeeds. Transfers exclude partial files
and do not delete history. The receiver uses the config's `snapshot_root`.

## Recovery and maintenance

Scientific fit failures produce diagnostic records; extraction/IO failures also
make the job fail after other available sources are published. Inspect diagnostics
and worker logs, fix the cause, then explicitly rerun the failed FireWork with
`lpad rerun_fws -i <fw_id>`. The watcher does not repeatedly resubmit a failed
unchanged revision. `--max-file-age`, if used, applies to the newest dependency.
For intentional reprocessing after a change preserving size/mtime, change the
config's `reprocess_token`; never silently replace calibration contents in place.

Use the default shared-filesystem lock on NERSC. Keep the worker wall limit at
most six hours, below the 24-hour lock lease. Set `ARCUBE_NEARLINE_LOCAL_OUTPUT=1`
only for a single-host local/WSL setup; never mix lock backends for one state
namespace. Size worker RAM for complete native timing arrays. State, fit caches
and old plot generations are retained, so review storage use periodically.

Publication renders a complete `.lifetime-generations/` directory and atomically
switches `.lifetime-current`; stable PNG/HTML/status aliases follow it. Existing
regular outputs are backed up in `.previous/` on initial migration. Before
promotion, verify the served preview, symlink handling, permissions and worker
behavior, then stop the old publisher. The first alias migration is sequential;
subsequent generation switches are atomic. Keep the existing output filenames to
preserve shifter URLs. For rollback, stop writes and switch to a retained
generation, or restore the legacy backups and publisher.

For code changes, run `python -m pytest tests -q` in the analysis environment
(with `pytest` installed), and `bash -n` on the modified shell wrappers.
