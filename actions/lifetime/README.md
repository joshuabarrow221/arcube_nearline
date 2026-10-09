# Electron-lifetime monitoring

The lifetime action combines track, packet and slow-control measurements into
one recent-history monitoring window. `watcher.py` schedules updates as completed
inputs arrive; shifters normally use the plot and its linked diagnostics page.

## Reading the plot

Open `elifetime_time_series.png.html` to zoom, inspect points and toggle curves.
The PNG is 3000 × 2000 pixels. `display_days` controls how much recent history is
shown; zooming changes the view, not the averaging windows or fitted values.
The display follows the latest available measurements, so check their timestamps
when assessing freshness, not just when the page was generated.

The deployment example displays UTC for comparison across detector systems;
averaging windows remain Chicago calendar days/six-hour intervals. PRM source
wall times are localized to Chicago before UTC normalization. The operator guide
describes correcting older snapshots that mislabeled those wall times as UTC.
The detailed study PNGs and interactive companions also display UTC. Gas
equivalents are placed at the **center of their six-hour window**, in every view;
for example, 11:00–17:00 UTC is plotted at 14:00 UTC.

| Series | Input | Summary |
|---|---|---|
| **All tracks** | Selected through-going MIP segments in FLOW | Six-hour pooled fit |
| **Cosmic-enriched tracks** | The selected zero-external-trigger subset | Six-hour pooled fit |
| **Packets** | Detector-wide native packet charge objects | Six-hour pooled fit |
| **PRM** | Recorded purity-monitor lifetimes | Daily arithmetic mean |
| **Gas equivalents** | O₂ and H₂O concentrations | Six-hour means, then conversion |

Filled markers represent accepted numerical results. Open packet markers show
finite candidates that failed monitoring checks; consult diagnostics before
interpreting them. Missing points mean no usable estimate, not zero lifetime.
Colored lines connect available points across gaps without filling missing data.
Current windows can change as more inputs arrive.

Error bars show conditional track/packet fit uncertainties or PRM standard errors
of the mean (SEM). Gas uncertainties are not estimated. These bars exclude shared
calibration and selection systematics; an absent bar does not mean zero error.

## Track selections and trigger categories

**All tracks includes cosmic-enriched tracks.** Both use the same upstream
through-going candidate selection. All tracks adds no external-trigger cut;
cosmic-enriched requires the parent event to have `n_ext_trigs == 0`. The
complementary `n_ext_trigs > 0` subset has the compatibility key `beam` and is
retained in diagnostics. These are event-timing labels, not confirmed particle
origins or IFBeam beam-quality selections.

The inclusive fit uses the union of selected segments, not an average of subset
lifetimes. The All tracks and cosmic-enriched estimates share data and are
correlated. “All” refers to available selected candidates, not every reconstructed
track; missing upstream selections or timing associations cannot be inferred.

## Running the monitor

Use Python 3.10+ and [the analysis requirements](requirements.txt). Queued operation
also needs the existing nearline FireWorks environment. Configure the input,
state and publication directories using
[`config/lifetime.example.json`](../../config/lifetime.example.json) and the
[operator guide](../../docs/lifetime-integration.md).

From the repository root:

```bash
# Inspect discovery without connecting to MongoDB.
python3 watcher.py actions/lifetime.sh --lifetime-config /path/to/monitor.json --dry-run --once

# Queue incoming FLOW, native packets and slow-control snapshots.
python3 watcher.py actions/lifetime.sh --lifetime-config /path/to/monitor.json

# Refresh existing state/snapshots without waiting for a new detector file.
bash actions/lifetime.sh --monitor-config /path/to/monitor.json
```

Run workers through the normal `rlaunch` workflow. FLOW and native packet inputs
must have matching relative acquisition paths. Late packets update the same
acquisition without adding duplicate weight. Updates also work with only PRM
and gas data available.

**DAQ05 transfer is a separate step.** The watcher consumes local snapshots; it
does not query DAQ05. [export_purity.sh](../../extra/data_mgmt/export_purity.sh)
queries the configured databases on a reachable host;
[sync_purity.sh](../../extra/data_mgmt/sync_purity.sh) receives those snapshots with
`rsync` and refreshes the monitor. Operators configure access and schedule both
steps as described in the [operator guide](../../docs/lifetime-integration.md#slow-control-export-and-transfer).

## Processing and diagnostics

The action reads source files without modifying them, stores reusable selected
segments/charge objects, and refits compatible six-hour pools when membership
changes. It never averages failed per-file fits or imports measurements from a
published plot. PRM readings are deduplicated before averaging. Gas quality checks
reject sustained plateaus, inadequate coverage and configured calibration limits.

Publication preserves the existing filenames:

- `elifetime_time_series.png` / `.png.html`: the monitoring window;
- `elifetime_diagnostics.html`: fit failures, coverage and input-processing errors;
- `elifetime_status.json`: values, uncertainties and provenance.

If a series stops updating, check diagnostics, input timestamps and worker logs.
A failed fit remains a status record; other available methods can still update.
The last complete web generation survives a publication failure. See the
[operator guide](../../docs/lifetime-integration.md#recovery-and-maintenance) for
retries, permissions and rollback.

## Further reference

- [Methods and implementation](../../docs/lifetime-methods.md): pooling, errors,
  quality gates, code map and standalone analysis tools.
- [Packet analysis](README_packet_lifetime.md): extraction inputs and calibration.
- `python actions/lifetime/lifetime.py --help`: all supported analysis/export options.

The configured watcher profile owns its persistent pooled state. Older per-file
and raw-points study commands remain available separately; their histories are
not interchangeable with the pooled monitoring history.
