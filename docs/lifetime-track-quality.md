# Hit-level quality checks for lifetime segments

The rolling track CLI can veto segments using their original FLOW hits before
pooling them. These cuts address repeated hot pixels and spatially broad,
coherent electronics firing. **The default thresholds are provisional.** The
synthetic tests verify their implementation; they do not establish signal
efficiency, noise rejection on detector data, or an unbiased lifetime.

See [the rolling CLI guide](lifetime-recent-tracks.md) for commands, acquisition
time windows, output provenance, and fit results.

## Selection modes

`recent_tracks.py --selection quality-stored` is the default. It follows the
stored segment → track → event and segment → prompt-hit references. It reads
the original event's prompt hits, including hits outside the selected track,
and rejects a whole segment if any of its member hits fails quality. Surviving
segment rows retain their original `dQ`, `dx`, `nhits`, and drift time. Missing
or inconsistent associations are errors, rather than evidence that a segment
is clean. This mode cannot recover tracks that the original selection omitted.

`--selection reselect` reruns the installed `ndlar_flow` through-going track
selection and writes a separate sidecar. It evaluates electronics quality
before PCA cleaning, then vetoes affected whole candidate segments. The local
adapter partitions each io_group along its spatial principal axis in nominal
2-cm intervals, using the actual length of a short final interval. Each hit
belongs to at most one interval. This fixes the legacy segment builder's
constant-z failure and inclusive-boundary duplication. This is a different
segment definition, so its outputs must not be mixed with legacy segmentation.

`--selection stored` pools existing segment arrays without these new hit-level
checks. It supports the colleague's gather-and-fit workflow, including a
custom `--dataset`. Its results explicitly retain unverified upstream noise
quality. Segment-only arrays contain too little information to establish
pixel rates or chip/tile/module coincidences.

Source FLOW files are read-only. Selected segments, event diagnostics, quality
settings, and code identities are written to separate outputs. **Rejected
hits are not subtracted from a retained segment's charge.** Such subtraction
could reduce charge and hit density in a way correlated with drift time.
The inherited track candidate's PCA cleaning remains part of reconstruction;
the new electronics veto does not introduce another charge-removal step.

The existing single-file entry point also accepts `--quality-select-segments`:

```bash
python actions/lifetime/lifetime.py --input_file /data/run/packet-2026_11_03_06_00_00_CST.FLOW.hdf5 \
  --quality-select-segments --output_file_plot /data/study/fit.png \
  --output_file_json /data/study/history.json
```

The configured watcher continues to consume its existing stored-segment
workflow. These new cuts are active through the rolling CLI's default mode or
the explicit single-file option, not retrospectively applied to published
monitor values.

## Clock, geometry, and electronics identities

The coincidence clock is `charge/calib_prompt_hits/data.ts_pps`, not
`t_drift` or segment `t`. FLOW's `run_info/crs_ticks` supplies microseconds per
charge tick. At the standard 0.1 µs per tick, **50 µs is 500 ticks**. The lifetime
fitter's historical `0..1960` drift edges are ticks, corresponding to
`0..196` µs; they must not be interpreted as 1960 µs.

The FLOW adapters use the stored tile lookup and pair adjacent io_groups into
the four 2x2 modules. A physical pixel is identified by
`(io_group, tile_id, chip_id, channel_id)`, so alternate UART paths to one pixel
do not inflate distinct-channel occupancy. The standalone mask function also
supports explicit tile/module identities. Its documented fallback maps four
consecutive io_channels to one tile and two io_groups to one module; that
fallback is specific to 2x2 and is not a generic detector geometry.

Unmapped hits and hits outside their own module's xyz readout bounds, allowing
2.1 cm of reconstruction tolerance, cannot enter retained segments. Finite
out-of-volume hits remain available when evaluating electronics activity, so
trimming a long noise streak to the detector volume cannot hide its repeated
firing. Nonfinite required coordinates or arrival times fail hit quality.

When FLOW records `rollover_ticks`, timestamps are unwrapped **separately for
each io_group**, cutting at the largest circular gap. Each group's interval
must fit within half a rollover period. A group that cannot satisfy this
condition has ambiguous arrival timing: its hits and associated segments are
vetoed, with `ambiguous_clock_hits`, `ambiguous_clock_groups`, and `clock_domains`
recorded in the event diagnostics. Other evaluable groups and events continue;
one ambiguous event does not invalidate the entire file. This timing veto is
recorded separately from detected hot-pixel or coherent noise.

The io_group is the clock domain used by FLOW's timestamp correction and SYNC
unrolling. Combining independent groups into one short event interval can fail
even when each group is compact. Sustained activity can also produce genuinely
broad events within a group. Compact PPS phases alone do not establish the
absolute epoch of hits in a multi-second event; recovering ambiguous groups
would require packet/SYNC information beyond the prompt-hit phase field.

Without rollover metadata, no rollover is inferred. Source clock and drift
axis metadata take precedence over conflicting CLI overrides.

## Default cuts

All thresholds below can be set through a JSON object passed with
`--quality-config`. The complete effective configuration is saved with the
selection. A changed configuration or quality implementation separates pools.

For each physical pixel, the hot-pixel test scans sliding windows of at most
`window_us = 50`. It requires all of the following:

| Setting | Default | Requirement |
| --- | ---: | --- |
| `hot_min_hits` | 12 | At least this many hits from one physical pixel |
| `hot_min_rate_hz` | 100,000 Hz | `(N−1)/(latest−earliest)` reaches this rate; time denominator is at least one charge tick |
| `hot_min_drift_span_cm` | 2 cm | Sufficient extent along the configured drift axis |
| `hot_max_transverse_span_cm` | 0.5 cm | Upper bound on the norm of the other two coordinate extents |
| `hot_min_axis_cosine` | 0.995 | Drift extent divided by the norm of all three coordinate extents |

The direction test uses coordinate extents, not a blanket track-angle cut.
Separate physical pixels are not combined to manufacture the repetition count.
Each tested window starts on a hit and extends through the last hit within
50 µs; rate is evaluated over that window's occupied span. All hits belonging
to a qualifying window are marked. A real track confined to a rapidly firing
single pixel can satisfy the same conditions: **this topology is ambiguous and
the cut deliberately sacrifices it**. The rate is an event-local rate; no
run-wide hot-pixel history is inferred.

Coherent firing uses the same sliding arrival-time windows. It requires hit
multiplicity, distinct physical channels, and (for larger electronics regions)
activity across multiple smaller regions:

| Region | Minimum hits | Minimum distinct channels | Additional breadth |
| --- | ---: | ---: | --- |
| Chip | `chip_min_hits = 32` | `chip_min_channels = 24` | One physical chip |
| Tile | `tile_min_hits = 128` | `tile_min_channels = 96` | `tile_min_chips = 4` |
| Module, within one io_group | `module_min_hits = 512` | `module_min_channels = 384` | `module_min_tiles = 4` |

All coincidence tests stay within one io_group. The module test combines its
tiles on that readout side; it does not sum subthreshold activity across the
two independently timed io_groups. Independently shifting both clocks to zero
and then combining them would manufacture simultaneity.

Every coherent-firing test also requires
`burst_min_transverse_rms_cm = 0.5`: the minor-axis RMS of the distinct pixels'
coordinates in the plane perpendicular to drift must reach this value. Each
physical pixel contributes once to this spatial test, regardless of its hit
count. This protects dense, approximately straight projected tracks from an
occupancy-only veto. Conversely, noise arranged in a thin line, or firing below
the configured multiplicities, can escape this burst test.

Windows slide across arbitrary 50-µs boundaries, and the mask is the union of
all qualifying windows. Repeated or continuous coherent firing can therefore
remain vetoed across several windows. Single isolated hits do not trigger a
hot-pixel or coherent-firing veto. There is no added low-charge threshold.

## Module discontinuities and interpretation

The new quality veto does not require tracks in different modules to join
spatially. Fresh segmentation is local to an io_group, and stored segments are
checked through their own hit memberships. The inherited upstream candidate
selection still requires a through-going, sufficiently straight track. Its
global PCA-distance and face-crossing requirements may already have discarded
candidates with large module offsets; a later segment veto cannot restore them.

A spatial discontinuity alone need not spoil a local charge-density
measurement. A wrong module timing offset can still bias the segment's drift
time, which is the independent variable of the attenuation fit. These checks
do not determine or correct module timing offsets. The permission to retain
locally clean segments is not evidence that their absolute drift times agree.

The checks operate on one FLOW event at a time. Noise divided across event
boundaries is not combined, and long-term pixel rates are not measured. Their
efficiency and effects on drift acceptance should be checked on labeled noisy
and clean samples before treating the fitted uncertainty as a precision purity
measurement. A numerically successful pooled fit does not perform that check.

## Implementation

- [track_quality.py](../actions/lifetime/track_quality.py): configurable mask and
  event-local electronics diagnostics.
- [stored_quality.py](../actions/lifetime/stored_quality.py): authoritative FLOW
  references, clock/geometry context, and stored-segment veto.
- [track_selection.py](../actions/lifetime/track_selection.py): fresh candidate
  selection, local projection segmentation, and whole-segment veto.
- [test_track_quality.py](../tests/test_track_quality.py): synthetic hot-pixel,
  aliasing, coincidence-clock, sliding-window, sustained-noise, geometry, and
  straight-track checks.
- [test_track_selection.py](../tests/test_track_selection.py): orientation,
  module separation, charge conservation, and unique-membership checks.
