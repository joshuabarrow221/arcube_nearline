# Lifetime monitor: methods and implementation

This reference describes the deployed calculation contract. For routine use see
the [monitoring README](../actions/lifetime/README.md); setup and recovery are in
the [operator guide](lifetime-integration.md).

## Inputs, pooling and provenance

The configured profile uses fixed local calendar windows: 00–06, 06–12, 12–18
and 18–24 in `timezone` (default America/Chicago). PRM uses the calendar day.
DST changes elapsed duration. A file's entire selection belongs to its filename
start-time window; events are not split across boundaries. Pools are not widened
to obtain an acceptable fit. Points use mean member-file start time or mean PRM
reading time; gas points use the window midpoint. These do not imply continuous
exposure. Incomplete windows update as inputs arrive.

Detector observations are pooled **before fitting**. Default cohorts are full
acquisition directories; reconstruction/selection/calibration fingerprints must
also agree. Explicit `cohort_rules` can merge independently verified compatible
runs; ambiguous matches fail. `diagnostic_patterns` keeps named hardware/noise
studies in diagnostics. No file is removed because its fitted lifetime is high
or low, and per-file fit success is not a membership requirement.

Source FLOW/native files remain read-only. Compact array objects are immutable
and hash-checked; each input has one active descriptor, with old objects and fits
retained for audit. Source identity includes the full relative path. Stored
history records membership, configuration, selection provenance and failure
reasons. Plot PNG/JSON from another publication is never a measurement input.

## Tracks

The fitter consumes selected segment `dQ`, `dx`, `nhits` and drift `t`. It rejects
nonfinite values, nonpositive length/charge/hit count, and uses the upstream
through-going candidate selection. Existing FLOW selections are required by the
configured profile; the optional standalone selection adapter can create sidecars.

`all_mip` (All tracks) is the union of two event-timing subsets: `cosmic` has
`n_ext_trigs == 0`, and the legacy `beam` key has `n_ext_trigs > 0`. Both apply the
same candidate cuts. Explicit event→track→segment references partition legacy
segments; stored split datasets must have matching selection settings and both
categories. Recorded group/dataset settings enter the pool fingerprint; counts
and previous lifetime results do not. Explicit timing labels must be `beam` and
`off_beam`, respectively. Legacy files without recorded settings retain an empty
metadata record: equal missing metadata does not verify equal selector versions,
so independently reviewed acquisition cohorts remain necessary. Missing
associations allow only the combined stored selection.
No branch restores candidates excluded upstream, identifies particle origin,
or applies IFBeam/POT cuts.

All tracks includes cosmic-enriched segments, so their fits are correlated.
If all selected events have zero external triggers, their inputs coincide; if
none does, the cosmic-enriched fit is unavailable. The inclusive lifetime is
refitted from the union, not averaged from the two subset lifetimes. Different
subset statistics/drift coverage can change the result or prevent a fit.

DeMario's estimator in [lifetime_funcs.py](../actions/lifetime/lifetime_funcs.py)
uses 19 drift slices over 0–1960 ticks at 0.1 µs/tick. Every slice needs at least
five valid segments and nonempty fixed charge/hit-density histograms. The
convolved proxy peak uses its historical fit window, and attenuation uses the
historical subset of drift slices. Missing early support or a failed peak model
remains a failure even if later slices are populated.

Each slice MPV error comes from the covariance of an unweighted SciPy
`curve_fit`, including its residual scaling. The exponential attenuation fit
uses those MPV errors with `absolute_sigma=True`; the lifetime error is
`sqrt(C_tau,tau)` in µs. These are conditional fit errors, not MINUIT estimates,
file-to-file SEMs or total experimental uncertainties. Shared calibration,
selection, track/event correlations and within-window variation are excluded.
Errors of the overlapping track curves cannot be combined as independent errors.

## Packets

The configured extractor uses the exact static pedestal JSON named in matching
FLOW, verified against linked charge/ADC pairs. It groups signed `ADC - pedestal`
charge by physical channel using the existing 28-tick successive-packet rule,
then fits detector-wide charge objects relative to EXT. No TPC-by-TPC lifetime
average is taken. EXT settings and calibration/extraction hashes separate pools.
The grouping and timing convention remain provisional; see
[packet inputs and tools](../actions/lifetime/README_packet_lifetime.md).

The current prospective monitoring gates are:

- ≥50 objects per 10-µs drift slice; positive Landau–Gaussian MPV and slice
  chi-square/ndf <3;
- ≥80% successful MPV bootstrap trials, with at least 10 successes;
- ≥8 accepted fit slices within 20–180 µs spanning ≥100 µs;
- interior, nonsingular exponential fit with chi-square/ndf ≤3 and positive
  inverse lifetime resolved beyond one covariance sigma;
- ≥95% pedestal coverage in **every** member file, not average coverage.

By default each slice uses 100 object-resampling bootstrap trials with fixed
seeds. Its MPV error is the sample SD of successful replica MPVs (`ddof=1`), not
that SD divided by sqrt(number of trials). The final weighted residual Jacobian
uses covariance `pinv(J.T @ J)` without reduced-chi-square rescaling. With inverse
lifetime alpha in ms⁻¹, `tau_us = 1000/alpha` and
`sigma_tau_us = 1000*sigma_alpha/alpha**2`. There is no further division by file
count. Object resampling assumes independence; shared EXT/event/channel,
calibration and model effects are not included.

Finite positive failed candidates can be displayed open while keeping their
failed status. A conditional error is recovered only when a valid attenuation
covariance survives a later pedestal-coverage veto. Bad-model/unresolved
candidates are not given invented error bars; nonfinite/nonpositive estimates
have no plotted point. Gates are not tuned to agreement with another method.

## PRM and gas

Raw snapshot import normalizes timestamps/units, removes exact duplicates and
rejects conflicting readings. PRM averages distinct accepted **recorded
lifetimes**, without first inverting them. Daily uncertainty is
`SEM = sample_SD / sqrt(N)` for N>1; a singleton has no estimated error. Calibration
systematics and correlations between readings are excluded.

PRM source timestamps are timezone-free Chicago wall-clock readings; localize
them before UTC normalization and calendar-day grouping. This source convention
is operator-specified, rather than inferred from the database session timezone.
Ignition's epoch-based timestamps require no corresponding clock shift.
The detailed PNG and interactive plots display UTC; their `--timezone` option
still defines grouping boundaries, not the display clock.

Gas quality is evaluated across the retained archive before six-hour averaging.
Defaults require ≥3 hours of coverage, no sampling gap >30 minutes, nonnegative
concentrations and no ≥3-hour exactly flat run. Per-source tolerance and calibration
bounds are configurable. Plateau detection crosses window/export boundaries;
gaps reset it. There is no forward fill. A rejected analyzer does not contribute
a lifetime equivalent.

Every gas-equivalent marker and connecting line uses the **midpoint between the
averaging window's start and end**, in both detailed views and the shifter view.
The value represents mean concentrations over that window, then one conversion;
it is not assigned to the left/right bin edge or latest readout. For example,
06:00–12:00 Chicago during daylight time maps to 11:00–17:00 UTC and is plotted at
14:00 UTC. Compute the midpoint from timezone-aware instants, including DST
transitions. Partial windows retain that nominal center and their actual
coverage metadata; neither centering nor connecting lines implies full exposure.

Mean O₂ and H₂O concentrations are converted once using
`tau_us = 1 / (k_O2 * mean_O2_ppb + k_H2O * mean_H2O_ppb)`.
The supplied `eva_20261006` config uses a field-dependent oxygen coefficient and
water coefficient from [plot_lifetime_history.py](../actions/lifetime/plot_lifetime_history.py),
with a 212.125-V/cm **reference** field. Rows retain model/field/coefficient
provenance. This does not establish uniform TPC or gas/liquid conditions; no
unverified offset or phase correction is applied. Without an explicit conversion
config, the legacy coefficients are `1/299` and `1/17000` µs⁻¹ ppb⁻¹.

Each O₂ analyzer remains separate. Missing/rejected water is not treated as zero;
O₂-only equivalents are recorded separately and can be enabled with `gas_samples`.
N₂ is archived without conversion. Gas instrument/conversion errors are not
estimated, so no gas uncertainty bars are supplied.

## Code map and standalone tools

| Module | Responsibility |
|---|---|
| `lifetime.py` | Central CLI, track fitting and legacy study interface |
| `monitor_config.py`, `monitor.py` | Input pairing/revisions, serialized state and publication |
| `track_shards.py`, `track_pooling.py` | Selected-segment storage, provenance and joint fits |
| `packet_lifetime.py`, `packet_pooling.py` | Charge extraction, fixed-window fits and bootstrap checks |
| `packet_pedestal.py`, `calibration_pedestal.py` | Linked FLOW/static calibration adapters |
| `purity_sources.py` | Database export, normalization, PRM means and gas quality/conversion |
| `monitor_plot.py`, `lifetime_io.py` | Plot/diagnostics rendering, locks and atomic outputs |
| `flow_input.py`, `track_selection.py`, `track_audit.py` | Optional remote reads, sidecar selection and per-file contribution review |

Run `python actions/lifetime/lifetime.py --help` for per-file fits, raw-point
plots, reference-only corroboration, review rules, database export and optional
metadata writes. `lifetime_timeseries.py` retains its compatibility entry point.
`packet_pooling.py --help` exposes separate fixed 6/24/48-hour studies.

These expert paths are separate from `--monitor-config`: legacy track summaries
average per-file lifetimes, whereas the monitor fits pooled segments. Legacy
raw-point companions are not regenerated by the configured profile. Do not mix
those histories or interpret their different aggregates as identical products.
