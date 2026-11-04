# Zero-drift charge-scale diagnostics

`actions/lifetime/charge_scale.py` reports two **different observables**. Neither
changes the existing lifetime result or its selection.

* **Saved lifetime-proxy intercept:** reads the accepted exponential amplitude
  and covariance already recorded in `summary.json`. The legacy estimator fits
  the MPV of a discrete convolution of the hit-density and charge-density
  histograms. Its amplitude has **legacy proxy units**, not raw charge/cm. This
  monitor needs only the summary files; there is no refit.
* **Direct dQ/ds intercept:** optionally fits Langau MPVs of `dQ / dx` from the
  **entire pooled segment dataset**, then fits those MPVs to
  `A exp(-t / tau)`. Here `dx` denotes segment path length, so dQ/ds avoids
  confusion with the detector X coordinate. The units are the recorded FLOW
  charge units per cm. This is a separate diagnostic estimator, not a
  replacement for the legacy lifetime result or an absolute charge calibration.

For either exponential, `A` is the extrapolated value at drift time zero. The
saved errors are `sqrt(cov[A,A])` and `sqrt(cov[tau,tau])`; the payload also
contains their correlation and a curve with a band propagated from the full
two-parameter covariance. These are conditional fit uncertainties: detector
calibration, selection changes, fit-model choice and correlations among segments
from the same track are not included. The displayed fit range begins above zero,
so the segment from zero to the first fitted bin is an extrapolation.

The direct diagnostic uses the existing prospective time/charge grids: 20 time
edges between 0 and 196 us, 70 charge edges between 0 and 250, and nominal time
bins 4 through 17 inclusive. It requires at least five in-range segments and five
populated histogram centers in each fitted peak window, a resolved four-parameter
Langau covariance, and at least ten usable time bins. Peak fitting uses
unweighted histogram least squares with the window 0.5–1.5 times the coordinate
of the histogram maximum. The final exponential fit uses MPV covariance errors
as absolute uncertainties. Time bins are half-open, except for the final closed
edge, so boundary segments are never counted twice.

Nonfinite, boundary, non-positive-definite covariance and oscillatory-width
solutions are rejected. The direct diagnostic uses the same 20,000 us lifetime
acceptance ceiling as the lifetime workflow. A failed saved lifetime fit never
produces an accepted proxy intercept. A failed direct fit never produces an
accepted direct intercept; individual usable MPVs and failure reasons remain
available. Large samples do not guarantee an identifiable peak: the fixed charge
grid can leave too few populated centers, especially for small late-drift charge.
Do not treat a missing value as zero or tune bins until an individual result
passes.

## Monitor six-hour windows

Run this in the repository with the same Python environment used for lifetime
fits. Keep the wildcard quoted so Python can expand it recursively:

```bash
python3 actions/lifetime/charge_scale.py \
  '/pscratch/sd/j/jbarrow/lifetime_analysis/six_hour_*/**/summary.json' \
  --output /pscratch/sd/j/jbarrow/lifetime_analysis/charge_scale_windows.html
```

Add `--direct` to also calculate direct dQ/ds intercepts. That option requires the
original full pooled HDF5 referenced by each receipt, its completion marker, the
recorded segment count and dataset SHA-256, and the recorded drift tick period.
Changing or moving a receipt without updating its relative pool path can make
that diagnostic unavailable. Missing pools do not discard an otherwise valid
saved proxy intercept. There is no fallback to fitting the viewer's small track
sample.

The output is self-contained browser HTML plus a neighboring JSON file (override
with `--json-output`). The two observables use separate axes. Series are grouped
by selection mode and recorded run/configuration/estimator identity, and placed
at the UTC midpoint of each acquisition window. Equal recorded configuration
does not establish equal live detector conditions. Failed fits are retained as
null values with reasons, and copied identical receipts are deduplicated.
The plot's PNG export includes the earliest window start and latest window end
in its filename, for example
`charge_scale_20261102T180000Z_to_20261104T060000Z.png`.
Exit status 2 means an input error occurred or no accepted intercept was
available; the report is still written when individual receipts fail.

## Code interface

```python
from charge_scale import charge_scale_from_summary, fit_direct_charge_scale

proxy = charge_scale_from_summary(summary)
direct = fit_direct_charge_scale(full_pooled_segments, tick_us=0.1)
```

Both return `status`, `intercept`, `intercept_error`, `tau_us`, `tau_error_us`,
`covariance`, `correlation`, `fit_message`, `fit_range_us`, `fitpoints` and
`curve`. `fitpoints` contains `x`, `y`, `error`, `bins`; `curve` contains `x`, `y`,
`error`. Time is in microseconds. Accepted scalar results and curves are absent
when `status != 'ok'`. `charge_scale_payload(summary, segments=None, tick_us=0.1)`
combines these under `proxy` and `direct`.

## Related existing work

The nearline lifetime code already fits a pure exponential with an amplitude.
In the neighboring `ndlar_flow` repository,
`src/module0_flow/reco/combined/electron_lifetime.py` also stores `dqdx0` and
`dqdx0_err`, providing a precedent for monitoring charge scale. Its anode
truncation response modifies the exponential, so its stored amplitude must not
be equated with the literal value of that modified curve at zero. Local history
attributes those implementations to DeMario Ross and Peter Madigan,
respectively; no implementation attributed to Roberto Mandujano was found in
the inspected local repositories. The requested monitor follows Roberto's
stated idea without asserting that these earlier implementations are his work.
