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
