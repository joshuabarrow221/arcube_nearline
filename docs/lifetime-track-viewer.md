# Inspect the tracks behind a pooled lifetime fit

`actions/lifetime/track_viewer.py` exports an interactive, self-contained HTML
display from a `recent_tracks.py` receipt. It shows **parent tracks**, rather
than treating each short charge segment as an independent track. The default
export contains up to 100 parent tracks, with 15 visible per browser page;
choose 10, 15, or 20 in the display. Plotly is embedded, so the exported HTML
opens locally without a server, SSH tunnel, CDN, or internet connection.

## Export on NERSC

Activate the same Python 3.10+ analysis environment used for the fit. All
dependencies are already in `actions/lifetime/requirements.txt`.

```bash
cd /global/cfs/cdirs/dune/users/jbarrow/arcube_nearline
python3 actions/lifetime/track_viewer.py /path/to/recent-XXXX/summary.json \
  --output /path/to/track_views/lifetime_tracks__20261102T120000-0600_to_20261102T180000-0600__quality-stored.html \
  --json-output /path/to/track_views/lifetime_tracks__20261102T120000-0600_to_20261102T180000-0600__quality-stored.viewer.json \
  --max-tracks 100 --seed 7
```

Copy the HTML to your laptop and open it in a browser. The optional viewer
JSON can be loaded with **Open another viewer export** in any exported display.
It is a portable coordinate/lineage export, not the original summary receipt.
The HTML already includes its own data; its JSON companion is not required.

**A summary alone is insufficient.** The exporter reads the original FLOW,
saved selection sidecars, and pooled `segments.h5` referenced by the receipt.
Generate the HTML while those files are available on NERSC. A summary copied
into the batch's `summaries/` directory works because its file pointers remain
recorded. For relocated inputs, repeat `--path-map OLD_PREFIX=NEW_PREFIX` as
needed; matching respects path-component boundaries and uses the longest
matching prefix. Source files are opened read-only.

## What the browser shows

- Rotate and zoom the 3D view, or switch to XY, XZ, or YZ projections.
- Color by parent track, segment decision, or segment drift time.
- Toggle original hit points, independent segment pieces, and module bounds.
- Select a parent track or click a segment/hit to inspect its event, source
  file, original row IDs, charge, length, drift time, and electronics address.
- Inspect the selected parent's segment charge-versus-drift distribution.
- Enable **Include rejected context** to show segments removed from those
  same parent tracks. Context points inherit a whole-segment decision; they
  are not individually identified as noise.
- Read the saved **event-level quality evidence**. Those counters include
  other hits in the event and do not assign specific veto reasons to segments.

Coordinates are FLOW x/y/z in **cm**. Drift times are converted using the
recorded source `crs_ticks` to **microseconds**. Hit PPS values are phases of
independent IO-group clocks, not absolute event times. Individual segment
pieces are never connected across gaps or modules. The default view fits the
detector bounds; **Fit all displayed** includes distant reconstructed outliers.
Out-of-frame hit counts are shown when applicable.
When extreme outliers require unequal axis scales in the all-data view, a
visible warning states that angles and length ratios are distorted. Return to
**Fit detector** for equal spatial scales. The plot download buttons include
both recorded time boundaries, selection mode, page and view in PNG filenames.

This is an overlay of different events and files. Apparent crossings do not
imply a single physical interaction. `beam` and `cosmic` retain the upstream
external-trigger timing labels; they are not new particle-origin decisions.

## Sampling and comparisons

The default `--population retained` samples uniformly from original parent
tracks having at least one numerically valid segment in the saved pool. It
does not preferentially choose long or visually straight tracks. The default
`--seed 0` makes sampling repeatable for unchanged inputs. Browser pagination
only navigates that exported sample; increase `--max-tracks` (maximum 1000) or
change `--seed` to inspect additional parents.

Use `--population rejected` for parents with no retained segment and at least
one quality-vetoed segment, or `--population all` for all original selected
parents. These options do not recover tracks rejected by the upstream FLOW
track selector. In `stored` mode no new quality veto is inferred.

For paired parent comparisons, export both `stored` and `quality-stored` with
`--population all` and the same seed, using receipts with identical ordered
source files and original segment layouts. This samples the same parent
population in both exports. Default retained-only samples generally differ
between the two selections. Within one quality export, rejected context always
belongs to the displayed original parent tracks.

At most 2000 referenced hits per parent are displayed by default. Sampling is
uniform and occurs before reading hit payloads, and the display reports the
exported and total referenced hit counts. All segment pieces remain present.
Use `--max-hits-per-track 10000` for more detail, or `0` for segment geometry
only. Nonfinite hit positions are omitted and counted separately in the JSON.
No whole-event hit clouds are loaded or displayed, so unrelated noise elsewhere
in an event is represented only by its saved diagnostic counters.

## Provenance and limits

The exporter verifies each selected dataset's recorded count and SHA-256,
maps sidecar rows back to the original segment rows, verifies byte-identical
segment values, and checks the pooled payload and its member/dataset/row
lineage. Parent identities include file membership, dataset, and row; repeated
file-local IDs are not merged. Original track/event and displayed-hit
associations are checked. Receipts do not contain hashes of the complete raw
hit payloads, so this is not a full-file integrity audit.

`Retained` means included in the pooled **input sample**. It does not imply
that a segment entered a particular fitted histogram bin, or that the physical
track is clean. A successful attenuation fit is not a validation of the cuts;
this display is an inspection aid, not an efficiency measurement.

A failed attenuation fit (including an invalid covariance) can still be
displayed because its pooled inputs exist. Receipts with input errors or no
pool are rejected rather than presented as a complete fit sample. The current
viewer supports `stored` and `quality-stored`; `reselect` uses different lineage
and is explicitly unsupported.
