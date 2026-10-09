#!/usr/bin/env python3
"""One entry point for direct FLOW fits, slow controls and purity publication.

Keep three grains distinct: per-file track fits, individual PRM observations,
and time-window summaries. Published reference values never enter these
calculations. The legacy per-file path fits before taking its publication lock.
The configured six-hour profile instead serializes ingestion, pooled fits and
publication together, so retries cannot change membership under another fit.
"""
import argparse
import json
import os
import re
import sys
import matplotlib
matplotlib.use("Agg")
from pathlib import Path
from lifetime_io import public_permissions

import h5py
import matplotlib.pyplot as plt
import numpy as np
from lifetime_io import output_lock as Lock

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nearline_util import date_from_filename


# These keys are an event-timing partition of the *same selected through-going
# candidate population*, not separate particle-identification algorithms.
# ``beam`` is a legacy key for n_ext_trigs > 0, not proof of beam origin.
# ``cosmic`` means n_ext_trigs == 0 and is therefore a subset of ``all_mip``.
# Neither category adds an IFBeam/POT cut or a cosmic-direction classifier.
SPLIT_SAMPLES = (
    ('beam', 'analysis/beam_rock_muon_segments/data'),
    ('cosmic', 'analysis/cosmic_muon_segments/data'),
)
LEGACY_SAMPLE = ('all_mip', 'analysis/rock_muon_segments/data')

# Keep the short legend labels while exposing their precise scope in hover
# and diagnostics. "All" is limited to available upstream-selected candidates,
# not every reconstructed track, hit, event, or physical muon in the detector.
TRACK_SAMPLE_DESCRIPTIONS = {
    'all_mip': ('All available selected through-going MIP candidates, with no additional '
                'external-trigger cut. Includes the cosmic-enriched subset.'),
    'cosmic': ('The same selected candidates, restricted to events with n_ext_trigs == 0. '
               'A subset of All tracks; cosmic origin and beam-off status are not established.'),
    'beam': ('Selected candidates in events with n_ext_trigs > 0. The legacy beam key '
             'means externally triggered, not an IFBeam-confirmed beam-muon sample.'),
}


def parse_sample_spec(spec):
    """Parse a ``LABEL=HDF5_PATH`` command-line sample definition."""
    try:
        label, dset_path = spec.split('=', 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f'sample must have the form LABEL=HDF5_PATH; got {spec!r}'
        ) from exc

    if not label or not dset_path:
        raise argparse.ArgumentTypeError(
            f'sample must have the form LABEL=HDF5_PATH; got {spec!r}'
        )
    return label, dset_path.lstrip('/')


def find_samples(h5_file, requested_samples=None):
    """Resolve explicit datasets or the standard overlapping display samples.

    The two timing datasets are disjoint. Their union is ``all_mip``; the
    ``cosmic`` curve reuses its no-external-trigger subset. Include each segment
    once in that union, then refit it, rather than averaging category lifetimes.
    A legacy-only dataset means all candidates *stored by its upstream selector*;
    it cannot recover events excluded upstream or establish a cosmic partition.
    """
    if requested_samples:
        missing = [path for _, path in requested_samples if path not in h5_file]
        if missing:
            raise KeyError(f'missing requested segment datasets: {missing}')
        return requested_samples

    split_samples = [(label, path) for label, path in SPLIT_SAMPLES if path in h5_file]
    if split_samples:
        return split_samples + ([('all_mip', [path for _, path in split_samples])] if len(split_samples) == 2 else [])
    if LEGACY_SAMPLE[1] in h5_file:
        return [LEGACY_SAMPLE]
    raise KeyError('no beam/cosmic or legacy rock-muon segment dataset found')


def valid_segments(segments):
    """Remove unusable segment rows before histogramming and fitting."""
    mask = (segments['dx'] > 0) & (segments['dQ'] > 0) & (segments['nhits'] > 0)
    for field in ('dx', 'dQ', 't', 'nhits'):
        mask &= np.isfinite(segments[field])
    return segments[mask]


def sample_plot_path(output_file_plot, sample, multiple_samples):
    path = Path(output_file_plot)
    if not multiple_samples:
        return path
    safe_sample = re.sub(r'[^A-Za-z0-9_.-]+', '_', sample)
    return path.with_name(f'{path.stem}.{safe_sample}{path.suffix}')


def classify_fit_quality(lifetime_us, error_us):
    """Add a diagnostic flag without turning policy thresholds into fit failures."""
    if not np.isfinite(lifetime_us) or not np.isfinite(error_us) or lifetime_us <= 0:
        return 'nonphysical'
    if error_us < 0 or error_us >= lifetime_us:
        return 'poorly_constrained'
    if lifetime_us > 10_000:
        return 'outside_monitor_range'
    return 'ok'


def fit_samples(input_file, output_file_plot, requested_samples=None, timestamp=None):
    """Fit selected FLOW arrays and retain provenance, including failed fits.

    The internal ``beam`` key is a trigger category, not an IFBeam decision.
    The combined sample pools disjoint segments before fitting; it is not an
    average of the two category lifetimes. Cosmic-enriched segments also enter
    the combined fit, so those two estimates and their errors are correlated.
    Both use the same valid-segment filtering and attenuation estimator; the
    differing event membership can change drift coverage, precision and tau.
    """
    from purity_sources import aware_time
    timestamp = aware_time(timestamp) if timestamp else date_from_filename(input_file)
    results = []

    print(f'Opening file: {input_file}')
    from flow_input import open_flow
    import hashlib
    with open_flow(input_file) as (h5_file, provenance):
        samples = find_samples(h5_file, requested_samples)
        multiple_samples = len(samples) > 1

        for sample, dset_path in samples:
            paths = dset_path if isinstance(dset_path, list) else [dset_path]
            raw_segments = [h5_file[path][:] for path in paths]
            dataset_provenance = [dict(path=path, shape=list(array.shape), dtype=str(array.dtype),
                                       sha256=hashlib.sha256(array.tobytes()).hexdigest())
                                  for path, array in zip(paths, raw_segments)]
            segments = valid_segments(np.concatenate(raw_segments))
            print(f'Extracting {sample} lifetime using {len(segments)} segments')

            result = {
                'timestamp': timestamp.isoformat(),
                'sample': sample,
                'method': 'track',
                'calculation_source': 'flow_segments',
                'uncertainty': 'DeMario fit covariance; excludes selection/calibration systematics',
                'input_file': h5_file.attrs.get('source_flow_url', provenance['input_file']),
                'segments_dset': dset_path,
                'n_segments': len(segments),
                'flow_provenance': provenance,
                'flow_datasets': dataset_provenance,
                'fit_status': 'failed',
                'lifetime_us': None,
                'error_us': None,
            }

            if 'source_flow_url' in h5_file.attrs:
                result['flow_segment_cache'] = os.path.abspath(input_file)
                result['flow_source_etag'] = h5_file.attrs['source_etag']
                result['flow_source_size_bytes'] = int(h5_file.attrs['source_size_bytes'])
                result['flow_datasets'] = json.loads(h5_file.attrs['source_datasets'])
            try:
                if not len(segments):
                    raise ValueError('no valid segments')
                import lifetime_funcs as LT
                lifetime, lifetime_error, fig = LT.langau_lifetime(
                    nhits=segments['nhits'] / segments['dx'],
                    dqdx=segments['dQ'] / segments['dx'],
                    time_drifted=segments['t'],
                    time_bins=np.linspace(0, 1960, 20),
                    dqdx_bins=np.linspace(0, 250, 70),
                    nhits_bins=np.linspace(0, 30, 42),
                    wanted_title=f'{timestamp} ({sample})',
                )
                plot_path = sample_plot_path(output_file_plot, sample, multiple_samples)
                plot_path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(plot_path)
                plt.close(fig)

                lifetime = float(lifetime)
                lifetime_error = float(lifetime_error)
                if not np.isfinite(lifetime) or not np.isfinite(lifetime_error) or lifetime <= 0 or lifetime_error < 0:
                    raise ValueError('nonphysical or nonfinite fit result')
                result.update(
                    fit_status='ok',
                    fit_quality=classify_fit_quality(lifetime, lifetime_error),
                    lifetime_us=lifetime,
                    error_us=lifetime_error,
                    plot_file=str(plot_path),
                )
                if hasattr(fig, 'lifetime_diagnostics'):
                    result['fit_diagnostics'] = fig.lifetime_diagnostics
                print(
                    f'Timestamp: {timestamp}, sample: {sample}, '
                    f'lifetime: {lifetime} +- {lifetime_error} us'
                )
            except Exception as exc:
                result['fit_message'] = f'{type(exc).__name__}: {exc}'
                print(f'Lifetime fit failed for {sample}: {result["fit_message"]}')

            results.append(result)

    return results


def write_hdf5_metadata(input_file, results):
    """Attach per-file fit results to each selected-segment HDF5 group."""
    with h5py.File(input_file, 'r+') as h5_file:
        for result in results:
            group_path = str(Path(result['segments_dset']).parent)
            attrs = h5_file[group_path].attrs
            attrs['electron_lifetime_fit_status'] = result['fit_status']
            attrs['electron_lifetime_fit_sample'] = result['sample']
            attrs['electron_lifetime_fit_timestamp'] = result['timestamp']
            attrs['electron_lifetime_fit_n_segments'] = result['n_segments']
            attrs['electron_lifetime_fit_quality'] = result.get('fit_quality', 'fit_failed')
            attrs['electron_lifetime_fit_quality_definition'] = (
                'diagnostic only: positive finite result, error < lifetime, lifetime <= 10000 us'
            )
            attrs['electron_lifetime_fit_method'] = (
                'Landau-Gaussian proxy MPV versus drift time; exponential attenuation fit'
            )
            attrs['electron_lifetime_us'] = (
                np.nan if result['lifetime_us'] is None else result['lifetime_us']
            )
            attrs['electron_lifetime_error_us'] = (
                np.nan if result['error_us'] is None else result['error_us']
            )
            if 'fit_message' in result:
                attrs['electron_lifetime_fit_message'] = result['fit_message']
            elif 'electron_lifetime_fit_message' in attrs:
                del attrs['electron_lifetime_fit_message']


def record_key(entry):
    """Stable replacement identity: retrying a file/window must not add weight."""
    method = entry.get('method', 'track')
    if method in ('prm', 'gas', 'gas_o2', 'packet_pool'):
        return method, entry.get('source'), entry.get('period_start', entry['timestamp'])
    return method, entry.get('input_file', entry['timestamp']), entry.get('sample', 'mixed')


def valid_lifetime(entry):
    """Numerical fit validity, independent of any monitoring exclusion review."""
    value, error = entry.get('lifetime_us'), entry.get('error_us')
    return (entry.get('fit_status', 'ok') == 'ok' and value is not None
            and np.isfinite(value) and value > 0
            and (error is None or (np.isfinite(error) and error >= 0)))


def packet_candidate_points(entries):
    """Expose finite failed-check fits for display without promoting their status.

    The user requested diagnostic packet markers even when monitoring gates
    fail. Keep lifetime_us null and fit_status unchanged in both the original
    history and these copies; display_* fields make the plotting exception
    explicit. Never substitute a value for a window with no finite fit.
    """
    points = []
    for row in entries:
        value = row.get('candidate_lifetime_us')
        if (row.get('method') != 'packet_pool' or valid_lifetime(row)
                or row.get('monitoring_excluded') or value is None
                or not np.isfinite(value) or value <= 0):
            continue
        error = None
        alpha, sigma = row.get('alpha_per_ms'), row.get('alpha_error_per_ms')
        # A retained covariance means the attenuation fit passed its checks
        # before the later pedestal-coverage veto. Recover only that existing
        # conditional statistical error; do not estimate missing systematics
        # or manufacture an error for an unresolved/bad-model candidate.
        if (row.get('covariance') is not None and alpha is not None and sigma is not None
                and np.isfinite(alpha) and np.isfinite(sigma) and alpha > sigma > 0):
            propagated = 1000*sigma/alpha**2
            if np.isfinite(propagated):
                error = float(propagated)
        points.append(dict(row, display_lifetime_us=float(value), display_error_us=error,
                           display_quality='diagnostic candidate; monitoring checks failed'))
    return points


def is_flow_track(entry):
    """Only direct FLOW fits (including documented pre-schema-4 fits) qualify."""
    return (entry.get('method', 'track') == 'track'
            and bool(entry.get('input_file'))
            and (entry.get('calculation_source') == 'flow_segments'
                 or (bool(entry.get('segments_dset')) and 'n_segments' in entry)))


def compare_reference_lifetimes(derived, references):
    """Exact-timestamp corroboration only: references never supply plot values."""
    from purity_sources import aware_time
    comparisons = []
    for reference in references:
        sample = reference.get('sample', 'all_mip')
        if sample == 'mixed':
            sample = 'all_mip'
        candidates = [r for r in derived if is_flow_track(r)
                      and r.get('sample') == sample
                      and aware_time(r['timestamp']) == aware_time(reference['timestamp'])]
        basename = reference.get('published_input_basename')
        if basename:
            candidates = [r for r in candidates if Path(r['input_file']).name == basename]
        comparison = dict(timestamp=reference['timestamp'], sample=sample,
            reference_lifetime_us=reference.get('lifetime_us'), reference_error_us=reference.get('error_us'),
            reference_file=reference.get('reference_history'),
            match_status='unmatched' if not candidates else 'ambiguous' if len(candidates)>1 else 'matched')
        if len(candidates) == 1:
            row = candidates[0]
            comparison.update(input_file=row['input_file'], derived_fit_status=row.get('fit_status'),
                derived_lifetime_us=row.get('lifetime_us'), derived_error_us=row.get('error_us'))
            if valid_lifetime(row) and valid_lifetime(reference):
                delta = row['lifetime_us']-reference['lifetime_us']
                comparison.update(difference_us=delta, relative_difference_percent=100*delta/reference['lifetime_us'])
        comparisons.append(comparison)
    return comparisons


def aggregate_track_lifetimes(entries, timezone_name='America/Chicago'):
    """Equal-weight means of per-file track fits in local six-hour windows.

    Preserve method/selection separation and member provenance. This is a
    summary of fitted lifetimes, not a pooled-segment refit or exposure mean.
    """
    from collections import defaultdict
    import pandas as pd
    from purity_sources import aware_time
    groups, output = defaultdict(list), []
    unique = {record_key(row): row for row in entries}
    for row in unique.values():
        if row.get('plot_role') == 'diagnostic_only':
            continue
        if row.get('method', 'track') != 'track':
            output.append(row)
            continue
        if not is_flow_track(row):
            continue
        stamp = aware_time(row['timestamp']).tz_convert(timezone_name)
        start = (stamp.tz_localize(None).normalize() + pd.Timedelta(hours=6*(stamp.hour//6))).tz_localize(timezone_name)
        groups[(row.get('sample', 'mixed'), start)].append(row)
    for (sample, start), members in sorted(groups.items()):
        end = (start.tz_localize(None)+pd.Timedelta(hours=6)).tz_localize(timezone_name)
        # A deliberate noise study can have a perfectly finite fit. Keep its
        # scientific result in history while excluding it from routine means
        # only through an explicit, reasoned review decision.
        accepted = [r for r in members if valid_lifetime(r) and not r.get('monitoring_excluded', False)]
        values = np.array([r['lifetime_us'] for r in accepted], dtype=float)
        errors = [r.get('error_us') for r in accepted]
        n = len(accepted)
        # Equal file weights imply sqrt(sum(sigma_i^2))/N, not inverse-variance
        # weighting. The larger of this and the observed SEM is a monitoring
        # convention, not a calibrated confidence interval or systematic error.
        propagated = float(np.linalg.norm(errors)/n) if n and all(e is not None for e in errors) else None
        scatter_sem = float(values.std(ddof=1)/np.sqrt(n)) if n > 1 else None
        available_errors = [e for e in (propagated, scatter_sem) if e is not None]
        output.append(dict(timestamp=(start+(end-start)/2).isoformat(),
            period_start=start.isoformat(), period_end=end.isoformat(),
            sample=sample, method='track', aggregation='arithmetic mean of per-file lifetimes, six-hour local windows',
            lifetime_us=float(values.mean()) if n else None,
            error_us=max(available_errors) if available_errors else None,
            propagated_fit_error_us=propagated, between_file_sem_us=scatter_sem,
            uncertainty='larger of propagated independent fit errors and between-file SEM, where available; excludes shared systematics',
            fit_status='ok' if n else 'unavailable', n_measurements=n,
            n_rejected=len(members)-n,
            n_excluded=sum(bool(r.get('monitoring_excluded')) for r in members),
            first_observed_at=min((r['timestamp'] for r in accepted), key=aware_time) if n else None,
            last_observed_at=max((r['timestamp'] for r in accepted), key=aware_time) if n else None,
            members=[dict(timestamp=r['timestamp'], input_file=r.get('input_file') or r.get('published_input_basename'),
                          published_plot_url=r.get('published_plot_url'), run_context=r.get('run_context'),
                          fit_status=r.get('fit_status', 'ok'), lifetime_us=r.get('lifetime_us'),
                          monitoring_excluded=r.get('monitoring_excluded', False),
                          monitoring_review_reason=r.get('monitoring_review_reason'),
                          error_us=r.get('error_us')) for r in members]))
    return output


def connected_values(rows):
    """Connect successive valid observations, including across missing windows.

    Lines are visual guides only; rejected observations remain absent.
    """
    from purity_sources import aware_time
    valid = sorted((r for r in rows if valid_lifetime(r)), key=lambda r: aware_time(r['timestamp']))
    return [r['timestamp'] for r in valid], [r['lifetime_us'] for r in valid]


def apply_track_review(entries, review):
    """Apply an explicit full-source-identity review without deleting any fits.

    A review is a complete replacement policy, not an accumulating blacklist.
    Empty rules restore all numerically valid fits. Exact file identities avoid
    accidentally rejecting unrelated files with similar basenames. No lifetime
    thresholds or automatic outlier clipping are used.
    """
    rules = review.get('rules', [])
    for rule in rules:
        if (not rule.get('input_file') or rule.get('action') not in ('include', 'exclude')
                or not isinstance(rule.get('reason'), str) or not rule['reason'].strip()):
            raise ValueError('track review requires input_file, include/exclude action, and reason')
    output = []
    for entry in entries:
        row = dict(entry)
        if is_flow_track(row):
            row.pop('monitoring_excluded', None)
            row.pop('monitoring_review_reason', None)
            matches = [r for r in rules if r['input_file'] == row['input_file']
                       and (not r.get('sample') or r['sample'] == row.get('sample'))]
            if len(matches) > 1:
                raise ValueError(f"overlapping track-review rules for {row['input_file']}")
            if matches:
                row.update(monitoring_excluded=matches[0]['action'] == 'exclude',
                           monitoring_review_reason=matches[0]['reason'])
        output.append(row)
    return output


def raw_plot_path(output_file):
    """Keep the existing shifter URL and publish its raw companion alongside it."""
    path = Path(output_file)
    return path.with_name(path.stem + '_raw' + path.suffix)


def resolve_gas_conversion_config(config_path, history_path):
    """Keep later source refreshes on the explicitly selected gas model."""
    from purity_sources import gas_conversion_rates
    if config_path:
        config = json.loads(Path(config_path).read_text())
    elif Path(history_path).exists():
        config = json.loads(Path(history_path).read_text()).get('gas_conversion_config')
    else:
        config = None
    gas_conversion_rates(config)
    return config


def update_json(output_file_json, results, output_timeseries=None, timezone_name='America/Chicago', annotations=None, reference_histories=None,
                raw_prm_measurements=None, track_review=None, packet_window_hours=None, packet_snapshot=None, gas_conversion_config=None):
    """Serialize read/update/render so concurrent workers cannot publish stale history."""
    from datetime import timedelta
    from purity_sources import atomic_json
    path = Path(output_file_json)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Lock(str(path) + '.lock', default_timeout=timedelta(seconds=60), lifetime=timedelta(minutes=10)):
        data = json.loads(path.read_text()) if path.exists() else {'lifetimes': []}
        # Migrate unprovenanced/externally published track values out of the calculation history.
        derived, references = [], list(data.get('reference_lifetimes', []))
        for row in [*data['lifetimes'], *results]:
            if row.get('method', 'track') == 'track' and not is_flow_track(row):
                references.append(row)
            else:
                derived.append(row)
        for reference_path in reference_histories or []:
            payload = json.loads(Path(reference_path).read_text())
            references.extend(dict(row, reference_history=str(Path(reference_path).resolve()))
                              for row in payload['lifetimes'])
        indexed = {}
        for row in derived:
            key = record_key(row)
            previous = indexed.get(key)
            # A slower worker must not overwrite a newer evaluation of a
            # growing packet window after waiting for the publication lock.
            if (previous and row.get('method') == 'packet_pool'
                    and previous.get('pool_snapshot_at',previous.get('evaluated_at',''))
                    > row.get('pool_snapshot_at',row.get('evaluated_at',''))):
                continue
            indexed[key] = row
        reference_index = {}
        for row in references:
            key = (row['timestamp'], row.get('sample', 'all_mip'), row.get('published_input_basename'))
            # Retain provenance when the plain published JSON is also supplied.
            match = next((k for k in reference_index if k[:2] == key[:2]
                          and (not k[2] or not key[2] or k[2] == key[2])), key)
            reference_index[match] = dict(reference_index.get(match, {}), **row)
        # Persist the policy so a later watcher refresh/refit cannot silently
        # reinstate a reviewed noise-study file. An explicit empty policy clears
        # it. Keep original fit_status/error/value untouched for reproducibility.
        if track_review is not None:
            data['track_review'] = track_review if isinstance(track_review, dict) else json.loads(Path(track_review).read_text())
        derived = apply_track_review(list(indexed.values()), data.get('track_review', {}))
        # Older histories have only PRM means. Never fabricate raw points from
        # them; importing the original snapshots backfills this separate table.
        from purity_sources import aware_time
        prm_index = {}
        for row in [*data.get('raw_prm_measurements', []), *(raw_prm_measurements or [])]:
            key = (row['source'], aware_time(row['timestamp']).isoformat())
            if key in prm_index and prm_index[key]['lifetime_us'] != row['lifetime_us']:
                raise ValueError(f'conflicting raw PRM observation: {key}')
            prm_index[key] = row
        data.update(schema_version=5, lifetimes=derived, reference_lifetimes=list(reference_index.values()),
                    raw_prm_measurements=list(prm_index.values()))
        if packet_window_hours is not None:
            data['packet_window_hours'] = packet_window_hours
        if packet_snapshot is not None:
            inventories = data.setdefault('packet_pool_inventories', {})
            period = str(packet_snapshot['window_hours'])
            if packet_snapshot['evaluated_at'] >= inventories.get(period, {}).get('evaluated_at',''):
                inventories[period] = packet_snapshot
        active = data.get('packet_pool_inventories', {}).get(str(data.get('packet_window_hours',24)))
        publication = [r for r in data['lifetimes'] if r.get('method') != 'packet_pool'
                       or (r['window_hours'] == data.get('packet_window_hours',24)
                           and (active is None or list(record_key(r)) in active['keys']))]
        # Per-file packet attempts are diagnostic inputs once pooled windows
        # exist; do not plot both as independent measurements of the same data.
        if any(r.get('method') == 'packet_pool' for r in publication):
            publication = [r for r in publication if r.get('method') != 'packet']
        data['plot_lifetimes'] = aggregate_track_lifetimes(publication, timezone_name)
        # A separate display table keeps these visible exceptions inspectable
        # without letting downstream consumers mistake them for accepted fits.
        data['packet_candidate_points'] = packet_candidate_points(publication)
        data['reference_comparisons'] = compare_reference_lifetimes(data['lifetimes'], data['reference_lifetimes'])
        if annotations is not None:
            data['annotations'] = annotations if isinstance(annotations, list) else json.loads(Path(annotations).read_text())
        if gas_conversion_config is not None:
            data['gas_conversion_config'] = gas_conversion_config
        # Study plots use UTC on both exported surfaces. Keep the grouping
        # clock explicit: switching an axis must never reassign observations.
        data['timezone'] = timezone_name
        data['display_timezone'] = 'UTC'
        atomic_json(path, data)
        if output_timeseries:
            draw_overlay(publication, output_timeseries, timezone_name, data.get('annotations'))
            draw_overlay(publication, raw_plot_path(output_timeseries), timezone_name,
                         data.get('annotations'), raw_points=True, raw_prm_measurements=data['raw_prm_measurements'])
    return data


LABELS = {
    'gas_o2': ('O₂-only gas equivalent · 6 h (provisional)', '#AAB965', 'x'),
    'prm': ('Purity monitor · daily mean', '#2878B5', 'o'),
    'gas': ('O₂ + H₂O gas equivalent · 6 h (provisional)', '#819B28', 's'),
    'beam': ('Externally triggered muon candidates · 6 h mean', '#CB7B21', '^'),
    'cosmic': ('Off-beam / cosmic-enriched candidates · 6 h mean', '#AA4D83', 'v'),
    'all_mip': ('All selected through-going MIP candidates · 6 h mean', '#444444', 'D'),
    'mixed': ('FLOW track candidates · 6 h mean', '#777777', 'o'),
    'packet': ('Packets · whole detector (provisional)', '#28A1A1', 'P'),
}


def overlay_groups(entries, raw_points=False, raw_prm_measurements=None, timezone_name='America/Chicago'):
    """Build identical source-separated series for the PNG and HTML renderers.

    Lines always consume window means. In the diagnostic view, only the PRM
    and track markers change grain: actual PRM timestamps and per-FLOW fits.
    Excluded track fits remain visible as separate markers, never in the line.
    """
    from collections import defaultdict
    groups = defaultdict(lambda: dict(means=[], points=[]))
    def source_key(row):
        sample = row.get('sample', 'mixed')
        return sample, row.get('source', '') if sample in ('prm', 'gas', 'gas_o2', 'packet') else ''
    for row in aggregate_track_lifetimes(entries, timezone_name):
        group = groups[source_key(row)]
        group['means'].append(row)
        if not raw_points or row.get('method', 'track') not in ('track', 'prm'):
            group['points'].append(row)
    if raw_points:
        # record_key removes retries. Reference-only track rows cannot pass
        # is_flow_track, even when their lifetime is positive and plausible.
        for row in {record_key(r): r for r in entries if is_flow_track(r)}.values():
            groups[source_key(row)]['points'].append(row)
        for row in raw_prm_measurements or []:
            groups[source_key(row)]['points'].append(row)
        # The diagnostic view places a mean at the centroid of the actual
        # acquisition times contributing to it, not at a nominal calendar noon.
        # Average epoch seconds explicitly: pandas datetime integer resolution
        # can be us or ns depending on version, so integer casts are unsafe here.
        from purity_sources import aware_time
        import pandas as pd
        for (sample, _), group in groups.items():
            if sample in ('gas', 'gas_o2', 'packet'):
                continue
            positioned = []
            for mean in group['means']:
                mean = dict(mean)
                points = [r for r in group['points'] if valid_lifetime(r) and not r.get('monitoring_excluded')
                          and aware_time(mean['period_start']) <= aware_time(r['timestamp']) < aware_time(mean['period_end'])]
                if points:
                    seconds = [aware_time(r['timestamp']).timestamp() for r in points]
                    mean['bin_midpoint'] = mean['timestamp']
                    mean['timestamp'] = pd.Timestamp(float(np.mean(seconds)), unit='s', tz='UTC').round('us').isoformat()
                    mean['first_observed_at'] = pd.Timestamp(min(seconds), unit='s', tz='UTC').isoformat()
                    mean['last_observed_at'] = pd.Timestamp(max(seconds), unit='s', tz='UTC').isoformat()
                    mean['time_position'] = 'mean acquisition time; horizontal span is first-to-last contributing observation'
                positioned.append(mean)
            group['means'] = positioned
    return dict(groups)


def draw_overlay(entries, output_file, timezone_name='America/Chicago', annotations=None,
                 raw_points=False, raw_prm_measurements=None):
    """Publish a 3000x2000 PNG and matching self-contained interactive plot.

    Raw points are never joined. Connecting only the mean rows makes file-
    level scatter and isolated PRM spikes visible without implying that the
    line traces their acquisition sequence. Gas series keep their six-hour
    conversion grain in both views and stay at the averaging-window midpoint.
    ``timezone_name`` controls grouping only. Both PNG and HTML display true
    UTC instants, including raw readings, mean positions and annotations.
    """
    from collections import defaultdict
    from datetime import datetime, timezone, timedelta
    from zoneinfo import ZoneInfo
    import tempfile
    import matplotlib.dates as mdates
    from purity_sources import aware_time
    import plotly.graph_objects as go

    failed = sum(not valid_lifetime(entry) and entry.get('plot_role') != 'diagnostic_only' for entry in entries)
    series = overlay_groups(entries, raw_points, raw_prm_measurements, timezone_name)
    packet_candidates = packet_candidate_points(entries)
    excluded = sum(bool(r.get('monitoring_excluded')) for r in entries if is_flow_track(r))
    fig, ax = plt.subplots(figsize=(15, 10), dpi=200)
    interactive = go.Figure()
    # Match the UTC instants already emitted to Plotly; previously only the
    # PNG ticks/footer were converted back to the grouping timezone.
    tz = ZoneInfo('UTC')
    packet_legend_shown = False
    for (sample, source), group in sorted(series.items()):
        all_rows = group['means']
        rows = sorted((r for r in group['points'] if valid_lifetime(r)), key=lambda r: aware_time(r['timestamp']))
        if not rows and not any(valid_lifetime(r) for r in all_rows):
            continue
        label, color, marker = LABELS.get(sample, (sample, '#777777', 'o'))
        if sample == 'packet' and rows:
            label = f"Packets · whole detector · {rows[0].get('window_hours','?')} h pooled (provisional)"
        raw_series = raw_points and sample not in ('gas', 'gas_o2', 'packet')
        if raw_series:
            label = label.replace('daily mean', 'readings; □/line: daily mean').replace('6 h mean', 'per-file fits; □/line: 6 h mean')
        if sample in ('gas', 'gas_o2'):
            label = ('O₂ + H₂O' if sample == 'gas' else 'O₂ only') + ' equivalent · 6 h'
            label += f' ({source}'
            if sample == 'gas':
                label += '; ' + str(rows[0].get('water_source') or 'H₂O source unspecified')
            label += ')'
        elif source and sample != 'packet':
            label += f' ({source})'
        # Accepted packet cohorts also remain independent traces. A shared
        # method legend keeps a many-cohort overlay legible; hover text below
        # carries each trace's actual cohort, membership and calibration checks.
        show_series_legend = sample != 'packet' or not packet_legend_shown
        if sample == 'packet':
            packet_legend_shown = True
        open_marker = '1874' in source
        if sample in ('gas', 'gas_o2') and open_marker:
            color = '#8060A6' if sample == 'gas' else '#AD92CC'
        if sample == 'gas_o2' and open_marker:
            marker = '+'
        excluded_rows = [r for r in rows if r.get('monitoring_excluded')]
        rows = [r for r in rows if not r.get('monitoring_excluded')]
        times = [aware_time(row['timestamp']).to_pydatetime().astimezone(tz) for row in rows]
        values = [row['lifetime_us'] for row in rows]
        errors = [row.get('error_us') for row in rows]
        line_times, line_values = connected_values(all_rows)
        linestyle = '--' if sample == 'gas_o2' else '-'
        ax.plot([aware_time(t).to_pydatetime().astimezone(tz) for t in line_times],
                [np.nan if v is None else v for v in line_values],
                linestyle=linestyle, color=color, linewidth=1.5, alpha=0.85)
        ax.plot(times, values, linestyle='none', marker=marker, color=color,
                markersize=3 if raw_series and sample == 'prm' else 5, alpha=0.45 if raw_series and sample == 'prm' else 0.9, label=label if show_series_legend else '_nolegend_',
                markerfacecolor='none' if open_marker else color)
        indices = [i for i, error in enumerate(errors) if error is not None]
        if indices:
            ax.errorbar([times[i] for i in indices], [values[i] for i in indices],
                        yerr=[errors[i] for i in indices], fmt='none', color=color, capsize=2, alpha=0.6)
        symbols = {'beam':'triangle-up', 'cosmic':'triangle-down', 'all_mip':'diamond', 'packet':'cross'}
        symbol = (('square-open' if open_marker else 'square') if sample == 'gas' else
                  ('cross' if open_marker else 'x') if sample == 'gas_o2' else symbols.get(sample, 'circle'))
        from html import escape
        details = [f"{row.get('uncertainty', 'fit covariance uncertainty')}<br>"
                   f"n={row.get('n_measurements', row.get('n_segments', 'unknown'))}<br>"
                   f"{row.get('period_start', '')} — {row.get('period_end', '')}<br>"
                   f"{escape(str(row.get('input_file', '')))}<br>{escape(str(row.get('run_context', '')))}<br>"
                   f"{escape(str(row.get('monitoring_review_reason', '')))}" for row in rows]
        if sample == 'packet':
            details = [detail + f"<br>{escape(source)}<br>Files: {row.get('n_files','?')}; "
                       f"objects: {row.get('n_objects','?')}<br>Minimum pedestal coverage: "
                       f"{format(row['minimum_pedestal_coverage'], '.2%') if row.get('minimum_pedestal_coverage') is not None else 'unavailable'}<br>"
                       f"{escape(row.get('note',''))}" for detail,row in zip(details,rows)]
        interactive.add_trace(go.Scatter(x=[aware_time(t).isoformat() for t in line_times],
            y=line_values, mode='lines', line=dict(color=color, width=1.5, dash='dash' if sample=='gas_o2' else 'solid'),
            legendgroup=label, showlegend=False, hoverinfo='skip', connectgaps=False))
        interactive.add_trace(go.Scatter(x=[aware_time(r['timestamp']).isoformat() for r in rows],
            y=values, mode='markers', name=label, legendgroup=label, showlegend=show_series_legend,
            marker=dict(color=color, symbol=symbol, size=4 if raw_series and sample=='prm' else 7,
                        opacity=0.45 if raw_series and sample=='prm' else 0.9),
            error_y=dict(type='data', array=errors, visible=True), text=details,
            hovertemplate='%{x}<br>%{y:.3f} µs<br>%{text}<extra>%{fullData.name}</extra>'))
        if raw_series:
            # Open squares identify actual means. Horizontal bars are observed
            # coverage, not timing uncertainty or a claim of continuous readout.
            means = [r for r in all_rows if valid_lifetime(r) and r.get('time_position')]
            centers = [aware_time(r['timestamp']) for r in means]
            left = [(t-aware_time(r['first_observed_at'])).total_seconds() for t,r in zip(centers,means)]
            right = [(aware_time(r['last_observed_at'])-t).total_seconds() for t,r in zip(centers,means)]
            ax.errorbar(mdates.date2num([t.to_pydatetime().astimezone(tz) for t in centers]), [r['lifetime_us'] for r in means],
                xerr=np.array([left,right])/86400, fmt='s', markersize=7, markerfacecolor='white',
                markeredgecolor=color, ecolor=color, capsize=3, linestyle='none', zorder=5)
            interactive.add_trace(go.Scatter(x=[t.isoformat() for t in centers],y=[r['lifetime_us'] for r in means],
                mode='markers', name=label+' (mean)', legendgroup=label, showlegend=False,
                marker=dict(color=color, symbol='square-open', size=10),
                error_x=dict(type='data', symmetric=False, array=np.array(right)*1000, arrayminus=np.array(left)*1000),
                text=[f"Mean of {r.get('n_measurements','?')} observations<br>Window: {r['period_start']} — {r['period_end']}<br>"
                      f"Observed: {r['first_observed_at']} — {r['last_observed_at']}<br>Mean error: {r.get('error_us')} µs" for r in means],
                hovertemplate='%{x}<br>Mean %{y:.3f} µs<br>%{text}<extra></extra>'))
        if excluded_rows:
            # Preserve excluded measurements for inspection. A hollow diamond
            # marks a reviewed omission, not a failed fit or a zero lifetime.
            excluded_label = label.split(' ·')[0] + ' · excluded from mean'
            ax.plot([aware_time(r['timestamp']).to_pydatetime().astimezone(tz) for r in excluded_rows],
                    [r['lifetime_us'] for r in excluded_rows], linestyle='none', marker='D',
                    markerfacecolor='none', markeredgecolor=color, markersize=8, label=excluded_label)
            ax.errorbar([aware_time(r['timestamp']).to_pydatetime().astimezone(tz) for r in excluded_rows if r.get('error_us') is not None],
                        [r['lifetime_us'] for r in excluded_rows if r.get('error_us') is not None],
                        yerr=[r['error_us'] for r in excluded_rows if r.get('error_us') is not None],
                        fmt='none', color=color, capsize=2, alpha=0.6)
            interactive.add_trace(go.Scatter(x=[aware_time(r['timestamp']).isoformat() for r in excluded_rows],
                y=[r['lifetime_us'] for r in excluded_rows], mode='markers', name=excluded_label,
                marker=dict(color=color, symbol='diamond-open', size=10),
                error_y=dict(type='data', array=[r.get('error_us') for r in excluded_rows], visible=True),
                text=[escape(r.get('input_file', ''))+'<br>'+escape(r['monitoring_review_reason']) for r in excluded_rows],
                hovertemplate='%{x}<br>%{y:.3f} µs<br>%{text}<extra>Excluded from mean</extra>'))

    # Show failed-check packet candidates as a separate visual class. Source
    # keys already encode cohort, calibration/selection, and window length,
    # so dashed guides cannot join incompatible samples or different periods.
    candidate_groups = defaultdict(list)
    for row in packet_candidates:
        candidate_groups[row['source']].append(row)
    for group_index, (source, rows) in enumerate(sorted(candidate_groups.items())):
        from html import escape
        rows.sort(key=lambda r: aware_time(r['timestamp']))
        color = LABELS['packet'][1]
        # Many configuration cohorts can share a method/date. Give the method
        # one compact legend entry, while preserving separate traces/lines and
        # complete cohort identities in hover text. Never join cohorts to make
        # a continuous-looking packet curve or crowd the figure with hashes.
        label = f"Packet fit candidates · {rows[0]['window_hours']} h · checks failed"
        times = [aware_time(r['timestamp']).to_pydatetime().astimezone(tz) for r in rows]
        values = [r['display_lifetime_us'] for r in rows]
        errors = [r['display_error_us'] for r in rows]
        ax.plot(times, values, '--', color=color, linewidth=1.5, alpha=.85)
        # PathCollection supports dashed marker outlines in the PNG. Plotly
        # uses its native open hexagon plus a dashed connecting line instead.
        ax.scatter(times, values, marker='h', s=90, facecolors='none', edgecolors=color,
                   linewidths=1.5, linestyles='--', label=label if group_index == 0 else '_nolegend_', zorder=6)
        measured = [i for i,e in enumerate(errors) if e is not None]
        if measured:
            ax.errorbar([times[i] for i in measured], [values[i] for i in measured],
                        yerr=[errors[i] for i in measured], fmt='none', ecolor=color, capsize=3, alpha=.7)
        details = []
        for row in rows:
            coverage = row.get('minimum_pedestal_coverage')
            coverage_text = f'{coverage:.2%}' if coverage is not None else 'unavailable'
            uncertainty = ('Conditional fit-covariance error only; missing-calibration/selection effects excluded'
                           if row['display_error_us'] is not None else 'Candidate uncertainty unavailable')
            details.append(f"Diagnostic candidate — checks failed<br>{escape(row.get('fit_message',''))}<br>"
                f"Minimum pedestal coverage: {coverage_text}<br>Files: {row.get('n_files','?')}; objects: {row.get('n_objects','?')}<br>"
                f"Window: {row.get('period_start','')} — {row.get('period_end','')}<br>"
                f"File starts: {row.get('first_observed_at','')} — {row.get('last_observed_at','')}<br>"
                f"{uncertainty}<br>{escape(source)}")
        interactive.add_trace(go.Scatter(x=[aware_time(r['timestamp']).isoformat() for r in rows],
            y=values, mode='lines+markers', name=label, legendgroup='packet-candidates', showlegend=group_index == 0,
            line=dict(color=color, dash='dash', width=1.5),
            marker=dict(color=color, symbol='hexagon-open', size=13, line=dict(width=2)),
            error_y=dict(type='data', array=errors, visible=True), text=details,
            hovertemplate='%{x}<br>Candidate %{y:.3f} µs<br>%{text}<extra>%{fullData.name}</extra>'))
    title = 'DUNE ND Prototype 2×2 · liquid argon purity' + (' · individual observations' if raw_points else '')
    fig.suptitle(title, fontsize=15)
    ax.set(xlabel='Date (UTC)', ylabel='Electron lifetime [µs]')
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18)
    # Include individual timestamps when setting the axis: otherwise an
    # excluded latest file or an early raw PRM reading could be clipped away.
    groups = {key: [r for r in [*g['means'], *g['points']] if valid_lifetime(r)] for key, g in series.items()}
    groups.update({('packet_candidate',source):rows for source,rows in candidate_groups.items()})
    groups = {key: rows for key, rows in groups.items() if rows}
    if groups:
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.015), fontsize=8,
                  framealpha=0.9, ncol=2 if len(groups)>3 else 1)
        observed = [aware_time(r['timestamp']).to_pydatetime().astimezone(tz)
                    for rows in groups.values() for r in rows]
        first = min(observed).replace(hour=0, minute=0, second=0, microsecond=0)
        last = max(observed).replace(hour=0, minute=0, second=0, microsecond=0)
        days = (last.date()-first.date()).days
        step = max(1, int(np.ceil(days/6)))
        ticks = [first+timedelta(days=d) for d in range(0, days+1, step)]
        if ticks[-1] != last:
            if len(ticks)>1 and (last-ticks[-1]).total_seconds()/86400 < step/2:
                ticks.pop()
            ticks.append(last)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d', tz=tz))
    else:
        ax.text(0.5, 0.5, 'No valid lifetime measurements available', transform=ax.transAxes, ha='center')
    if annotations:
        events = annotations if isinstance(annotations, list) else json.loads(Path(annotations).read_text())
        for event in events:
            stamp = aware_time(event['timestamp']).to_pydatetime()
            if 'end' in event:
                end = aware_time(event['end']).to_pydatetime()
                ax.axvspan(stamp, end, color='grey', alpha=0.12)
                center = stamp+(end-stamp)/2
                ax.text(center, 0.96, event['label'], transform=ax.get_xaxis_transform(),
                        ha='center', va='top', fontsize=9, color='#555555')
                interactive.add_shape(type='rect', xref='x', yref='paper', x0=stamp.isoformat(),
                    x1=end.isoformat(), y0=0, y1=1, fillcolor='grey', opacity=0.12, line_width=0)
                interactive.add_annotation(x=center.isoformat(), y=0.96, xref='x', yref='paper',
                    text=event['label'], showarrow=False, yanchor='top')
            else:
                ax.axvline(stamp, color='#333333', linestyle='--', linewidth=1.2, alpha=0.8)
                ax.text(stamp, 0.89, event['label'], transform=ax.get_xaxis_transform(), rotation=90,
                        va='top', ha='right', fontsize=9, color='#333333')
                interactive.add_shape(type='line', xref='x', yref='paper', x0=stamp.isoformat(),
                    x1=stamp.isoformat(), y0=0, y1=1, line=dict(color='#333333', dash='dash'))
                interactive.add_annotation(x=stamp.isoformat(), y=0.89, xref='x', yref='paper',
                    text=event['label'], showarrow=False, textangle=-90, yanchor='top')
    accepted_groups = {key: [r for r in rows if valid_lifetime(r) and not r.get('monitoring_excluded')] for key, rows in groups.items()}
    latest = '; '.join(f"{LABELS.get(s, (s,))[0].split(' ·')[0]}: "
                       f"{max(aware_time(r.get('last_observed_at') or r['timestamp']) for r in rows):%m-%d %H:%M %Z}"
                       for (s, _), rows in sorted(accepted_groups.items()) if rows)
    point_note = ('PRM: individual readings (errors unavailable). Tracks: individual FLOW fits with fit errors. Lines join daily PRM / 6 h track means.\n'
                  if raw_points else
                  'PRM bars: daily SEM. Tracks: 6 h file means; bars = max(propagated fit error, between-file SEM), where available.\n')
    note = (point_note +
            'Gas conversion/packet selections provisional. Lines connect displayed points across gaps as visual guides, not interpolated measurements.\n'
            'External triggers do not establish beam origin. PRM means use recorded positive lifetimes; waveform quality is not verified.\n'
            f'Failed/unavailable entries: {failed}; reviewed track exclusions: {excluded} (retained in JSON / raw plot). Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}.')
    gas_fields = sorted({r['gas_conversion']['field_strength_vpercm'] for r in entries if r.get('gas_conversion')})
    if gas_fields:
        note += '\nGas: Eva 2026-10-06 rates at ' + ', '.join(f'{v:g}' for v in gas_fields) + ' V/cm reference field; raw DB concentrations; offset/phase mapping unverified.'
    packet_windows = [r for r in entries if r.get('method') == 'packet_pool']
    if packet_windows and not any(valid_lifetime(r) for r in packet_windows):
        note += f"\nPackets: no accepted {packet_windows[0]['window_hours']} h pooled estimate; drift-bin fit and calibration gates are recorded in JSON."
    if packet_candidates:
        note += ('\nOpen teal hexagons / dashed guides: packet fit candidates with failed checks. Bars, where present: conditional fit covariance only.\n'
                 'Candidate bars exclude missing-calibration and selection effects; finite candidates are displayed without changing fit acceptance.')
    if any(r.get('method') == 'packet_pool_refresh' and r.get('fit_status') == 'failed' for r in entries):
        note += '\nPacket refresh failed: retained packet history may be stale; inspect the refresh diagnostic in JSON.'
    if raw_points and not raw_prm_measurements:
        note += '\nRaw PRM readings unavailable: import the source snapshots; daily means are not raw observations.'
    if raw_points:
        note += '\nOpen squares: means at mean acquisition times. Horizontal bars: first–last contributing observations, not full-window coverage.'
    import textwrap
    fig.text(0.08, 0.045, note + '\nLatest accepted data: ' + '\n'.join(textwrap.wrap(latest or 'none', 145)), fontsize=8, va='bottom')
    fig.tight_layout(rect=(0, (0.22 if gas_fields else 0.20) if packet_candidates else 0.16, 1, 0.96))
    target = Path(output_file)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(suffix='.png', dir=target.parent)
    os.close(fd)
    try:
        fig.savefig(temporary, dpi=200)
        public_permissions(temporary)
        os.replace(temporary, target)
    finally:
        plt.close(fig)
        if os.path.exists(temporary):
            os.unlink(temporary)
    interactive.update_layout(title=title,
        xaxis_title='Timestamp (UTC)', yaxis_title='Electron lifetime [µs]',
        yaxis_rangemode='tozero', template='plotly_white', width=1500, height=1000)
    interactive.add_annotation(text=note.replace('\n', '<br>'), x=0, y=-0.20, xref='paper', yref='paper',
                               showarrow=False, align='left', xanchor='left', yanchor='top', font=dict(size=11))
    interactive.update_layout(margin=dict(b=310 if packet_candidates else 220), legend=dict(font=dict(size=11)))
    fd, temporary = tempfile.mkstemp(suffix='.html', dir=target.parent)
    os.close(fd)
    try:
        interactive.write_html(temporary, include_plotlyjs=True)
        public_permissions(temporary)
        os.replace(temporary, str(target)+'.html')
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def packet_result(args):
    """Extract packet objects in a fresh directory for subsequent window fits.

    A new attempt directory prevents a failed subprocess from reusing a stale
    summary. Packet failure is independent of track/slow-controls availability.
    Legacy standalone summaries remain importable with null lifetime error;
    new pooled fits derive a conditional error from attenuation covariance.
    """
    import subprocess
    import tempfile
    from purity_sources import aware_time
    stamp = args.timestamp or date_from_filename(args.input_file or args.packet_file).isoformat()
    result = dict(timestamp=aware_time(stamp).isoformat(), sample='packet', method='packet',
        input_file=os.path.abspath(args.input_file or args.packet_file or args.packet_summary), fit_status='failed',
        lifetime_us=None, error_us=None, fit_quality='provisional_selection',
        uncertainty='not estimated by the reference packet fitter')
    try:
        if args.packet_summary:
            summary_path = Path(args.packet_summary)
        else:
            output = Path(args.output_file_json).parent / 'packet_diagnostics'
            output.mkdir(parents=True, exist_ok=True)
            # A unique attempt directory avoids accepting stale results after a failed retry.
            output = Path(tempfile.mkdtemp(prefix=Path(args.packet_file).stem+'.', dir=output))
            command = [sys.executable, '-u', str(Path(__file__).with_name('packet_lifetime.py')),
                args.packet_file, '--outdir', str(output), '--detector-wide',
                '--charge-mode', 'adc-minus-ped', '--successive-min-ticks', str(args.successive_ticks),
                '--successive-max-ticks', str(args.successive_ticks), '--mpv-bootstrap', str(args.mpv_bootstrap),
                '--sample-directory', args.packet_samples, '--extract-only', '--sample-timestamp', stamp]
            if args.packet_cohort:
                command += ['--cohort', args.packet_cohort]
            if args.pedestal:
                command += ['--ped', args.pedestal, '--ped-source', args.packet_pedestal_source or 'panel']
                if args.packet_pedestal_source == 'calibration':
                    if not args.input_file: raise ValueError('static calibration requires matching FLOW')
                    command += ['--flow-file', args.input_file]
            elif args.input_file:
                command += ['--ped-source', 'flow', '--flow-file', args.input_file]
            else:
                raise ValueError('packets require a matching FLOW file or explicit panel pedestal')
            if args.ext_io is not None:
                command += ['--ext-io', str(args.ext_io)]
            with (output/'analysis.log').open('w') as log:
                subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
            summary_path = output/'summary.json'
        summary = json.loads(summary_path.read_text())
        if summary.get('extraction_only'):
            result.update(fit_status='extracted', plot_role='diagnostic_only',
                          summary_file=str(summary_path.resolve()), n_objects=summary['sample']['n_objects'],
                          fit_message='Charge objects saved; lifetime estimated in fixed pooled windows')
            return result
        if summary.get('aggregation') != 'detector-wide pooled charge':
            raise ValueError('packet summary is not a detector-wide pooled fit')
        tau = summary.get('weighted_tau_ms') or summary.get('unweighted_tau_ms')
        if tau is None or not np.isfinite(tau) or tau <= 0:
            raise ValueError('no positive finite packet lifetime')
        result.update(lifetime_us=1000*float(tau), fit_status='ok',
            summary_file=str(summary_path.resolve()), packet_file=summary['packet_file'],
            n_objects=summary['n_fit_objects'], estimator='weighted' if summary.get('weighted_tau_ms') else 'unweighted',
            packet_configuration={k: summary.get(k) for k in ('ext_io_group', 'ext_trigger_type',
                'charge_convention', 'successive_hit_grouping', 'pedestal_coverage', 'aggregation')})
    except Exception as exc:
        result['fit_message'] = f'{type(exc).__name__}: {exc}'
        if 'summary_path' in locals():
            result['summary_file'] = str(summary_path.resolve())
        if 'summary' in locals():
            result['pedestal_coverage'] = summary.get('pedestal_coverage')
            result['n_objects'] = summary.get('n_fit_objects')
        print('Packet lifetime unavailable:', result['fit_message'])
    return result


def main():
    """Dispatch export or analysis, then publish both views of one history."""
    parser = argparse.ArgumentParser(description='Fit, aggregate and overlay all 2x2 purity methods.')
    parser.add_argument('--monitor-config', help='Configured six-hour pooled watcher profile; owns state and publication paths')
    parser.add_argument('--input_file', '--input-file', help='FLOW file; optional for slow-controls-only refresh')
    parser.add_argument('--output_file_plot', '--output-file-plot', help='Per-file track diagnostic plot')
    parser.add_argument('--output_file_json', '--output-file-json', help='Shared lifetime history JSON')
    parser.add_argument('--output-timeseries', help='Mean overlay PNG; also publishes _raw.png and both .html companions')
    parser.add_argument('--sample', action='append', type=parse_sample_spec)
    parser.add_argument('--flow-cache', help='Optional persistent segment cache for remote FLOW; default streams into memory')
    parser.add_argument('--reference-history', action='append', default=[], help='Published JSON for corroboration only; never supplies plotted estimates')
    parser.add_argument('--select-muons', action='store_true', help='Run installed ndlar_flow selection into a sidecar when FLOW has no segments')
    parser.add_argument('--write-hdf5-metadata', action='store_true')
    parser.add_argument('--slow-controls', action='append', default=[], help='Raw snapshot JSON; repeat for multiple files')
    parser.add_argument('--timezone', default='America/Chicago', help='Calendar-day and six-hour bin timezone')
    parser.add_argument('--gas-quality-config', help='JSON overrides for multi-hour coverage, plateaus and calibration limits')
    parser.add_argument('--gas-conversion-config', help='Explicit Eva field-dependent gas model JSON; saved in history for subsequent refreshes')
    parser.add_argument('--annotations', help='JSON list of timestamp/label and optional end spans')
    parser.add_argument('--track-review', help='JSON rules with exact input_file, optional sample, include/exclude action and reason; persists in history')
    parser.add_argument('--timestamp', help='Explicit offset-aware timestamp if filename has no timestamp')
    parser.add_argument('--packet-file', help='Native raw packet HDF5 paired with the input FLOW')
    parser.add_argument('--packet-summary', help='Previously generated detector-wide packet summary.json')
    parser.add_argument('--packet-samples', help='Persistent compact packet shard directory; defaults beside history')
    parser.add_argument('--packet-window-hours', type=int, choices=[6,24,48], default=24,
                        help='Fixed packet pooling interval; 24 h is an operational default, not a validated optimum')
    parser.add_argument('--packet-pool-bootstrap', type=int, default=100)
    parser.add_argument('--packet-cohort', help='Acquisition/configuration cohort; incompatible cohorts never pooled')
    parser.add_argument('--packet-pedestal-source', choices=['panel','calibration'],
                        help='Interpret --pedestal as panel means or the exact static calibration named by FLOW')
    parser.add_argument('--pedestal', help='Pedestal JSON interpreted by --packet-pedestal-source; default without JSON recovers FLOW Q_raw baseline')
    parser.add_argument('--ext-io', type=int)
    parser.add_argument('--successive-ticks', type=int, default=28)
    parser.add_argument('--mpv-bootstrap', type=int, default=20)
    parser.add_argument('--export-slow-controls', metavar='OUTPUT', help='Export database observations atomically; no plotting')
    parser.add_argument('--source-config', help='Database source config; credentials via environment variables')
    parser.add_argument('--start', help='Inclusive export start, ISO timestamp with offset')
    parser.add_argument('--end', help='Exclusive export end, ISO timestamp with offset')
    args = parser.parse_args()
    if args.monitor_config:
        # Configured and legacy publication have different state/membership
        # contracts. Reject conflicting knobs instead of silently ignoring an
        # export or legacy selection request when the profile is enabled.
        if any(value != parser.get_default(name) for name, value in vars(args).items()
               if name not in ('monitor_config', 'input_file')):
            parser.error('--monitor-config owns paths and pairing; supply only an optional --input-file')
        from monitor import run_monitor
        run_monitor(args.monitor_config, args.input_file)
        return
    from purity_sources import aggregate_measurements, export_snapshot, aware_time, load_measurements, prm_observations
    if args.export_slow_controls:
        if not all((args.source_config, args.start, args.end)):
            parser.error('export requires --source-config, --start and --end')
        export_snapshot(args.source_config, args.start, args.end, args.export_slow_controls)
        return
    if not args.output_file_json:
        parser.error('--output_file_json is required')
    args.packet_samples = args.packet_samples or str(Path(args.output_file_json).parent/'packet_samples')
    if args.packet_pool_bootstrap < 10:
        parser.error('--packet-pool-bootstrap requires at least 10 trials')
    if args.input_file and not args.output_file_plot:
        parser.error('FLOW fitting requires --output_file_plot')
    if args.packet_summary and not (args.timestamp or args.input_file):
        parser.error('--packet-summary without FLOW requires --timestamp')
    from flow_input import is_remote
    if args.input_file and is_remote(args.input_file) and (args.packet_file or args.write_hdf5_metadata):
        parser.error('remote FLOW segment access supports track fits only; packet recovery/metadata writes need a local full FLOW')
    if args.packet_file and args.packet_summary:
        parser.error('choose --packet-file or --packet-summary')
    results = []
    if args.input_file:
        try:
            fit_input = args.input_file
            if is_remote(args.input_file) and args.flow_cache:
                from flow_input import cache_flow_segments
                fit_input = cache_flow_segments(args.input_file, args.flow_cache or str(Path(args.output_file_json).parent/'flow_segments'),
                                               [path for _, path in args.sample] if args.sample else None)
            if args.select_muons:
                from flow_input import open_flow
                with open_flow(fit_input) as (source, _):
                    needs_selection = not any(path in source for _,path in (*SPLIT_SAMPLES, LEGACY_SAMPLE))
                if needs_selection:
                    if is_remote(fit_input):
                        raise ValueError('remote FLOW has no segments; full event selection requires a local FLOW file')
                    from track_selection import select_muons
                    fit_input = str(Path(args.output_file_plot).with_suffix('.segments.h5'))
                    select_muons(args.input_file, fit_input)
            stamp = args.timestamp or date_from_filename(args.input_file).isoformat()
            fitted = fit_samples(fit_input, args.output_file_plot, args.sample, stamp)
            if fit_input != args.input_file and not is_remote(args.input_file):
                for row in fitted:
                    row['selection_file'] = os.path.abspath(fit_input)
                    row['input_file'] = os.path.abspath(args.input_file)
            results.extend(fitted)
        except Exception as exc:
            stamp = args.timestamp or date_from_filename(args.input_file).isoformat()
            results.append(dict(timestamp=aware_time(stamp).isoformat(), sample='all_mip', method='track',
                input_file=args.input_file if is_remote(args.input_file) else os.path.abspath(args.input_file),
                calculation_source='flow_segments', fit_status='unavailable',
                lifetime_us=None, error_us=None, fit_message=f'{type(exc).__name__}: {exc}'))
            print('Track sample unavailable:', exc)
        if args.write_hdf5_metadata:
            write_hdf5_metadata(args.input_file, [r for r in results if isinstance(r.get('segments_dset'), str) and 'selection_file' not in r])
    if args.packet_file or args.packet_summary:
        results.append(packet_result(args))
    # Refresh packet windows even when this call was triggered only by a new
    # slow-controls snapshot. Unchanged windows reuse their audited fit cache.
    from packet_pooling import pool_samples
    from datetime import datetime, timezone
    packet_snapshot = {}
    try:
        results.extend(pool_samples(args.packet_samples, args.packet_window_hours, args.timezone,
                                    bootstrap=args.packet_pool_bootstrap, snapshot=packet_snapshot))
        if packet_snapshot:
            # Replace a previous refresh failure for this directory after a
            # successful recovery. This status is never a lifetime estimate.
            results.append(dict(method='packet_pool_refresh', sample='packet',
                input_file=str(Path(args.packet_samples).resolve()), timestamp=datetime.now(timezone.utc).isoformat(),
                fit_status='refreshed', plot_role='diagnostic_only', lifetime_us=None, error_us=None))
    except Exception as exc:
        # A corrupt shard or pooling failure must remain explicit while the
        # independent PRM/gas/track sources can still be refreshed for shifters.
        packet_snapshot = None
        results.append(dict(method='packet_pool_refresh', sample='packet',
            input_file=str(Path(args.packet_samples).resolve()), timestamp=datetime.now(timezone.utc).isoformat(),
            fit_status='failed', plot_role='diagnostic_only', lifetime_us=None, error_us=None,
            fit_message=f'{type(exc).__name__}: {exc}'))
        print('Packet window refresh unavailable:', exc)
    raw_prm = None
    conversion_config = resolve_gas_conversion_config(args.gas_conversion_config, args.output_file_json)
    if args.slow_controls:
        # Normalize once, feeding both the means and the individual PRM view
        # from the exact same accepted/deduplicated observations.
        measurements = load_measurements(args.slow_controls)
        raw_prm = prm_observations(measurements)
        results.extend(aggregate_measurements(args.slow_controls, args.timezone,
            json.loads(Path(args.gas_quality_config).read_text()) if args.gas_quality_config else None,
            measurements=measurements, conversion_config=conversion_config))
    update_json(args.output_file_json, results, args.output_timeseries, args.timezone, args.annotations, args.reference_history,
                raw_prm_measurements=raw_prm, track_review=args.track_review,
                packet_window_hours=args.packet_window_hours, packet_snapshot=packet_snapshot or None,
                gas_conversion_config=conversion_config)
    if not results:
        print('Refreshed existing history; no new observations supplied.')
    elif not any(r.get('fit_status') == 'ok' for r in results):
        print('No new valid lifetimes. Unavailable/failed results recorded; history plot refreshed.')


if __name__ == '__main__':
    main()
