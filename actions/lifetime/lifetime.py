#!/usr/bin/env python3
import argparse
import json
import os
import re
import sys
import matplotlib
matplotlib.use("Agg")
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from flufl.lock import Lock

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nearline_util import date_from_filename


SPLIT_SAMPLES = (
    ('beam', 'analysis/beam_rock_muon_segments/data'),
    ('cosmic', 'analysis/cosmic_muon_segments/data'),
)
LEGACY_SAMPLE = ('all_mip', 'analysis/rock_muon_segments/data')


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
    """Resolve explicit samples, split beam/cosmic samples, or the legacy sample."""
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
    from purity_sources import aware_time
    timestamp = aware_time(timestamp) if timestamp else date_from_filename(input_file)
    results = []

    print(f'Opening file: {input_file}')
    with h5py.File(input_file, 'r') as h5_file:
        samples = find_samples(h5_file, requested_samples)
        multiple_samples = len(samples) > 1

        for sample, dset_path in samples:
            paths = dset_path if isinstance(dset_path, list) else [dset_path]
            segments = valid_segments(np.concatenate([h5_file[path][:] for path in paths]))
            print(f'Extracting {sample} lifetime using {len(segments)} segments')

            result = {
                'timestamp': timestamp.isoformat(),
                'sample': sample,
                'method': 'track',
                'uncertainty': 'DeMario fit covariance; excludes selection/calibration systematics',
                'input_file': os.path.abspath(input_file),
                'segments_dset': dset_path,
                'n_segments': len(segments),
                'fit_status': 'failed',
                'lifetime_us': None,
                'error_us': None,
            }

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
    method = entry.get('method', 'track')
    if method in ('prm', 'gas', 'gas_o2'):
        return method, entry.get('source'), entry.get('period_start', entry['timestamp'])
    return method, entry.get('input_file', entry['timestamp']), entry.get('sample', 'mixed')


def update_json(output_file_json, results, output_timeseries=None, timezone_name='America/Chicago', annotations=None):
    """Serialize read/update/render so concurrent workers cannot publish stale history."""
    from datetime import timedelta
    from purity_sources import atomic_json
    path = Path(output_file_json)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Lock(str(path) + '.lock', default_timeout=timedelta(seconds=60), lifetime=timedelta(minutes=10)):
        data = json.loads(path.read_text()) if path.exists() else {'lifetimes': []}
        indexed = {record_key(row): row for row in data['lifetimes']}
        indexed.update({record_key(row): row for row in results})
        data.update(schema_version=2, lifetimes=list(indexed.values()))
        atomic_json(path, data)
        if output_timeseries:
            draw_overlay(data['lifetimes'], output_timeseries, timezone_name, annotations)
    return data


LABELS = {
    'gas_o2': ('O₂-only gas equivalent · 6 h (provisional)', '#AAB965', 'x'),
    'prm': ('Purity monitor · daily mean', '#2878B5', 'o'),
    'gas': ('Gas analyser · 6 h equivalent (provisional)', '#819B28', 's'),
    'beam': ('Externally triggered muon candidates', '#CB7B21', '^'),
    'cosmic': ('Off-beam / cosmic-enriched candidates', '#AA4D83', 'v'),
    'all_mip': ('All selected through-going MIP candidates', '#444444', 'D'),
    'mixed': ('Published mixed muon candidates', '#777777', 'o'),
    'packet': ('Packets · whole detector (provisional)', '#28A1A1', 'P'),
}


def draw_overlay(entries, output_file, timezone_name='America/Chicago', annotations=None):
    """Publish exact 3000x2000 PNG and a self-contained interactive companion."""
    from collections import defaultdict
    from datetime import datetime, timezone
    from zoneinfo import ZoneInfo
    import tempfile
    import matplotlib.dates as mdates
    from purity_sources import aware_time
    import plotly.graph_objects as go

    groups = defaultdict(list)
    failed = 0
    for entry in entries:
        value, error = entry.get('lifetime_us'), entry.get('error_us')
        if (entry.get('fit_status', 'ok') != 'ok' or value is None or
                not np.isfinite(value) or value <= 0 or
                (error is not None and (not np.isfinite(error) or error < 0))):
            failed += 1
            continue
        sample = entry.get('sample', 'mixed')
        groups[(sample, entry.get('source', '') if sample in ('prm', 'gas', 'gas_o2') else '')].append(entry)
    fig, ax = plt.subplots(figsize=(15, 10), dpi=200)
    interactive = go.Figure()
    tz = ZoneInfo(timezone_name)
    for (sample, source), rows in sorted(groups.items()):
        rows.sort(key=lambda r: aware_time(r['timestamp']))
        label, color, marker = LABELS.get(sample, (sample, '#777777', 'o'))
        if source:
            label += f' ({source})'
        open_marker = '1874' in source
        if sample == 'gas_o2' and open_marker:
            marker = '+'
        times = [aware_time(row['timestamp']).to_pydatetime().astimezone(tz) for row in rows]
        values = [row['lifetime_us'] for row in rows]
        errors = [row.get('error_us') for row in rows]
        # Markers do not interpolate through missing measurements or fit failures.
        ax.plot(times, values, linestyle='none', marker=marker, color=color,
                markersize=5, alpha=0.85, label=label,
                markerfacecolor='none' if open_marker else color)
        indices = [i for i, error in enumerate(errors) if error is not None]
        if indices:
            ax.errorbar([times[i] for i in indices], [values[i] for i in indices],
                        yerr=[errors[i] for i in indices], fmt='none', color=color, capsize=2, alpha=0.6)
        details = [f"{row.get('uncertainty', 'fit covariance uncertainty')}<br>"
                   f"{row.get('fit_quality', '')}<br>"
                   f"n={row.get('n_segments', row.get('n_measurements', 'unknown'))}" for row in rows]
        interactive.add_trace(go.Scatter(x=[aware_time(r['timestamp']).isoformat() for r in rows],
            y=values, mode='markers', name=label, marker=dict(color=color, symbol=('square-open' if open_marker else 'square') if sample == 'gas' else ('cross' if open_marker else 'x') if sample == 'gas_o2' else 'circle'),
            error_y=dict(type='data', array=errors, visible=True), text=details,
            hovertemplate='%{x}<br>%{y:.3f} µs<br>%{text}<extra>%{fullData.name}</extra>'))
    fig.suptitle('DUNE ND Prototype 2×2 · liquid argon purity', fontsize=15)
    ax.set(xlabel=f'Date ({timezone_name})', ylabel='Electron lifetime [µs]')
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.18)
    if groups:
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.015), fontsize=8,
                  framealpha=0.9, ncol=2 if len(groups)>3 else 1)
        locator = mdates.AutoDateLocator(tz=tz, minticks=4, maxticks=8)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator, tz=tz))
    else:
        ax.text(0.5, 0.5, 'No valid lifetime measurements available', transform=ax.transAxes, ha='center')
    if annotations:
        for event in json.loads(Path(annotations).read_text()):
            stamp = aware_time(event['timestamp']).to_pydatetime()
            if 'end' in event:
                ax.axvspan(stamp, aware_time(event['end']).to_pydatetime(), color='grey', alpha=0.15)
            else:
                ax.axvline(stamp, color='grey', linestyle=':', alpha=0.7)
            ax.text(stamp, 0.97, event['label'], transform=ax.get_xaxis_transform(), rotation=90, va='top', fontsize=8)
    latest = '; '.join(f"{LABELS.get(s, (s,))[0].split(' ·')[0]}: "
                       f"{max(aware_time(r.get('last_observed_at') or r['timestamp']) for r in rows).tz_convert(timezone_name):%m-%d %H:%M %Z}"
                       for (s, _), rows in sorted(groups.items()))
    note = ('PRM bars: daily SEM. Track bars: fit covariance. Gas/packet bars absent when unavailable.\n'
            'Gas conversion and packet selections are provisional. No interpolation across missing data.\n'
            f'Failed/unavailable entries retained in JSON: {failed}. Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}.')
    import textwrap
    fig.text(0.08, 0.045, note + '\nLatest accepted data: ' + '\n'.join(textwrap.wrap(latest, 145)), fontsize=8, va='bottom')
    fig.tight_layout(rect=(0, 0.16, 1, 0.96))
    target = Path(output_file)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(suffix='.png', dir=target.parent)
    os.close(fd)
    try:
        fig.savefig(temporary, dpi=200)
        os.chmod(temporary, 0o644)
        os.replace(temporary, target)
    finally:
        plt.close(fig)
        if os.path.exists(temporary):
            os.unlink(temporary)
    interactive.update_layout(title='DUNE ND Prototype 2×2 · liquid argon purity',
        xaxis_title='Timestamp (UTC)', yaxis_title='Electron lifetime [µs]',
        yaxis_rangemode='tozero', template='plotly_white', width=1500, height=1000)
    fd, temporary = tempfile.mkstemp(suffix='.html', dir=target.parent)
    os.close(fd)
    try:
        interactive.write_html(temporary, include_plotlyjs=True)
        os.chmod(temporary, 0o644)
        os.replace(temporary, str(target)+'.html')
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def packet_result(args):
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
                '--successive-max-ticks', str(args.successive_ticks), '--mpv-bootstrap', str(args.mpv_bootstrap)]
            if args.pedestal:
                command += ['--ped', args.pedestal, '--ped-source', 'panel']
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
    parser = argparse.ArgumentParser(description='Fit, aggregate and overlay all 2x2 purity methods.')
    parser.add_argument('--input_file', '--input-file', help='FLOW file; optional for slow-controls-only refresh')
    parser.add_argument('--output_file_plot', '--output-file-plot', help='Per-file track diagnostic plot')
    parser.add_argument('--output_file_json', '--output-file-json', help='Shared lifetime history JSON')
    parser.add_argument('--output-timeseries', help='3000x2000 overlay PNG (plus .html)')
    parser.add_argument('--sample', action='append', type=parse_sample_spec)
    parser.add_argument('--select-muons', action='store_true', help='Run installed ndlar_flow selection into a sidecar when FLOW has no segments')
    parser.add_argument('--write-hdf5-metadata', action='store_true')
    parser.add_argument('--slow-controls', action='append', default=[], help='Raw snapshot JSON; repeat for multiple files')
    parser.add_argument('--timezone', default='America/Chicago', help='Calendar-day and six-hour bin timezone')
    parser.add_argument('--gas-quality-config', help='JSON overrides for multi-hour coverage, plateaus and calibration limits')
    parser.add_argument('--annotations', help='JSON list of timestamp/label and optional end spans')
    parser.add_argument('--timestamp', help='Explicit offset-aware timestamp if filename has no timestamp')
    parser.add_argument('--packet-file', help='Native raw packet HDF5 paired with the input FLOW')
    parser.add_argument('--packet-summary', help='Previously generated detector-wide packet summary.json')
    parser.add_argument('--pedestal', help='Panel-pedestal JSON; default recovers static baseline from FLOW Q_raw')
    parser.add_argument('--ext-io', type=int)
    parser.add_argument('--successive-ticks', type=int, default=28)
    parser.add_argument('--mpv-bootstrap', type=int, default=20)
    parser.add_argument('--export-slow-controls', metavar='OUTPUT', help='Export database observations atomically; no plotting')
    parser.add_argument('--source-config', help='Database source config; credentials via environment variables')
    parser.add_argument('--start', help='Inclusive export start, ISO timestamp with offset')
    parser.add_argument('--end', help='Exclusive export end, ISO timestamp with offset')
    args = parser.parse_args()
    from purity_sources import aggregate_measurements, export_snapshot, aware_time
    if args.export_slow_controls:
        if not all((args.source_config, args.start, args.end)):
            parser.error('export requires --source-config, --start and --end')
        export_snapshot(args.source_config, args.start, args.end, args.export_slow_controls)
        return
    if not args.output_file_json:
        parser.error('--output_file_json is required')
    if args.input_file and not args.output_file_plot:
        parser.error('FLOW fitting requires --output_file_plot')
    if args.packet_summary and not (args.timestamp or args.input_file):
        parser.error('--packet-summary without FLOW requires --timestamp')
    if args.packet_file and args.packet_summary:
        parser.error('choose --packet-file or --packet-summary')
    results = []
    if args.input_file:
        try:
            fit_input = args.input_file
            if args.select_muons:
                with h5py.File(args.input_file, 'r') as source:
                    needs_selection = not any(path in source for _,path in (*SPLIT_SAMPLES, LEGACY_SAMPLE))
                if needs_selection:
                    from track_selection import select_muons
                    fit_input = str(Path(args.output_file_plot).with_suffix('.segments.h5'))
                    select_muons(args.input_file, fit_input)
            stamp = args.timestamp or date_from_filename(args.input_file).isoformat()
            fitted = fit_samples(fit_input, args.output_file_plot, args.sample, stamp)
            if fit_input != args.input_file:
                for row in fitted:
                    row['selection_file'] = os.path.abspath(fit_input)
                    row['input_file'] = os.path.abspath(args.input_file)
            results.extend(fitted)
        except Exception as exc:
            stamp = args.timestamp or date_from_filename(args.input_file).isoformat()
            results.append(dict(timestamp=aware_time(stamp).isoformat(), sample='all_mip', method='track',
                input_file=os.path.abspath(args.input_file), fit_status='unavailable',
                lifetime_us=None, error_us=None, fit_message=f'{type(exc).__name__}: {exc}'))
            print('Track sample unavailable:', exc)
        if args.write_hdf5_metadata:
            write_hdf5_metadata(args.input_file, [r for r in results if isinstance(r.get('segments_dset'), str) and 'selection_file' not in r])
    if args.packet_file or args.packet_summary:
        results.append(packet_result(args))
    if args.slow_controls:
        results.extend(aggregate_measurements(args.slow_controls, args.timezone, json.loads(Path(args.gas_quality_config).read_text()) if args.gas_quality_config else None))
    update_json(args.output_file_json, results, args.output_timeseries, args.timezone, args.annotations)
    if not results:
        print('Refreshed existing history; no new observations supplied.')
    elif not any(r.get('fit_status') == 'ok' for r in results):
        print('No new valid lifetimes. Unavailable/failed results recorded; history plot refreshed.')


if __name__ == '__main__':
    main()
