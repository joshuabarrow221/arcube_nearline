#!/usr/bin/env python3
import argparse
import json
import os
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from flufl.lock import Lock

import lifetime_funcs as LT
from nearline_util import date_from_filename


SPLIT_SAMPLES = (
    ('beam', 'analysis/beam_rock_muon_segments/data'),
    ('cosmic', 'analysis/cosmic_muon_segments/data'),
)
LEGACY_SAMPLE = ('mixed', 'analysis/rock_muon_segments/data')


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
        return split_samples
    if LEGACY_SAMPLE[1] in h5_file:
        return [LEGACY_SAMPLE]
    raise KeyError('no beam/cosmic or legacy rock-muon segment dataset found')


def valid_segments(segments):
    """Remove unusable segment rows before histogramming and fitting."""
    mask = segments['dx'] > 0
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


def fit_samples(input_file, output_file_plot, requested_samples=None):
    timestamp = date_from_filename(input_file)
    results = []

    print(f'Opening file: {input_file}')
    with h5py.File(input_file, 'r') as h5_file:
        samples = find_samples(h5_file, requested_samples)
        multiple_samples = len(samples) > 1

        for sample, dset_path in samples:
            segments = valid_segments(h5_file[dset_path][:])
            print(f'Extracting {sample} lifetime using {len(segments)} segments')

            result = {
                'timestamp': timestamp.isoformat(),
                'sample': sample,
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


def update_json(output_file_json, results):
    """Update the shared time-series store while making retries idempotent."""
    with Lock(output_file_json + '.lock'):
        if os.path.exists(output_file_json):
            with open(output_file_json) as stream:
                data = json.load(stream)
        else:
            data = {'lifetimes': []}

        entries = data.setdefault('lifetimes', [])
        new_keys = {(entry['input_file'], entry['sample']) for entry in results}
        entries[:] = [
            entry for entry in entries
            if (entry.get('input_file'), entry.get('sample', 'mixed')) not in new_keys
        ]
        entries.extend(results)

        with open(output_file_json, 'w') as stream:
            json.dump(data, stream, indent=2)


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


def main(
    input_file,
    output_file_plot,
    output_file_json,
    sample=None,
    write_hdf5_metadata_enabled=False,
):
    results = fit_samples(input_file, output_file_plot, sample)
    update_json(output_file_json, results)
    if write_hdf5_metadata_enabled:
        write_hdf5_metadata(input_file, results)

    if not any(result['fit_status'] == 'ok' for result in results):
        raise RuntimeError('all requested lifetime fits failed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input_file',
        required=True,
        help='Path to a flowed charge HDF5 file',
    )
    parser.add_argument(
        '--output_file_plot',
        required=True,
        help='Output plot path (sample labels are added when more than one sample is fit)',
    )
    parser.add_argument(
        '--output_file_json',
        required=True,
        help='Shared electron-lifetime JSON time-series file',
    )
    parser.add_argument(
        '--sample',
        action='append',
        type=parse_sample_spec,
        help=(
            'Optional LABEL=HDF5_PATH definition; repeat for multiple samples. '
            'By default, split beam/cosmic datasets are auto-detected, with fallback '
            'to analysis/rock_muon_segments/data.'
        ),
    )
    parser.add_argument(
        '--write-hdf5-metadata',
        dest='write_hdf5_metadata_enabled',
        action='store_true',
        help='Write each fit result as attributes on its segment dataset group',
    )

    args = parser.parse_args()
    print(args)
    main(**vars(args))
