#!/usr/bin/env python3
"""Inspect zero-drift charge scales without changing the lifetime estimator.

The saved lifetime amplitude is a LEGACY PROXY, not raw dQ/ds. The optional
direct diagnostic fits Langau MPVs of dQ/ds from the full pooled segment array.
These two observables have different units and must never share a fitted curve.
"""
import argparse
from datetime import datetime, timezone
import glob
import hashlib
import html
import json
from pathlib import Path
import warnings

import h5py
import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit

from lifetime_funcs import ElectronLifetimeFunc, langau_func, attenuation_covariance_diagnostics


DIRECT_VERSION = 'direct-dqds-langau-mpv-exponential-v1'
TIME_BINS_US = np.linspace(0., 196., 20)
CHARGE_BINS = np.linspace(0., 250., 70)
NOMINAL_BINS = tuple(range(4, 18))
MIN_USABLE_BINS = 10
MAX_TAU_US = 20_000.


def _finite(value):
    try:
        return float(value) if np.isfinite(float(value)) else None
    except (TypeError, ValueError, OverflowError):
        return None


def _empty(observable, units):
    return dict(status='unavailable', observable=observable, model='A exp(-t / tau)',
        units=units, intercept=None, intercept_error=None, tau_us=None, tau_error_us=None,
        covariance=None, correlation=None, fit_message=None, fit_range_us=None,
        fitpoints=dict(x=[], y=[], error=[], bins=[]), curve=dict(x=[], y=[], error=[]),
        uncertainty='Conditional fit covariance only; excludes calibration, selection and within-track correlations')


def _covariance_ok(covariance):
    """Require a finite positive-definite correlation matrix for all parameters."""
    cov = np.asarray(covariance, dtype=float)
    if cov.ndim != 2 or cov.shape[0] != cov.shape[1] or not np.isfinite(cov).all():
        return False
    diagonal = np.diag(cov)
    if np.any(diagonal <= 0):
        return False
    scale = np.sqrt(diagonal)
    corr = cov / scale[:, None] / scale[None, :]
    if not np.isfinite(corr).all():
        return False
    tolerance = np.finfo(float).eps * len(cov) * max(1., np.linalg.norm(corr, 2))
    return bool(np.allclose(corr, corr.T, rtol=0, atol=tolerance) and
                np.linalg.eigvalsh((corr + corr.T) / 2)[0] > tolerance)


def _accepted(result, amplitude, tau, covariance):
    covariance = np.asarray(covariance, dtype=float)
    amplitude, tau = float(amplitude), float(tau)
    result.update(status='ok', intercept=amplitude, intercept_error=float(np.sqrt(covariance[0, 0])),
        tau_us=tau, tau_error_us=float(np.sqrt(covariance[1, 1])), covariance=covariance.tolist(),
        correlation=float(covariance[0, 1] / np.sqrt(covariance[0, 0]) / np.sqrt(covariance[1, 1])))
    if result['fitpoints']['x']:
        first, last = min(result['fitpoints']['x']), max(result['fitpoints']['x'])
        result['fit_range_us'] = [first, last]
        times = np.linspace(0., last, 240)
        exponential = np.exp(-times / tau)
        curve = amplitude * exponential
        jacobian = np.column_stack((exponential, curve * times / tau**2))
        variance = np.einsum('ij,jk,ik->i', jacobian, covariance, jacobian)
        result['curve'] = dict(x=times.tolist(), y=curve.tolist(), error=np.sqrt(np.maximum(variance, 0.)).tolist())
    return result


def charge_scale_from_summary(summary):
    """Extract the existing accepted lifetime-fit proxy intercept; never refit."""
    result = _empty('legacy lifetime proxy MPV', 'legacy proxy units')
    result['method'] = 'saved lifetime fit; no refit'
    fit = summary.get('fit') or summary
    diagnostic = fit.get('fit_diagnostics') or summary.get('fit_diagnostics') or {}
    result['estimator_version'] = diagnostic.get('estimator_version')
    times, values, errors = (diagnostic.get(key, []) for key in ('drift_time_us', 'mpv', 'mpv_error'))
    for index in diagnostic.get('fitted_bins', []):
        if not isinstance(index, int) or index < 0 or index >= min(len(times), len(values), len(errors)):
            continue
        x, y, error = (_finite(array[index]) for array in (times, values, errors))
        if x is not None and y is not None and error is not None and error > 0:
            for key, value in [('x', x), ('y', y), ('error', error), ('bins', index)]:
                result['fitpoints'][key].append(value)
    result['chi2_ndf'] = _finite(diagnostic.get('chi2_ndf'))
    result['ndf'] = diagnostic.get('ndf')
    if (summary.get('errors') or summary.get('fit_status', fit.get('fit_status')) != 'ok'
            or fit.get('fit_status', 'ok') != 'ok'):
        result['fit_message'] = summary.get('fit_message') or fit.get('fit_message') or 'Original lifetime fit was not accepted; no accepted proxy intercept.'
        return result
    parameters = diagnostic.get('candidate_parameters') or {}
    amplitude, tau = _finite(parameters.get('amplitude')), _finite(parameters.get('lifetime_us'))
    covariance = diagnostic.get('attenuation_covariance')
    if amplitude is None or tau is None or amplitude <= 0 or tau <= 0 or tau > MAX_TAU_US:
        result['fit_message'] = 'Saved fit lacks physical accepted amplitude/lifetime parameters.'
        return result
    try:
        check = attenuation_covariance_diagnostics(covariance)
    except (ValueError, TypeError):
        check = dict(valid=False, reason='invalid covariance representation')
    if not check['valid']:
        result['fit_message'] = 'Saved amplitude covariance is unavailable or invalid: ' + str(check.get('reason'))
        return result
    return _accepted(result, amplitude, tau, covariance)


def _fit_direct_peak(hist, centers):
    """Unweighted histogram least squares, with a resolved four-parameter fit.

    The fitting window is fixed relative to the histogram maximum (0.5--1.5
    times its coordinate), matching the historical peak-window convention.
    This is a conditional diagnostic, not a Poisson likelihood measurement.
    """
    maximum = float(centers[np.argmax(hist)])
    selected = (centers >= .5 * maximum) & (centers <= 1.5 * maximum)
    x, y = centers[selected], hist[selected]
    if len(x) < 5 or np.count_nonzero(y) < 5:
        raise ValueError('fewer than five populated peak-window centers for four Langau parameters')
    std = float(np.std(x))
    lower = np.array([0., .01, 0., 0.])
    initial = [maximum, max(.02, std / 6), std / 4, float(np.max(y))]
    with warnings.catch_warnings():
        warnings.simplefilter('error', OptimizeWarning)
        # pylandau locates its maximum numerically. Machine-scale differences
        # can differentiate that numerical noise instead of the peak shape.
        params, covariance = curve_fit(langau_func, x, y, p0=initial,
            bounds=(lower, np.full(4, np.inf)), max_nfev=10000,
            diff_step=1e-4, x_scale='jac', ftol=1e-10, xtol=1e-10, gtol=1e-10)
    if not np.isfinite(params).all() or np.any(params <= lower + 1e-8 * np.maximum(1., np.abs(params))):
        raise ValueError('Langau parameter lies on a nonnegative/width boundary')
    if params[2] > 100 * params[1]:
        raise ValueError('Langau sigma exceeds 100 times eta, where pylandau warns of numerical oscillations')
    if not _covariance_ok(covariance):
        raise ValueError('unresolved four-parameter Langau covariance')
    return float(params[0]), float(np.sqrt(covariance[0, 0])), len(x)


def fit_direct_charge_scale(segments, tick_us=.1):
    """Fit full-pool dQ/ds MPVs versus drift; independent of the proxy estimate."""
    from recent_tracks import fit_warning_receipt
    with fit_warning_receipt() as fit_warnings:
        result = _fit_direct_charge_scale(segments, tick_us)
    result['fit_warnings'] = fit_warnings
    return result


def _fit_direct_charge_scale(segments, tick_us):
    result = _empty('direct dQ/ds Langau MPV', 'FLOW charge units / cm')
    result.update(estimator_version=DIRECT_VERSION, method='full-pool direct dQ/ds diagnostic; separate from the lifetime proxy',
        nominal_bins=list(NOMINAL_BINS), time_edges_us=TIME_BINS_US.tolist(), dqds_edges=CHARGE_BINS.tolist(),
        minimum_usable_bins=MIN_USABLE_BINS, minimum_segments_per_bin=5, minimum_populated_peak_centers=5,
        uncertainty='Conditional MPV/attenuation covariance; unweighted histogram slice fits; excludes calibration, selection and within-track correlations',
        bin_status=[], n_segments=len(segments))
    if not np.isfinite(tick_us) or tick_us <= 0:
        raise ValueError('tick_us must be finite and positive')
    required = {'dx', 'dQ', 'nhits', 't'}
    if not required <= set(segments.dtype.names or ()):
        raise ValueError('direct charge-scale fit requires dx, dQ, nhits and t')
    valid = np.logical_and.reduce([np.isfinite(segments[name]) for name in required] +
                                 [segments[name] > 0 for name in ('dx', 'dQ', 'nhits')])
    result['n_valid_segments'] = int(valid.sum())
    charge = segments['dQ'][valid] / segments['dx'][valid]
    time = segments['t'][valid] * tick_us
    centers = (CHARGE_BINS[:-1] + CHARGE_BINS[1:]) / 2
    for index, (start, stop) in enumerate(zip(TIME_BINS_US[:-1], TIME_BINS_US[1:])):
        mask = (time >= start) & ((time < stop) if index < len(TIME_BINS_US) - 2 else (time <= stop))
        values = charge[mask]
        record = dict(bin=index, n_segments=int(mask.sum()), status='excluded', reason='outside nominal attenuation range')
        result['bin_status'].append(record)
        if index not in NOMINAL_BINS:
            continue
        in_range = (values >= CHARGE_BINS[0]) & (values <= CHARGE_BINS[-1])
        record['n_histogram_segments'] = int(in_range.sum())
        if in_range.sum() < 5:
            record.update(status='insufficient_support', reason='fewer than five in-range segments')
            continue
        hist, _ = np.histogram(values, bins=CHARGE_BINS, density=True)
        try:
            value, error, n_centers = _fit_direct_peak(hist, centers)
            if not np.isfinite([value, error]).all() or value <= 0 or error <= 0:
                raise ValueError('nonpositive or nonfinite direct MPV/error')
        except (ValueError, RuntimeError, FloatingPointError, OptimizeWarning) as exc:
            record.update(status='peak_fit_failed', reason=f'{type(exc).__name__}: {exc}')
            continue
        record.update(status='usable', reason=None, peak_centers=n_centers)
        for key, number in [('x', float((start + stop) / 2)), ('y', value), ('error', error), ('bins', index)]:
            result['fitpoints'][key].append(number)
    result['usable_bins'] = len(result['fitpoints']['x'])
    if result['usable_bins'] < MIN_USABLE_BINS:
        result.update(status='failed', fit_message=f"Insufficient usable direct dQ/ds bins: {result['usable_bins']}; require {MIN_USABLE_BINS}.")
        return result
    x, y, error = (np.asarray(result['fitpoints'][key]) for key in ('x', 'y', 'error'))
    try:
        params, covariance = curve_fit(ElectronLifetimeFunc, x, y, sigma=error, absolute_sigma=True,
            p0=[float(y.max()), 2000.], bounds=([0., 0.], [np.inf, np.inf]), max_nfev=10000)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        result.update(status='failed', fit_message=f'Direct attenuation fit failed: {type(exc).__name__}: {exc}')
        return result
    check = attenuation_covariance_diagnostics(covariance)
    result['candidate_parameters'] = dict(intercept=_finite(params[0]), tau_us=_finite(params[1]))
    result['covariance_check'] = check
    if not np.isfinite(params).all() or np.any(params <= 1e-8):
        result.update(status='failed', fit_message='Direct attenuation parameter is nonfinite or at its positive boundary.')
        return result
    if not check['valid']:
        result.update(status='failed', fit_message='Direct attenuation covariance invalid: ' + check['reason'])
        return result
    if params[1] > MAX_TAU_US:
        result.update(status='rejected', fit_message=f'Direct diagnostic lifetime exceeds {MAX_TAU_US:g} us acceptance ceiling.')
        return result
    result['chi2'] = float(np.sum(((y - ElectronLifetimeFunc(x, *params)) / error)**2))
    result['ndf'] = len(x) - 2
    result['chi2_ndf'] = result['chi2'] / result['ndf']
    return _accepted(result, *params, covariance)


def charge_scale_payload(summary, segments=None, tick_us=.1):
    direct = (_empty('direct dQ/ds Langau MPV', 'FLOW charge units / cm') if segments is None
              else fit_direct_charge_scale(segments, tick_us))
    if segments is None:
        direct['fit_message'] = 'Direct diagnostic not requested; requires the full pooled segment array.'
    return dict(proxy=charge_scale_from_summary(summary), direct=direct)


def _time(value):
    result = datetime.fromisoformat(str(value).replace('Z', '+00:00'))
    if result.tzinfo is None:
        raise ValueError('window timestamps require an explicit UTC offset')
    return result.astimezone(timezone.utc)


def _json(value):
    if isinstance(value, dict): return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_json(v) for v in value]
    if isinstance(value, (np.ndarray, np.generic)): return _json(value.tolist())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def monitor_records(paths, direct=False):
    """Read every receipt, preserving failures and deduplicating copied receipts."""
    records, seen = [], {}
    for path in paths:
        path = Path(path).resolve()
        try:
            summary = json.loads(path.read_text())
            if not isinstance(summary, dict) or 'members' not in summary or 'selection' not in summary:
                raise ValueError('expected a recent_tracks summary receipt')
            start, stop = _time(summary['period_start']), _time(summary['period_end'])
            if stop <= start:
                raise ValueError('window end must follow its start')
            cohorts = sorted({Path(member['input_file']).parent.name for member in summary.get('members', [])})
            cohort_hash = hashlib.sha256(json.dumps(dict(runs=cohorts,
                reconstruction=summary.get('reconstruction_signature'),
                proxy_estimator=(summary.get('fit') or {}).get('fit_diagnostics', {}).get('estimator_version')),
                sort_keys=True).encode()).hexdigest()[:12]
            result = charge_scale_payload(summary)
            if direct and summary.get('errors'):
                result['direct']['fit_message'] = 'Full-pool direct diagnostic omitted because the input receipt contains errors.'
            if direct and not summary.get('errors') and summary.get('pooled_segments'):
                try:
                    pooled = Path(summary['pooled_segments'])
                    if not pooled.is_absolute(): pooled = path.parent / pooled
                    with h5py.File(pooled, 'r') as source:
                        if not bool(source.attrs.get('complete', False)):
                            raise ValueError('incomplete pooled segment file')
                        key = 'analysis/rock_muon_segments/data'
                        data = source[key][:]
                        if len(data) != summary['n_valid_segments']:
                            raise ValueError('pooled segment count disagrees with receipt')
                        hashes = [entry.get('sha256') for entry in (summary.get('fit') or {}).get('flow_datasets', []) if entry.get('path') == key]
                        if len(hashes) != 1 or hashes[0] != hashlib.sha256(data.tobytes()).hexdigest():
                            raise ValueError('pooled segment SHA-256 is missing or disagrees with fitted input receipt')
                        signature = summary.get('reconstruction_signature') or {}
                        ticks = (signature.get('reconstruction') or {}).get('run_info') or {}
                        tick_us = ticks.get('crs_ticks')
                        if tick_us is None:
                            raise ValueError('missing recorded pooled drift tick period')
                        result['direct'] = fit_direct_charge_scale(data, float(tick_us))
                except (OSError, ValueError, KeyError, TypeError) as exc:
                    result['direct'].update(input_error=True,
                        fit_message=f'Direct diagnostic input error: {type(exc).__name__}: {exc}')
            pool_identity = summary.get('pooled_segments')
            if pool_identity:
                pool_identity = Path(pool_identity)
                pool_identity = str((pool_identity if pool_identity.is_absolute() else path.parent / pool_identity).resolve())
            record = dict(receipt=str(path), period_start=start.isoformat(), period_end=stop.isoformat(),
                timestamp=(start + (stop - start) / 2).isoformat(), selection=summary['selection'],
                cohort=cohort_hash, runs=cohorts, n_segments=summary.get('n_valid_segments'),
                pooled_segments=pool_identity, **result)
            identity = (record['pooled_segments'] or str(path), record['period_start'], record['period_end'], record['selection'])
            comparison = json.dumps({key: value for key, value in record.items() if key != 'receipt'}, sort_keys=True)
            if identity in seen:
                if seen[identity] != comparison:
                    raise ValueError('conflicting duplicate receipt for the same pooled sample')
                continue
            seen[identity] = comparison
            records.append(record)
        except (OSError, ValueError, KeyError, TypeError) as exc:
            records.append(dict(receipt=str(path), error=f'{type(exc).__name__}: {exc}'))
    return records


def write_monitor(records, output, json_output):
    from plotly.offline import get_plotlyjs
    traces = []
    groups = sorted({(r['selection'], r['cohort'], kind) for r in records if 'error' not in r for kind in ('proxy', 'direct')})
    colors = {'stored': '#596c87', 'quality-stored': '#087f76', 'reselect': '#a56b15'}
    for selection, cohort, kind in groups:
        rows = sorted([r for r in records if r.get('selection') == selection and r.get('cohort') == cohort], key=lambda r: r['timestamp'])
        if not any(r[kind]['status'] == 'ok' for r in rows): continue
        traces.append(dict(type='scatter', mode='lines+markers', x=[r['timestamp'] for r in rows],
            y=[r[kind]['intercept'] if r[kind]['status'] == 'ok' else None for r in rows],
            error_y=dict(type='data', array=[r[kind]['intercept_error'] if r[kind]['status'] == 'ok' else None for r in rows], visible=True),
            text=['<br>'.join([html.escape(' + '.join(r['runs']) or 'unknown acquisition cohort'), html.escape(r['period_start'] + ' → ' + r['period_end']), 'tau='+str(r[kind]['tau_us'])+' µs', 'χ²/ndf='+str(r[kind].get('chi2_ndf'))]) for r in rows],
            hovertemplate='%{x}<br>A=%{y}<br>%{text}<extra>%{fullData.name}</extra>',
            name=f'{selection} · {cohort} · {kind}', connectgaps=False,
            line=dict(color=colors.get(selection, '#596c87'), dash='solid' if kind == 'proxy' else 'dot'),
            xaxis='x' if kind == 'proxy' else 'x2', yaxis='y' if kind == 'proxy' else 'y2'))
    accepted = {kind: any(row.get(kind, {}).get('status') == 'ok' for row in records)
                for kind in ('proxy', 'direct')}
    annotations = [dict(text=f'No accepted {label} fits<br>See receipt diagnostics below.',
        x=.5, y=y, xref='paper', yref='paper', showarrow=False, font=dict(size=14, color='#596c7e'))
        for kind, label, y in [('proxy', 'proxy', .79), ('direct', 'direct dQ/ds', .21)]
        if not accepted[kind]]
    layout = dict(height=750 if any(accepted.values()) else 400, margin=dict(l=80,r=30,t=45,b=60), hovermode='closest',
        xaxis=dict(title='Window midpoint [UTC]', type='date', anchor='y', domain=[0,1], visible=accepted['proxy']),
        yaxis=dict(title='Saved proxy intercept<br>[proxy units]', anchor='x', domain=[.58,1], visible=accepted['proxy']),
        xaxis2=dict(title='Window midpoint [UTC]', type='date', anchor='y2', domain=[0,1], visible=accepted['direct']),
        yaxis2=dict(title='Direct dQ/ds intercept<br>[FLOW charge units / cm]', anchor='x2', domain=[0,.42], visible=accepted['direct']),
        annotations=annotations,
        legend=dict(orientation='h', y=1.12), paper_bgcolor='white', plot_bgcolor='white')
    body = []
    for row in records:
        if 'error' in row:
            body.append('<tr><td colspan="6">'+html.escape(row['receipt']+': '+row['error'])+'</td></tr>')
        else:
            for kind in ('proxy','direct'):
                value=row[kind]
                body.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in
                    (row['timestamp'],row['selection'],row['cohort'],kind,value['status'],
                     f"{value['intercept']} ± {value['intercept_error']}" if value['status']=='ok' else value.get('fit_message')))+'</tr>')
    windows = [row for row in records if 'error' not in row]
    filename = 'charge_scale_windows_no_valid_receipts'
    if windows:
        start = min(_time(row['period_start']) for row in windows)
        stop = max(_time(row['period_end']) for row in windows)
        filename = f"charge_scale_{start:%Y%m%dT%H%M%SZ}_to_{stop:%Y%m%dT%H%M%SZ}"
    config = dict(responsive=True, displaylogo=False,
                  toImageButtonOptions=dict(format='png', filename=filename, scale=2))
    encoded=json.dumps(dict(traces=traces,layout=layout,config=config),allow_nan=False).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    document='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Zero-drift charge-scale monitor</title>
<style>body{font:14px system-ui;margin:28px;color:#172b40}h1{margin-bottom:8px}p{max-width:1000px}.receipt-table{overflow-x:auto}table{border-collapse:collapse;width:100%;min-width:720px;font-size:12px}td,th{text-align:left;border-bottom:1px solid #ddd;padding:8px;overflow-wrap:anywhere}#plot{width:100%}@media(max-width:600px){body{margin:16px}h1{font-size:25px}}</style>
<h1>Zero-drift charge-scale monitor</h1><p>Saved lifetime-proxy intercept and optional direct dQ/ds MPV intercept are separate observables. Error bars are conditional covariance uncertainties. Curves join only matching recorded acquisition/configuration cohorts; failed fits retain null values. Equal metadata does not establish stable live detector conditions.</p>
<div id="plot"></div><div class="receipt-table" role="region" aria-label="Receipt diagnostics; scroll horizontally on small screens" tabindex="0"><table><thead><tr><th>Midpoint UTC</th><th>Selection</th><th>Cohort</th><th>Observable</th><th>Status</th><th>Intercept / diagnostic</th></tr></thead><tbody>'''+''.join(body)+'''</tbody></table></div><script>'''+get_plotlyjs()+'''</script><script>const monitor='''+encoded+''';Plotly.newPlot('plot',monitor.traces,monitor.layout,monitor.config);</script></html>'''
    output, json_output = Path(output), Path(json_output)
    output.parent.mkdir(parents=True,exist_ok=True)
    json_output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(document)
    json_output.write_text(json.dumps(_json(dict(schema_version=1,records=records)),indent=2,allow_nan=False)+'\n')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summaries',nargs='+',help='summary files or quoted globs; copied receipts are deduplicated')
    parser.add_argument('--output',type=Path,required=True,help='offline HTML monitor')
    parser.add_argument('--json-output',type=Path,help='default: output stem with .json suffix')
    parser.add_argument('--direct',action='store_true',help='also fit direct dQ/ds MPVs from each full pooled HDF5')
    args=parser.parse_args(argv)
    try:
        paths=[]
        for specification in args.summaries:
            matched=glob.glob(specification,recursive=True)
            if not matched: raise ValueError(f'no receipts match {specification}')
            paths.extend(Path(p).resolve() for p in matched)
        paths=sorted(set(paths))
        json_output=args.json_output or args.output.with_suffix('.json')
        if args.output.suffix.lower()!='.html' or json_output.suffix.lower()!='.json':
            raise ValueError('output requires .html and json-output requires .json')
        protected=set(paths)
        for path in paths:
            try:
                summary=json.loads(path.read_text())
                referenced=[summary.get('pooled_segments')]
                for member in summary.get('members',[]):referenced.extend([member.get('input_file'),member.get('selection_file')])
                for reference in referenced:
                    if reference:
                        candidate=Path(reference)
                        protected.add((candidate if candidate.is_absolute() else path.parent/candidate).resolve())
            except (OSError,ValueError,AttributeError,TypeError):
                pass  # The per-receipt error is retained in monitor_records.
        if args.output.resolve()==json_output.resolve() or any(p in protected for p in (args.output.resolve(),json_output.resolve())):
            raise ValueError('outputs must be distinct from every input summary, FLOW, sidecar and pool file')
        records=monitor_records(paths,args.direct)
        write_monitor(records,args.output,json_output)
    except (OSError,ValueError,TypeError) as exc:
        parser.exit(2,f'Charge-scale monitor failed: {exc}\n')
    good=sum(r.get(kind,{}).get('status')=='ok' for r in records for kind in ('proxy','direct'))
    print(f'{len(records)} receipts, {good} accepted intercepts: {args.output}')
    return 0 if good and not any('error' in r or r.get('direct',{}).get('input_error') for r in records) else 2


if __name__=='__main__':
    raise SystemExit(main())
