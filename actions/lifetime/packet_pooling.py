"""Fixed-window, detector-wide packet pooling with auditable acceptance gates.

Pool independent charge objects BEFORE drift-bin fitting. Failed per-file tau
estimates are never averaged. Every tried window, including failures, remains
inspectable. Cohorts and selection/calibration fingerprints must match.
"""
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from lifetime_io import output_lock as Lock
from scipy.optimize import least_squares

from purity_sources import atomic_json, aware_time

# Prospective feasibility gates, not detector-validated performance claims.
# Do not loosen them based on agreement with another lifetime method.
DEFAULT_GATES = dict(minimum_fit_bins=8, minimum_drift_span_us=100.,
                     maximum_attenuation_chi2_ndf=3., minimum_pedestal_coverage=.95)


def digest_json(value):
    """Stable configuration identity independent of JSON whitespace/order."""
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def store_sample(frame, metadata, directory):
    """Atomically publish compact, non-pickle objects even when a file fit fails.

    Array payloads are immutable/content-addressed; the small descriptor is
    replaced on a retry. Filename identity prevents double-counting copied
    versions of the same native file. Selection provenance stays in metadata.
    """
    directory = Path(directory); directory.mkdir(parents=True, exist_ok=True)
    dtype = [('q_sum', 'f8'), ('time_us', 'f8'), ('ext', 'i8'), ('io_group', 'i4')]
    objects = np.empty(len(frame), dtype=dtype)
    for field in objects.dtype.names:
        objects[field] = frame[field].to_numpy()
    content_hash = hashlib.sha256(objects.tobytes()).hexdigest()
    # The grouping/selection code version is explicit in the fingerprint.
    metadata = dict(metadata, object_sha256=content_hash, n_objects=len(objects),
                    sample_schema=1, created_at=datetime.now(timezone.utc).isoformat())
    identity = hashlib.sha256(Path(metadata['packet_file']).name.encode()).hexdigest()[:24]
    metadata['sample_id'] = identity
    metadata['configuration_id'] = digest_json(metadata['configuration'])
    array_path = directory/(identity+'-'+content_hash+'.npz')
    metadata['objects_file'] = array_path.name
    with Lock(str(directory/'.pool.lock'), default_timeout=timedelta(hours=6), lifetime=timedelta(hours=24)):
        if not array_path.exists():
            fd, temporary = tempfile.mkstemp(suffix='.npz', dir=directory); os.close(fd)
            try:
                np.savez_compressed(temporary, objects=objects)
                os.replace(temporary, array_path)
            finally:
                if os.path.exists(temporary): os.unlink(temporary)
        atomic_json(directory/(identity+'.sample.json'), metadata)
    return metadata


def window_start(timestamp, hours, timezone_name):
    """Fixed local calendar bins; never slide a boundary to rescue a fit."""
    if hours not in (6, 24, 48):
        raise ValueError('packet windows must be predeclared: 6, 24 or 48 hours')
    local = aware_time(timestamp).tz_convert(timezone_name).tz_localize(None)
    if hours == 6:
        start = local.normalize() + pd.Timedelta(hours=6*(local.hour//6))
    elif hours == 24:
        start = local.normalize()
    else:
        days = (local.normalize()-pd.Timestamp('1970-01-01')).days
        start = pd.Timestamp('1970-01-01') + pd.Timedelta(days=2*(days//2))
    return start.tz_localize(timezone_name)


def fit_objects(objects, configuration, gates=None, bootstrap=100):
    """Apply fixed slice criteria, then test the attenuation model.

    MPV bootstraps resample charge objects, as in the reference fitter. Final
    covariance treats drift-bin MPVs as independent; it excludes shared EXT,
    channel and calibration effects. Thus a passing result remains provisional.
    Keep every failed slice reason; more counts cannot fix a wrong peak model.
    """
    from packet_lifetime import fit_langau, bootstrap_mpv
    if bootstrap < 10:
        raise ValueError('at least 10 MPV bootstrap trials are required')
    gates = dict(DEFAULT_GATES, **(gates or {}))
    diagnostics = []; accepted = []
    width = configuration['time_bin_us']
    for i in range(int(configuration['drift_max_us']/width)):
        q = objects['q_sum'][(objects['time_us'] >= i*width) & (objects['time_us'] < (i+1)*width)]
        row = dict(time_bin=i, time_us=(i+.5)*width, n_objects=len(q), status='rejected')
        try:
            if len(q) < configuration['minimum_bin_objects']:
                raise ValueError('insufficient objects')
            fit = fit_langau(q)
            row.update(mpv=float(fit['mpv']) if np.isfinite(fit['mpv']) else None,
                       slice_chi2_ndf=float(fit['red_chi2']) if np.isfinite(fit['red_chi2']) else None)
            if row['mpv'] is None or row['mpv'] <= 0:
                raise ValueError('nonpositive or nonfinite slice MPV')
            if not np.isfinite(fit['red_chi2']) or fit['red_chi2'] >= configuration['maximum_slice_chi2_ndf']:
                raise ValueError('slice model disagreement')
            mpv, error, successful = bootstrap_mpv(q, bootstrap, 241000+i)
            row.update(mpv=float(mpv) if np.isfinite(mpv) else None, bootstrap_successes=successful,
                       mpv_error=float(error) if np.isfinite(error) else None)
            if (row['mpv'] is None or row['mpv'] <= 0 or successful < max(10, int(np.ceil(.8*bootstrap)))
                    or not np.isfinite(error) or error <= 0):
                raise ValueError('unstable MPV bootstrap')
            row['status'] = 'accepted'
            if configuration['fit_tmin_us'] <= row['time_us'] <= configuration['fit_tmax_us']:
                accepted.append(row)
        except Exception as exc:
            row['reason'] = str(exc)
        diagnostics.append(row)
    result = dict(lifetime_us=None, error_us=None, fit_status='unavailable',
                  slice_diagnostics=diagnostics, n_accepted_fit_bins=len(accepted), gates=gates,
                  uncertainty='MPV object-bootstrap errors propagated through independent-bin attenuation covariance; shared systematics/correlations excluded')
    reasons = []
    if len(accepted) < gates['minimum_fit_bins']:
        reasons.append('too few accepted drift bins')
    span = max((r['time_us'] for r in accepted), default=0)-min((r['time_us'] for r in accepted), default=0)
    result['accepted_drift_span_us'] = span
    if span < gates['minimum_drift_span_us']:
        reasons.append('insufficient accepted drift span')
    if reasons:
        result['fit_message'] = '; '.join(reasons)
        return result
    t = np.array([r['time_us'] for r in accepted])/1000
    q = np.array([r['mpv'] for r in accepted]); e = np.array([r['mpv_error'] for r in accepted])
    # Match Ahmad's inverse-lifetime bounds and exponential parameterization.
    def residual(parameters):
        alpha, amplitude = parameters
        return (q-amplitude*np.exp(-alpha*t))/e
    fit = least_squares(residual, [.5, float(np.median(q))], bounds=([-5, 0], [5, 50]), max_nfev=10000)
    if not fit.success or np.any(fit.active_mask) or np.linalg.matrix_rank(fit.jac) < 2:
        result['fit_message'] = 'attenuation failed, singular, or at parameter boundary'
        return result
    alpha = float(fit.x[0]); chi2 = float(np.sum(fit.fun**2)); ndf = len(q)-2
    covariance = np.linalg.pinv(fit.jac.T@fit.jac)
    error_alpha = float(np.sqrt(covariance[0,0]))
    result.update(alpha_per_ms=alpha, alpha_error_per_ms=error_alpha, attenuation_chi2=chi2,
                  attenuation_ndf=ndf, attenuation_chi2_ndf=chi2/ndf,
                  candidate_lifetime_us=1000/alpha if alpha > 0 else None)
    if chi2/ndf > gates['maximum_attenuation_chi2_ndf']:
        reasons.append('attenuation model disagreement')
    if alpha <= error_alpha:
        reasons.append('finite positive attenuation not resolved at one covariance sigma')
    if reasons:
        result['fit_message'] = '; '.join(reasons)
        return result
    result.update(lifetime_us=1000/alpha, error_us=1000*error_alpha/alpha**2, fit_status='ok',
                  fit_quality='provisional_packet_selection', covariance=covariance.tolist())
    return result


def pool_samples(directory, hours=24, timezone_name='America/Chicago', gates=None, bootstrap=100, snapshot=None):
    """Refit fixed windows of compatible file shards; cache unchanged membership.

    Read/fit/cache under a separate lock. Evaluation timestamps let the central
    history reject a late-arriving older worker's result. Source extraction is
    also serialized with this lock, so membership is coherent during a fit.
    """
    directory = Path(directory)
    if not directory.exists(): return []
    gates = dict(DEFAULT_GATES, **(gates or {}))
    # The reference slice fitter can change independently of this module.
    # Include both implementations so an old fit cannot survive a code update.
    fit_code_hash = hashlib.sha256(Path(__file__).read_bytes()+
                                  Path(__file__).with_name('packet_lifetime.py').read_bytes()).hexdigest()
    # A large fixed-window replay can outlive a one-hour lease. Keep the shared
    # lock longer than the documented six-hour worker limit, as monitor.py does.
    with Lock(str(directory/'.pool.lock'), default_timeout=timedelta(hours=6), lifetime=timedelta(hours=24)):
        groups = {}
        for path in sorted(directory.glob('*.sample.json')):
            row = json.loads(path.read_text())
            start = window_start(row['timestamp'], hours, timezone_name)
            group = (row['cohort'], row['configuration_id'], start)
            groups.setdefault(group, []).append(row)
        # Capture one coherent inventory while holding the extraction lock.
        # Publication can then retire obsolete configurations after a file is
        # re-extracted, without a delayed worker restoring the old inventory.
        snapshot_at = datetime.now(timezone.utc).isoformat()
        output = []
        for (cohort, config_id, start), members in sorted(groups.items()):
            end = (start.tz_localize(None)+pd.Timedelta(hours=hours)).tz_localize(timezone_name)
            fingerprint = digest_json(dict(members=members, hours=hours, gates=gates, bootstrap=bootstrap,
                                            fit_code_sha256=fit_code_hash))
            cache = directory/(fingerprint+'.pool.json')
            cached = cache.exists()
            arrays = []
            for member in members:
                with np.load(directory/member['objects_file'], allow_pickle=False) as data:
                    array = data['objects']
                if hashlib.sha256(array.tobytes()).hexdigest() != member['object_sha256']:
                    raise ValueError('packet object shard hash mismatch')
                # Cached numerical results are not evidence that their source
                # arrays still exist and match the extraction receipt. Verify
                # every member on reuse, but retain arrays in memory only when
                # the fit actually needs to be rebuilt.
                if not cached:
                    arrays.append(array)
            if cached:
                output.append(dict(json.loads(cache.read_text()), pool_snapshot_at=snapshot_at)); continue
            objects = np.concatenate(arrays)
            # A sparse pilot may have identical membership in 6 h and 24 h
            # bins. Reuse the same mathematical fit, while retaining separate
            # window records; nominal bin width does not create extra data.
            fit_fingerprint = digest_json(dict(objects=[m['object_sha256'] for m in members],
                configuration=config_id, gates=gates, bootstrap=bootstrap, code=fit_code_hash))
            fit_cache = directory/(fit_fingerprint+'.fit.json')
            if fit_cache.exists():
                result = json.loads(fit_cache.read_text())
            else:
                result = fit_objects(objects, members[0]['configuration'], gates, bootstrap)
                result['fit_code_sha256'] = fit_code_hash
                atomic_json(fit_cache,result)
            # Pedestal coverage is assessed across every input, not averaged in
            # a way that could hide one file's severely missing channels.
            coverage = min(r['pedestal_coverage'] for r in members)
            if coverage < gates['minimum_pedestal_coverage']:
                result.update(fit_status='unavailable', lifetime_us=None, error_us=None,
                    fit_message='; '.join(filter(None, [result.get('fit_message'), 'insufficient pedestal coverage'])))
            result.update(timestamp=(start+(end-start)/2).isoformat(), period_start=start.isoformat(), period_end=end.isoformat(),
                source=cohort+' / '+config_id[:12]+' / '+str(hours)+' h', sample='packet', method='packet_pool',
                aggregation='detector-wide charge objects pooled before drift-bin fits', window_hours=hours,
                configuration=members[0]['configuration'], n_files=len(members), n_objects=len(objects),
                minimum_pedestal_coverage=coverage, cohort=cohort,
                first_observed_at=min(r['timestamp'] for r in members),
                last_observed_at=max(r['timestamp'] for r in members),
                summed_file_duration_s=sum(r.get('duration_s',0) for r in members),
                members=[{k:r.get(k) for k in ('sample_id','packet_file','flow_file','timestamp','n_objects','pedestal_coverage','object_sha256','duration_s')} for r in members],
                evaluated_at=datetime.now(timezone.utc).isoformat(),
                note='Nominal bin width is not live exposure. File-start assignment; unknown within-window stability. Passing fit remains provisional.')
            atomic_json(cache, result)
            output.append(dict(result, pool_snapshot_at=snapshot_at))
        if snapshot is not None:
            snapshot.update(evaluated_at=snapshot_at, window_hours=hours,
                            keys=[[r['method'],r['source'],r['period_start']] for r in output])
    return output


def main():
    """Report all prespecified candidate windows; do not select by fitted tau."""
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', required=True)
    parser.add_argument('--hours', type=int, nargs='+', default=[6,24,48])
    parser.add_argument('--output', required=True)
    parser.add_argument('--bootstrap', type=int, default=100)
    args = parser.parse_args()
    results = [r for h in args.hours for r in pool_samples(args.samples,h,bootstrap=args.bootstrap)]
    atomic_json(args.output, dict(windows=results, candidate_hours=args.hours,
        conclusion='Feasibility comparison only; no automatic choice of averaging period or physics validation.'))


if __name__ == '__main__':
    main()
