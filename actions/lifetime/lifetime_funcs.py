"""DeMario's historical segment-based lifetime estimator.

The histogram convolution is a legacy proxy observable, not a newly calibrated
charge distribution. Preserve its binning and fit choices for comparison with
past monitoring. Covariance errors describe the fitted model only; they omit
shared calibration, track-selection and within-track correlations.
"""
import pylandau
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import warnings


# Record the estimator behavior separately from input reconstruction metadata.
# Pooled-result caches already hash this file, so changing this implementation
# invalidates old fits without changing the identity of measured FLOW arrays.
TRACK_ESTIMATOR_VERSION = 'demario-proxy-active-bins-v3-min10'
MIN_USABLE_DRIFT_BINS = 10


class TrackFitError(ValueError):
    """Retain pre-fit/fit evidence without publishing an accepted lifetime."""

    def __init__(self, message, diagnostics, fit_quality='fit_failed'):
        super().__init__(message)
        self.fit_diagnostics = diagnostics
        self.fit_quality = fit_quality


class InvalidAttenuationCovariance(TrackFitError):
    """Keep the numerical attempt inspectable without publishing a lifetime."""

    def __init__(self, message, diagnostics):
        super().__init__(message, diagnostics, 'invalid_attenuation_covariance')


def attenuation_covariance_diagnostics(covariance):
    """Check whether the two-parameter covariance represents a resolved fit.

    Bounded scipy curve_fit uses a pseudoinverse for a rank-deficient Jacobian.
    For a near-constant exponential this can return a tiny positive tau variance
    even when tau is unconstrained, rather than a useful uncertainty. Test the
    dimensionless correlation matrix: rank tests on the raw covariance depend
    on the incomparable amplitude/tau units. The threshold below is floating-
    point numerical rank, not a lifetime range or statistical-significance cut.
    """
    cov = np.asarray(covariance, dtype=float)
    result = dict(valid=False, reason=None, correlation=None,
                  eigenvalues=None, rank_tolerance=None)
    if cov.shape != (2, 2):
        result['reason'] = 'expected a 2 by 2 attenuation covariance'
        return result
    if not np.all(np.isfinite(cov)):
        result['reason'] = 'nonfinite attenuation covariance'
        return result
    variances = np.diag(cov)
    if np.any(variances <= 0):
        result['reason'] = 'nonpositive attenuation parameter variance'
        return result
    scale = np.sqrt(variances)
    # Divide sequentially to avoid overflow/underflow in an outer product of
    # very different parameter scales. Nonfinite normalization is a failure.
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        correlation = cov / scale[:, None] / scale[None, :]
    if not np.all(np.isfinite(correlation)):
        result['reason'] = 'nonfinite normalized attenuation covariance'
        return result
    result['correlation'] = correlation.tolist()
    norm = float(np.linalg.norm(correlation, 2))
    if not np.isfinite(norm):
        result['reason'] = 'nonfinite normalized attenuation covariance norm'
        return result
    tolerance = np.finfo(float).eps * len(cov) * max(1., norm)
    result['rank_tolerance'] = float(tolerance)
    if not np.allclose(correlation, correlation.T, rtol=0, atol=tolerance):
        result['reason'] = 'asymmetric normalized attenuation covariance'
        return result
    eigenvalues = np.linalg.eigvalsh(correlation / 2 + correlation.T / 2)
    if not np.all(np.isfinite(eigenvalues)):
        result['reason'] = 'nonfinite normalized attenuation covariance eigenvalues'
        return result
    result['eigenvalues'] = eigenvalues.tolist()
    if eigenvalues[0] <= tolerance:
        result['reason'] = 'singular or non-positive-definite normalized attenuation covariance'
        return result
    result['valid'] = True
    return result


def finite_diagnostic(value):
    """A failed numerical diagnostic must still be representable in strict JSON."""
    return float(value) if np.isfinite(value) else None

def convolution(hist_a, edges_a, hist_b, edges_b):
    """Preserve the legacy discrete-bin convolution and first-grid axis.

    The caller's hit-density and charge-density bin widths differ. This is the
    historical proxy construction; do not interpret its horizontal axis as a
    physical sum of two observables or silently replace it during maintenance.
    """
    
    conv = np.convolve(hist_a, hist_b)

    
    bin_width = edges_a[1] - edges_a[0]

    
    new_edges = np.arange(
        edges_a[0] + edges_b[0],
        edges_a[0] + edges_b[0] + (len(conv) * bin_width) + 1e-9,
        bin_width
    )
    new_edges = new_edges[:len(conv)+1]

    return conv, new_edges

def ElectronLifetimeFunc(x, A, tau):
    "Electron Lifetime function."
    return A*np.exp(-x/tau)

def chi_squared1(params, function, x, y, uncertainties):
    """Computes Chi Squared."""
    y_pred = function(x, *params)

    top = y-y_pred
    nom = top**2
    dom = uncertainties**2
    chi_array = nom/dom
    chi2 = np.sum(chi_array)

    return chi2

def langau_func(x, mpv, eta, sigma, A):
    """Landau distribution function."""
    return pylandau.langau(x, mpv, eta, sigma, A)
def langau_fit(hist, bin_centers):
    """Fit a slice's proxy peak; return the MPV covariance error.

    This historical least-squares fit supplies no per-bin counting errors;
    scipy therefore scales its covariance by the fitted residual variance.
    It is not a Poisson likelihood or a bootstrap uncertainty estimate.
    """
    # Ignore specific warnings
    warnings.filterwarnings("ignore", message="sigma > 100 * eta can lead to oszillations. Check result.")

    std_bin_centers = np.std(bin_centers)
    eta_guess, sigma_guess = std_bin_centers/6, std_bin_centers/4
    
    guess_mpv = np.average(bin_centers, weights = hist)
    mpv_guess, A_guess = bin_centers[np.argmax(hist)], np.max(hist)

    # Initial guess for the parameters: MPV, eta (Landau width), sigma (Gaussian width), A (amplitude)
    initial_guess = [mpv_guess, eta_guess,sigma_guess, A_guess] # Example values

    # Set lower and upper bounds for sigma, and other parameters
    bounds = ([0, .01, 0, 0], [np.inf, np.inf,np.inf, np.inf])  # Nonnegative parameters; no finite upper bound.
    
    try:
        # Perform the curve fit with bounds
        params, cov = curve_fit(langau_func, bin_centers, hist, p0=initial_guess, bounds=bounds, max_nfev = 1000000)
        mpv_uncertainty = np.sqrt(cov[0,0])
        mpv,eta,sigma, A = params[0], params[1], params[2], params[3]
        mpv_uncertainty = np.sqrt(cov[0,0])


    except RuntimeError as e:
        print(f'Fitting failed: {e}')
        # PyLandau can raise an MPV evaluation error as well as scipy reporting
        # optimizer exhaustion. Preserve the actual cause instead of labeling
        # every RuntimeError as nonconvergence or substituting guessed peaks.
        raise ValueError(f"Landau-Gaussian slice fit failed: {e}") from e


    return mpv, eta, sigma, A, mpv_uncertainty, mpv_guess


def langau_lifetime(nhits, dqdx, time_drifted, time_bins, dqdx_bins, nhits_bins, wanted_title):
    """Fit proxy MPV attenuation in drift time, returning tau/error in µs.

    Inputs nhits and dqdx are already divided by segment length. Drift input
    and time-bin edges are in 0.1-µs ticks. Only the historically included
    attenuation slices are fitted: with the standard 19 bins these are 4–17.
    Excluded slices still have time/count diagnostics, but their MPVs/errors
    are null. An empty or pathological unused slice cannot veto this fit or
    influence its initial amplitude. Nominally included slices are usable only
    after the existing support/peak/finite-error checks. At least ten usable
    slices are required. A failed slice is recorded and omitted before the
    attenuation fit; neither its residual nor the resulting lifetime selects
    which bins survive. The attenuation covariance requirement is unchanged.
    """

    time_bins = np.asarray(time_bins)
    fit_points = sum(time_bins <= np.max(time_bins)*.1) + 2
    times = (time_bins[:-1] + time_bins[1:]) / 20  # midpoint, then ticks → µs
    nominal = np.zeros(len(times), dtype=bool)
    nominal[fit_points:-1] = True
    fitted = np.zeros(len(times), dtype=bool)

    # Preserve the original indexing for downstream diagnostics. NaN is only
    # an internal placeholder; finite_diagnostic converts it to strict-JSON
    # null below. Uncomputed excluded peaks must never appear as measurements.
    mpvs = np.full(len(times), np.nan)
    mpv_uncertainties = np.full(len(times), np.nan)
    masks = [(time_drifted >= start) & (time_drifted <= end)
             for start, end in zip(time_bins[:-1], time_bins[1:])]
    slice_status = [dict(bin=int(i), n_segments=int(mask.sum()),
                         status='pending' if nominal[i] else 'excluded',
                         reason=None if nominal[i] else 'outside nominal attenuation range')
                    for i, mask in enumerate(masks)]

    # This loop deliberately cannot call a peak fitter for an excluded slice.
    # Computing all peaks and discarding some afterward allowed an unused
    # bin's empty/one-center histogram to abort otherwise testable data.
    for i in np.flatnonzero(nominal):
        mask = masks[i]
        dq_dx = dqdx[mask]
        status = slice_status[i]
        
        hits = nhits[mask]
        if len(hits) < 5 or not np.any((hits >= nhits_bins[0]) & (hits <= nhits_bins[-1])) or not np.any((dq_dx >= dqdx_bins[0]) & (dq_dx <= dqdx_bins[-1])):
            status.update(status='insufficient_histogram_support',
                          reason='fewer than five segments or no support on a histogram axis')
            continue

        hist_a, edges_a = np.histogram(hits, bins=nhits_bins, density=True)
        hist_b, edges_b = np.histogram(dq_dx, bins=dqdx_bins, density=True)      
        conv, new_edges = convolution(hist_a,  edges_a, hist_b, edges_b)
        conv_bin_centers = (new_edges[:-1] + new_edges[1:])/2
        
        max_value = conv_bin_centers[np.argmax(conv)]
        
        # Historical fit window is 0.5--1.5 times the peak's x coordinate.
        # Despite the old variable name, this is not a measured peak FWHM.
        half_max = max_value/2
        fwhm_range = [max_value-half_max, max_value+half_max]
        
        mask_fit_range = (conv_bin_centers >= fwhm_range[0]) & (conv_bin_centers <= fwhm_range[1])
        
        bin_centers_fit = conv_bin_centers[mask_fit_range]
        h_fit = conv[mask_fit_range]
        status['fit_histogram_centers'] = int(len(bin_centers_fit))

        try:
            mpv, eta, sigma, A, mpv_uncertainty, mpv_guess = langau_fit(h_fit, bin_centers_fit)
        except (ValueError, RuntimeError, FloatingPointError) as exc:
            # Retain a failed included slice without guessing its peak or
            # uncertainty. Eligibility is decided before attenuation fitting,
            # so this cannot remove a bin for disagreeing with a desired tau.
            status.update(status='peak_fit_failed', reason=f'{type(exc).__name__}: {exc}')
            if exc.__cause__ is not None:
                status['inner_exception'] = f'{type(exc.__cause__).__name__}: {exc.__cause__}'
            continue

        if not np.isfinite(mpv) or not np.isfinite(mpv_uncertainty) or mpv_uncertainty <= 0:
            # The full-index MPV/error arrays remain null for an unusable peak.
            # Finite candidate values are still available in its audit record.
            status.update(status='invalid_peak_or_uncertainty',
                          reason='MPV must be finite and its uncertainty finite and positive',
                          candidate_mpv=finite_diagnostic(mpv),
                          candidate_mpv_error=finite_diagnostic(mpv_uncertainty))
            continue

        mpvs[i] = mpv
        mpv_uncertainties[i] = mpv_uncertainty
        fitted[i] = True
        status.update(status='usable', reason=None)

    # Preserve coverage evidence even when too few bins survive to attempt a
    # lifetime fit. All entries are JSON-safe, including nonfinite candidates.
    diagnostics = dict(estimator_version=TRACK_ESTIMATOR_VERSION,
        nominal_fitted_bins=np.flatnonzero(nominal).tolist(),
        fitted_bins=np.flatnonzero(fitted).tolist(),
        minimum_usable_drift_bins=MIN_USABLE_DRIFT_BINS,
        usable_drift_bin_count=int(fitted.sum()),
        excluded_bins=np.flatnonzero(~nominal).tolist(),
        excluded_bin_policy='not peak-fitted; MPV/error null; counts and time centers retained',
        bin_selection_policy='pre-attenuation support and peak validity only; no residual or lifetime selection',
        skipped_bins=[row for row in slice_status if row['status'] not in ('usable', 'excluded')],
        drift_bin_status=slice_status,
        drift_time_us=times.tolist(), mpv=[finite_diagnostic(v) for v in mpvs],
        mpv_error=[finite_diagnostic(v) for v in mpv_uncertainties],
        segments_per_bin=[int(m.sum()) for m in masks],
        uncertainty='weighted attenuation covariance; proxy MPV errors from unweighted slice-fit covariance')
    if fitted.sum() < MIN_USABLE_DRIFT_BINS:
        raise TrackFitError(
            f'Insufficient usable drift bins: {int(fitted.sum())}; require at least {MIN_USABLE_DRIFT_BINS}',
            diagnostics, 'insufficient_usable_drift_bins')

    # Initialization is now determined only by measurements actually fitted.
    # The lifetime initial guess, positive model bounds and absolute-sigma
    # weighting are unchanged. This amplitude change is recorded explicitly.
    initial_guess = [float(np.max(mpvs[fitted])), 2000]
    diagnostics['attenuation_initial_guess'] = initial_guess
    try:
        params, cov = curve_fit(ElectronLifetimeFunc, times[fitted], mpvs[fitted],
                                bounds=([0, 0.], [np.inf, np.inf]), p0=initial_guess,
                                sigma=mpv_uncertainties[fitted], absolute_sigma=True)
    except (ValueError, RuntimeError, FloatingPointError) as exc:
        raise TrackFitError(f'Attenuation fit failed: {exc}', diagnostics,
                            'attenuation_fit_failed') from exc

    # Report goodness-of-fit on exactly the fitted bins. The previous legend
    # accidentally counted the final (excluded) bin and did not subtract the
    # two fitted parameters. This correction changes diagnostics, not tau.
    chi2 = chi_squared1(params, ElectronLifetimeFunc, times[fitted], mpvs[fitted], mpv_uncertainties[fitted])
    ndf = int(fitted.sum()) - len(params)
    x_fit = np.linspace(times[fitted][0], times[fitted][-1], 1000)
    fit_curve = ElectronLifetimeFunc(x_fit, *params)
    
    # Validate before rounding the legacy public error or creating a plot. A
    # genuine sub-microsecond uncertainty may round to zero and is not rejected;
    # a rank-deficient covariance must not masquerade as tau +/- zero instead.
    check = attenuation_covariance_diagnostics(cov)
    diagnostics.update(chi2=finite_diagnostic(chi2), ndf=ndf, chi2_ndf=finite_diagnostic(chi2/ndf),
        candidate_parameters=dict(amplitude=finite_diagnostic(params[0]), lifetime_us=finite_diagnostic(params[1])),
        candidate_error_us_unrounded=(finite_diagnostic(np.sqrt(cov[1, 1]))
            if cov.shape == (2, 2) and np.isfinite(cov[1, 1]) and cov[1, 1] >= 0 else None),
        attenuation_covariance=[[finite_diagnostic(value) for value in row] for row in cov],
        attenuation_covariance_check=check)
    if not check['valid']:
        # The caller records this attempted fit with null accepted lifetime and
        # error. No diagnostic figure is saved that could imply a credible tau.
        raise InvalidAttenuationCovariance('Invalid attenuation covariance: ' + check['reason'], diagnostics)

    fig, ax = plt.subplots()
    # Attach machine-readable diagnostics without changing the public three-
    # value return signature used by the older lifetime entry points.
    fig.lifetime_diagnostics = diagnostics

    # Excluded slices have no fitted MPVs and must not be displayed as gray
    # measurements. Autoscaling includes full uncertainties of included peaks.
    ax.scatter(times[fitted], mpvs[fitted],
               marker=r"$\circ$", color='black', label='Fitted')
    ax.plot(x_fit, fit_curve, color='orange',label=f'$e^{{-}}$ lifetime = {params[1]:.4f} $\\pm$ {np.sqrt(cov[1,1]):.4f} [μs] \n $dQ_{{0}}/dx$ \
            = {params[0]:.4} $\\pm$ {np.sqrt(cov[0,0]):.4f} [proxy units] \n $\\chi^{2}/ndf = {chi2/ndf:.2f}$')
    
    #Axis titles
    ax.set_ylabel(r"MPV of ($\frac{dN}{dx}$ * $\frac{dQ}{dx}$)")
    ax.set_xlabel(r"Time Drifted [$\mu$s]")
    ax.set_title(f"Lifetime: {wanted_title}")

    # Add error bars using matplotlib
    ax.errorbar(x=times[fitted], y=mpvs[fitted], yerr=mpv_uncertainties[fitted],
                fmt='none', ecolor='black', elinewidth=1, capsize=3, zorder=2)
    ax.legend(edgecolor='black',  fontsize= 8, loc='upper right')

    ax.grid(axis='both', linestyle='--')

    return round(params[1], 3), round(np.sqrt(cov[1,1])), fig
