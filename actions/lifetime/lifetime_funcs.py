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
        raise ValueError("Landau-Gaussian slice fit did not converge") from e


    return mpv, eta, sigma, A, mpv_uncertainty, mpv_guess


def langau_lifetime(nhits, dqdx, time_drifted, time_bins, dqdx_bins, nhits_bins, wanted_title):
    """Fit proxy MPV attenuation in drift time, returning tau/error in µs.

    Inputs nhits and dqdx are already divided by segment length. Drift input
    and time-bin edges are in 0.1-µs ticks. Keep the historical early-bin and
    final-bin exclusions identical in the fit and its diagnostic plot.
    """

    method = " [PyLandau (Landau * Gaussian)]"
    fit_points = sum(time_bins <= np.max(time_bins)*.1) + 2
    
    mpvs = np.empty(len(time_bins)-1)
    times =  np.empty(len(time_bins)-1)

    mpv_uncertainties =  np.empty(len(time_bins)-1)

    masks = []
    convs = []
    edges = []

    for i in range(len(time_bins)-1):

        mask = (time_drifted >= time_bins[i]) & (time_drifted <= time_bins[i+1])
        masks.append(mask)
        dq_dx = dqdx[mask]
        
        hits = nhits[mask]
        if len(hits) < 5 or not np.any((hits >= nhits_bins[0]) & (hits <= nhits_bins[-1])) or not np.any((dq_dx >= dqdx_bins[0]) & (dq_dx <= dqdx_bins[-1])):
            raise ValueError(f"Insufficient histogram support in drift bin {i}")

        hist_a, edges_a = np.histogram(hits, bins=nhits_bins, density=True)
        hist_b, edges_b = np.histogram(dq_dx, bins=dqdx_bins, density=True)      
        conv, new_edges = convolution(hist_a,  edges_a, hist_b, edges_b)
        conv_bin_centers = (new_edges[:-1] + new_edges[1:])/2
        
        convs.append(conv)
        edges.append(new_edges)

        max_value = conv_bin_centers[np.argmax(conv)]
        
        # Historical fit window is 0.5--1.5 times the peak's x coordinate.
        # Despite the old variable name, this is not a measured peak FWHM.
        half_max = max_value/2
        fwhm_range = [max_value-half_max, max_value+half_max]
        
        mask_fit_range = (conv_bin_centers >= fwhm_range[0]) & (conv_bin_centers <= fwhm_range[1])
        
        bin_centers_fit = conv_bin_centers[mask_fit_range]
        h_fit = conv[mask_fit_range]

        mpv, eta, sigma, A, mpv_uncertainty, mpv_guess = langau_fit(h_fit, bin_centers_fit)

        mpvs[i] = mpv
        times[i] = (time_bins[i]+time_bins[i+1])/20  # midpoint (/2), then ticks → µs (/10)
        mpv_uncertainties[i] = mpv_uncertainty

    #Get Electron Lifetime
    if not np.all(np.isfinite(mpvs[fit_points:-1])) or not np.all(np.isfinite(mpv_uncertainties[fit_points:-1])) or np.any(mpv_uncertainties[fit_points:-1] <= 0):
        raise ValueError("Invalid drift-bin MPV uncertainty")
    params, cov = curve_fit(ElectronLifetimeFunc, times[fit_points:-1], mpvs[fit_points:-1], bounds = ([0, .0,], [np.inf, np.inf]), 
                            p0 = [np.max(mpvs), 2000], sigma=mpv_uncertainties[fit_points:-1], absolute_sigma=True)

    # Report goodness-of-fit on exactly the fitted bins. The previous legend
    # accidentally counted the final (excluded) bin and did not subtract the
    # two fitted parameters. This correction changes diagnostics, not tau.
    fitted = np.zeros(len(times), dtype=bool)
    fitted[fit_points:-1] = True
    chi2 = chi_squared1(params, ElectronLifetimeFunc, times[fitted], mpvs[fitted], mpv_uncertainties[fitted])
    ndf = int(fitted.sum()) - len(params)
    x_fit = np.linspace(times[0], times[-1],1000)
    fit_curve = ElectronLifetimeFunc(x_fit, *params)
    
    fig, ax = plt.subplots()
    # Attach machine-readable diagnostics without changing the public three-
    # value return signature used by the older lifetime entry points.
    fig.lifetime_diagnostics = dict(chi2=float(chi2), ndf=ndf, chi2_ndf=float(chi2/ndf),
        drift_time_us=times.tolist(), mpv=mpvs.tolist(), mpv_error=mpv_uncertainties.tolist(),
        fitted_bins=np.flatnonzero(fitted).tolist(), segments_per_bin=[int(m.sum()) for m in masks],
        uncertainty='weighted attenuation covariance; proxy MPV errors from unweighted slice-fit covariance')

    #Set y lim
    ax.set_ylim(min(mpvs)-2,max(mpvs)+2) 

    #Plot
    ax.scatter(times[fitted], mpvs[fitted],
               marker=r"$\circ$", color='black', label='Fitted')
    ax.scatter(times[~fitted], mpvs[~fitted],
               marker=r"$\circ$", color='grey', label='Not Fitted')
    ax.plot(x_fit, fit_curve, color='orange',label=f'$e^{{-}}$ lifetime = {params[1]:.4f} $\\pm$ {np.sqrt(cov[1,1]):.4f} [μs] \n $dQ_{{0}}/dx$ \
            = {params[0]:.4} $\\pm$ {np.sqrt(cov[0,0]):.4f} [proxy units] \n $\\chi^{2}/ndf = {chi2/ndf:.2f}$')
    
    #Axis titles
    ax.set_ylabel(r"MPV of ($\frac{dN}{dx}$ * $\frac{dQ}{dx}$)")
    ax.set_xlabel(r"Time Drifted [$\mu$s]")
    ax.set_title(f"Lifetime: {wanted_title}")

    # Add error bars using matplotlib
    ax.errorbar(x=times[fitted], y=mpvs[fitted], yerr=mpv_uncertainties[fitted],
                fmt='none', ecolor='black', elinewidth=1, capsize=3, zorder=2)
    ax.errorbar(x=times[~fitted], y=mpvs[~fitted], yerr=mpv_uncertainties[~fitted],
                fmt='none', ecolor='grey', elinewidth=1, capsize=3, zorder=2)

    ax.legend(edgecolor='black',  fontsize= 8, loc='upper right')

    ax.grid(axis='both', linestyle='--')

    return round(params[1], 3), round(np.sqrt(cov[1,1])), fig
