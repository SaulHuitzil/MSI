"""
Statistical analysis functions
"""

import numpy as np
from scipy.stats import linregress


def compute_variance_explained(eigenvals):
    """
    Compute cumulative variance explained and key thresholds

    Args:
        eigenvals: Array of eigenvalues (sorted descending)

    Returns:
        explained_variance: Cumulative variance explained
        idx_90, idx_95, idx_99: Indices for 90%, 95%, 99% thresholds
    """
    total_variance = np.sum(eigenvals)
    explained_variance = np.cumsum(eigenvals) / total_variance

    idx_90 = (
        np.argmax(explained_variance >= 0.90) + 1
        if np.any(explained_variance >= 0.90)
        else len(eigenvals)
    )
    idx_95 = (
        np.argmax(explained_variance >= 0.95) + 1
        if np.any(explained_variance >= 0.95)
        else len(eigenvals)
    )
    idx_99 = (
        np.argmax(explained_variance >= 0.99) + 1
        if np.any(explained_variance >= 0.99)
        else len(eigenvals)
    )

    return explained_variance, idx_90, idx_95, idx_99


def fit_power_law(eigenvals, fit_fraction=0.6, n_skip=0):
    """
    Fit power law to eigenvalue spectrum

    Args:
        eigenvals: Array of eigenvalues (sorted descending)
        fit_fraction: Fraction of data to use for fitting
        n_skip: Number of leading eigenvalues to exclude from the fit

    Returns:
        slope, r_squared, power_law_fit, indices, positive_eigenvals
    """
    # Filter positive eigenvalues
    positive_mask = eigenvals > 0
    positive_eigenvals = eigenvals[positive_mask]
    indices = np.arange(1, len(positive_eigenvals) + 1)

    if n_skip >= len(positive_eigenvals):
        raise ValueError(
            "n_skip is greater than or equal to number of positive eigenvalues"
        )

    # Log-log transform (skipping first n_skip)
    fit_indices = indices[n_skip:]
    fit_eigenvals = positive_eigenvals[n_skip:]

    log_indices = np.log10(fit_indices)
    log_eigenvals = np.log10(fit_eigenvals)

    # ==== 20251128 18:21
    n_total = len(log_indices)
    fit_range_fraction = int(fit_fraction * n_total)

    # Número de puntos deseados por límite superior
    max_points = 1000

    # Regla combinada
    if n_total >= max_points:
        fit_range = max_points
    else:
        fit_range = fit_range_fraction
    # ==== 20251128 18:21

    # Fit power law
    fit_range = int(fit_fraction * len(log_indices))
    slope, intercept, r_value, p_value, std_err = linregress(
        log_indices[:fit_range], log_eigenvals[:fit_range]
    )

    r_squared = r_value**2
    power_law_fit = 10**intercept * indices**slope

    return slope, r_squared, power_law_fit, indices, positive_eigenvals
