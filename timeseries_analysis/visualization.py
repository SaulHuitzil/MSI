"""
Visualization functions for eigenvalue analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from .statistics import compute_variance_explained, fit_power_law
from .clustering import hierarchical_clustering, get_clustered_order
from .io_utils import save_variable_order


def plot_eigenspectrum(
    results, metadata, output_file=None, show_variance=True, n_skip=10
):
    """
    Plot eigenvalue spectrum with power law fit and variance markers

    Args:
        results: Dictionary with eigenvalue analysis results
        metadata: Dictionary with metadata
        output_file: Output filename (if None, don't save)
        show_variance: Whether to show variance markers
        n_skip: Number of eigenvalues to skip in power law fit
    """
    eigenvalues = results["eigenvalues"]
    eigenvals = eigenvalues["eigenvalue"].values
    mode = results["mode"]
    K = results["K"]

    # Compute variance explained
    explained_var, idx_90, idx_95, idx_99 = compute_variance_explained(eigenvals)

    # Fit power law
    slope, r2, fit, indices, evals = fit_power_law(
        eigenvals, fit_fraction=0.9, n_skip=n_skip
    )

    # Create figure
    fig = plt.figure(figsize=(14, 6))
    ax = fig.add_subplot(111)

    # Plot eigenvalues
    ax.loglog(
        indices,
        evals,
        "o",
        markersize=4,
        alpha=0.6,
        label="Eigenvalues",
        color="steelblue",
    )

    # Power law fit
    ax.loglog(
        indices,
        fit,
        "r-",
        linewidth=2.5,
        alpha=0.8,
        label=f"Power law: λ ∝ k$^{{{slope:.2f}}}$ (R²={r2:.4f})",
    )

    # Show variance markers only for full analysis
    if mode == "full" and show_variance and False:
        for threshold, idx_val, color, label in [
            (0.90, idx_90, "green", "90% variance"),
            (0.95, idx_95, "orange", "95% variance"),
            (0.99, idx_99, "red", "99% variance"),
        ]:
            if idx_val <= len(evals):
                ax.axvline(
                    x=idx_val,
                    color=color,
                    linestyle="--",
                    linewidth=2,
                    alpha=0.7,
                    label=f"{label} (k={idx_val})",
                )

    ax.set_xlabel("Eigenvalue index (k)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Eigenvalue (λ)", fontsize=14, fontweight="bold")

    title = "Eigenvalue Spectrum"
    if mode == "partial":
        title += f" (K={K} largest eigenvalues)"
    ax.set_title(title, fontsize=15, fontweight="bold")

    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="best", fontsize=11, framealpha=0.9)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=250, bbox_inches="tight")
        print(f"Spectrum plot saved: {output_file}")

    return fig


def determine_eigenvector_sign(eigenvec_data):
    """
    Determine the sign correction for an eigenvector using majority voting

    Three criteria are used:
    1. First component should be positive
    2. Largest (in absolute value) component should be positive
    3. Average of all components should be positive

    Args:
        eigenvec_data: numpy array with eigenvector components

    Returns:
        sign_factor: +1 or -1 to multiply the eigenvector
    """
    votes = []

    # Criterion 1: First component should be positive
    if eigenvec_data[0] >= 0:
        votes.append(1)
    else:
        votes.append(-1)

    # Criterion 2: Largest component (in absolute value) should be positive
    max_abs_idx = np.argmax(np.abs(eigenvec_data))
    if eigenvec_data[max_abs_idx] >= 0:
        votes.append(1)
    else:
        votes.append(-1)

    # Criterion 3: Average should be positive
    avg = np.mean(eigenvec_data)
    if avg >= 0:
        votes.append(1)
    else:
        votes.append(-1)

    # Majority vote
    sign_factor = 1 if sum(votes) >= 0 else -1

    return sign_factor


def plot_eigenvectors(results, metadata, n_vecs=3, output_file=None):
    """
    Plot the first n eigenvectors with consistent sign convention

    Args:
        results: Dictionary with eigenvalue analysis results
        metadata: Dictionary with metadata
        n_vecs: Number of eigenvectors to plot
        output_file: Output filename (if None, don't save)
    """
    from scipy import stats

    eigenvectors = results["eigenvectors"]
    eigenvalues = results["eigenvalues"]
    eigenvals = eigenvalues["eigenvalue"].values
    total_variance = np.sum(eigenvals)

    # Adjust n_vecs if necessary
    n_vecs = min(n_vecs, len(eigenvals))

    fig, axes = plt.subplots(n_vecs, 1, figsize=(14, 4.0 * n_vecs))

    # Handle single subplot case
    if n_vecs == 1:
        axes = [axes]

    for i in range(n_vecs):
        ax = axes[i]
        eigenvec_data = eigenvectors[f"ev_{i+1}"].values

        # Apply sign convention using majority voting
        sign_factor = determine_eigenvector_sign(eigenvec_data)
        eigenvec_data_corrected = sign_factor * eigenvec_data

        components = np.arange(len(eigenvec_data_corrected))

        ax.bar(
            components,
            eigenvec_data_corrected,
            width=1.0,
            alpha=0.7,
            color=f"C{i}",
            edgecolor="none",
        )
        ax.axhline(y=0, color="k", linewidth=0.8, linestyle="-", alpha=0.5)
        ax.set_xlabel("Component", fontsize=13, fontweight="bold")
        ax.set_ylabel("Value", fontsize=13, fontweight="bold")

        # Calculate statistical moments
        mean = np.mean(eigenvec_data_corrected)
        variance = np.var(eigenvec_data_corrected, ddof=1)
        skewness = stats.skew(eigenvec_data_corrected)
        kurt = stats.kurtosis(eigenvec_data_corrected)

        # Main title with eigenvalue info
        title_line1 = f"Eigenvector {i+1} (λ = {eigenvals[i]:.4e})"

        # Subtitle with statistical moments
        title_line2 = (
            f"Stat. Moments = ({mean:.3e}, {variance:.3e}, {skewness:.3e}, {kurt:.3e})"
        )

        ax.set_title(
            f"{title_line1}\n{title_line2}",
            fontsize=13,
            fontweight="bold",
        )
        ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=250, bbox_inches="tight")
        print(f"Eigenvector plot saved: {output_file}")

    return fig


# def plot_eigenvectors(results, metadata, n_vecs=3, output_file=None):
#     """
#     Plot the first n eigenvectors with consistent sign convention

#     Args:
#         results: Dictionary with eigenvalue analysis results
#         metadata: Dictionary with metadata
#         n_vecs: Number of eigenvectors to plot
#         output_file: Output filename (if None, don't save)
#     """
#     eigenvectors = results["eigenvectors"]
#     eigenvalues = results["eigenvalues"]
#     eigenvals = eigenvalues["eigenvalue"].values
#     total_variance = np.sum(eigenvals)

#     # Adjust n_vecs if necessary
#     n_vecs = min(n_vecs, len(eigenvals))

#     fig, axes = plt.subplots(n_vecs, 1, figsize=(14, 3.5 * n_vecs))

#     # Handle single subplot case
#     if n_vecs == 1:
#         axes = [axes]

#     for i in range(n_vecs):
#         ax = axes[i]
#         eigenvec_data = eigenvectors[f"ev_{i+1}"].values

#         # Apply sign convention
#         sign_factor = determine_eigenvector_sign(eigenvec_data)
#         eigenvec_data_corrected = sign_factor * eigenvec_data

#         components = np.arange(len(eigenvec_data_corrected))

#         ax.bar(
#             components,
#             eigenvec_data_corrected,
#             width=1.0,
#             alpha=0.7,
#             color=f"C{i}",
#             edgecolor="none",
#         )
#         ax.axhline(y=0, color="k", linewidth=0.8, linestyle="-", alpha=0.5)
#         ax.set_xlabel("Component", fontsize=13, fontweight="bold")
#         ax.set_ylabel("Value", fontsize=13, fontweight="bold")

#         ax.set_title(
#             f"Eigenvector {i+1} (λ = {eigenvals[i]:.4e})",
#             fontsize=13,
#             fontweight="bold",
#         )
#         ax.grid(True, alpha=0.3, axis="y")

#     plt.tight_layout()

#     if output_file:
#         plt.savefig(output_file, dpi=250, bbox_inches="tight")
#         print(f"Eigenvector plot saved: {output_file}")

#     return fig


def plot_covariance_clustered(
    cov_matrix, variable_names=None, output_file=None, tipo_color="viridis"
):
    """
    Plot covariance matrix with hierarchical clustering reordering

    Args:
        cov_matrix: Covariance matrix as numpy array
        variable_names: Optional list of variable names
        output_file: Output filename (if None, don't save)
    """
    n = cov_matrix.shape[0]

    if variable_names is None:
        variable_names = [f"var_{i}" for i in range(n)]

    # Perform clustering
    Z, dist, corr_matrix = hierarchical_clustering(cov_matrix)
    ordered_labels, ordered_indices = get_clustered_order(Z, variable_names)

    # Reorder covariance matrix
    cov_reordered = cov_matrix[np.ix_(ordered_indices, ordered_indices)]

    # Rango de color sin usar la diagonal
    vmin, vmax = _limits_without_diagonal(cov_matrix, percentiles=(3, 97))

    # Determine how many labels to show
    n_labels_show = min(75, n)
    step = max(1, n // n_labels_show)
    tick_positions = list(range(0, n, step))

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Original
    # im0 = axes[0].imshow(cov_matrix, cmap=tipo_color)
    im0 = axes[0].imshow(cov_matrix, cmap=tipo_color, vmin=vmin, vmax=vmax)
    axes[0].set_title(
        "Covariance Matrix (Original Order)", fontsize=13, fontweight="bold"
    )
    axes[0].set_xlabel("Variable", fontsize=11)
    axes[0].set_ylabel("Variable", fontsize=11)
    axes[0].set_xticks(tick_positions)
    axes[0].set_yticks(tick_positions)
    axes[0].set_xticklabels(
        [variable_names[i] for i in tick_positions], rotation=90, fontsize=6
    )
    axes[0].set_yticklabels([variable_names[i] for i in tick_positions], fontsize=6)

    # Reordered
    # im1 = axes[1].imshow(cov_reordered, cmap=tipo_color)
    im1 = axes[1].imshow(cov_reordered, cmap=tipo_color, vmin=vmin, vmax=vmax)
    axes[1].set_title("Covariance Matrix (Clustered)", fontsize=13, fontweight="bold")
    axes[1].set_xlabel("Variable", fontsize=11)
    axes[1].set_ylabel("Variable", fontsize=11)
    axes[1].set_xticks(tick_positions)
    axes[1].set_yticks(tick_positions)
    axes[1].set_xticklabels(
        [ordered_labels[i] for i in tick_positions], rotation=90, fontsize=6
    )
    axes[1].set_yticklabels([ordered_labels[i] for i in tick_positions], fontsize=6)

    cbar = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label("Covariance", fontsize=11)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=250, bbox_inches="tight")
        print(f"Clustered covariance plot saved: {output_file}")

        # Save variable ordering
        order_file = output_file.replace(".png", "_variable_order.csv")
        # save_variable_order(variable_names, ordered_labels, order_file)

    return fig


def _limits_without_diagonal(matrix, percentiles=None):
    m = matrix.astype(float).copy()
    np.fill_diagonal(m, np.nan)
    if percentiles is None:
        vmin = np.nanmin(m)
        vmax = np.nanmax(m)
    else:
        p_low, p_high = percentiles
        vmin = np.nanpercentile(m, p_low)
        vmax = np.nanpercentile(m, p_high)
    return vmin, vmax


def plot_correlation_clustered(
    cor_matrix, variable_names=None, output_file=None, tipo_color="viridis"
):
    """
    Plot correlation matrix with hierarchical clustering reordering

    Args:
        cor_matrix: Correlation matrix as numpy array
        variable_names: Optional list of variable names
        output_file: Output filename (if None, don't save)
    """
    n = cor_matrix.shape[0]

    if variable_names is None:
        variable_names = [f"var_{i}" for i in range(n)]

    # Perform clustering
    Z, dist, corr_matrix = hierarchical_clustering(cor_matrix)
    ordered_labels, ordered_indices = get_clustered_order(Z, variable_names)

    # Reorder correlation matrix
    cor_reordered = cor_matrix[np.ix_(ordered_indices, ordered_indices)]

    # Determine how many labels to show
    n_labels_show = min(75, n)
    step = max(1, n // n_labels_show)
    tick_positions = list(range(0, n, step))

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))

    # Original
    im0 = axes[0].imshow(cor_matrix, cmap=tipo_color)
    axes[0].set_title(
        "Correlation Matrix (Original Order)", fontsize=13, fontweight="bold"
    )
    axes[0].set_xlabel("Variable", fontsize=11)
    axes[0].set_ylabel("Variable", fontsize=11)
    axes[0].set_xticks(tick_positions)
    axes[0].set_yticks(tick_positions)
    axes[0].set_xticklabels(
        [variable_names[i] for i in tick_positions], rotation=90, fontsize=6
    )
    axes[0].set_yticklabels([variable_names[i] for i in tick_positions], fontsize=6)

    # Reordered
    im1 = axes[1].imshow(cor_reordered, cmap=tipo_color)
    axes[1].set_title("Correlation Matrix (Clustered)", fontsize=13, fontweight="bold")
    axes[1].set_xlabel("Variable", fontsize=11)
    axes[1].set_ylabel("Variable", fontsize=11)
    axes[1].set_xticks(tick_positions)
    axes[1].set_yticks(tick_positions)
    axes[1].set_xticklabels(
        [ordered_labels[i] for i in tick_positions], rotation=90, fontsize=6
    )
    axes[1].set_yticklabels([ordered_labels[i] for i in tick_positions], fontsize=6)

    cbar = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    cbar.set_label("Correlation", fontsize=11)

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=250, bbox_inches="tight")
        print(f"Clustered correlation plot saved: {output_file}")

        # Save variable ordering
        order_file = output_file.replace(".png", "_variable_order.csv")
        # save_variable_order(variable_names, ordered_labels, order_file)

    return fig


def plot_dendrogram(cov_matrix, variable_names=None, output_file=None):
    """
    Plot dendrogram from hierarchical clustering of covariance matrix

    Args:
        cov_matrix: Covariance matrix as numpy array
        variable_names: Optional list of variable names
        output_file: Output filename (if None, don't save)
    """
    n = cov_matrix.shape[0]

    if variable_names is None:
        variable_names = [f"var_{i}" for i in range(n)]

    # Perform clustering
    Z, dist, corr_matrix = hierarchical_clustering(cov_matrix)

    # Plot dendrogram
    fig, ax = plt.subplots(figsize=(25, 8))

    dendrogram(
        Z,
        labels=variable_names,
        ax=ax,
        leaf_rotation=90,
        leaf_font_size=8,
        color_threshold=0.7 * max(Z[:, 2]),
    )

    ax.set_title("Hierarchical Clustering Dendrogram", fontsize=15, fontweight="bold")
    ax.set_xlabel("Variable", fontsize=13, fontweight="bold")
    ax.set_ylabel("Distance (1 - |correlation|)", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=250, bbox_inches="tight")
        print(f"Dendrogram saved: {output_file}")

    return fig
