"""
Time Series Eigenvalue Analysis Package
========================================
Modular package for eigenvalue analysis of time series data.
"""

__version__ = "1.0.0"
__author__ = "saulhuitzil@gmail.com"

# from .io_utils import read_csv_with_comments, save_covariance_matrix, save_matrix_for_c
from .io_utils import read_csv_with_comments, save_matrix_for_c
from .covariance import compute_covariance_matrix
from .eigenanalysis import run_eigenvalue_analysis_c
from .clustering import hierarchical_clustering, get_clustered_order
from .statistics import compute_variance_explained, fit_power_law
from .visualization import (
    plot_eigenspectrum,
    plot_eigenvectors,
    plot_covariance_clustered,
    plot_dendrogram,
)
from .cli import main, analyze_timeseries
from .gui import launch_gui
from .renormalization import renormalize_matrix, validate_block_sizes

__all__ = [
    "read_csv_with_comments",
    "save_covariance_matrix",
    "save_matrix_for_c",
    "compute_covariance_matrix",
    "run_eigenvalue_analysis_c",
    "hierarchical_clustering",
    "get_clustered_order",
    "compute_variance_explained",
    "fit_power_law",
    "plot_eigenspectrum",
    "plot_eigenvectors",
    "plot_covariance_clustered",
    "plot_dendrogram",
    "main",
    "analyze_timeseries",
    "renormalize_matrix",
    "validate_block_sizes",
]
