"""
Covariance matrix computation
"""

import numpy as np
import time


def compute_covariance_matrix(df, save_binary=False, output_prefix=None):
    """
    Compute covariance matrix from DataFrame

    Args:
        df: DataFrame with time series data
        save_binary: whether to save in binary format
        output_prefix: prefix for binary files (if save_binary=True)

    Returns:
        cov_matrix: Covariance matrix as numpy array
        cor_matrix: Correlation matrix as numpy array
    """
    print("\nComputing covariance matrix...")
    start_time = time.time()

    data = df.values
    cov_matrix = np.cov(data, rowvar=False)

    # Obtener matriz de correlación
    std = np.sqrt(np.diag(cov_matrix))
    corr_matrix = cov_matrix / np.outer(std, std)

    elapsed = time.time() - start_time
    print(f"Covariance matrix computed: {cov_matrix.shape[0]} × {cov_matrix.shape[1]}")
    print(f"Time elapsed: {elapsed:.2f} seconds")

    # Save binary if requested
    if save_binary and output_prefix:
        from .io_utils import save_covariance_binary

        save_covariance_binary(
            cov_matrix, f"{output_prefix}_covariance.bin", df.columns.tolist()
        )
        save_covariance_binary(
            corr_matrix, f"{output_prefix}_correlation.bin", df.columns.tolist()
        )

    return cov_matrix, corr_matrix


def renormalize_covariance_matrix(cov_matrix, block_size):
    """
    Renormalize covariance matrix by averaging blocks

    Args:
        cov_matrix: Original covariance matrix
        block_size: Block size for renormalization

    Returns:
        renorm_cov: Renormalized covariance matrix
    """
    from .renormalization import renormalize_matrix

    return renormalize_matrix(cov_matrix, block_size)
