"""
Renormalization utilities for covariance matrices
"""

import numpy as np


def renormalize_matrix(matrix, block_size):
    """
    Renormalize a matrix by averaging blocks

    Args:
        matrix: Input matrix (n × n)
        block_size: Size of blocks to average (b)

    Returns:
        renormalized_matrix: New matrix of size (n//b × n//b) or larger if n%b != 0
    """
    if block_size == 1:
        return matrix.copy()

    n = matrix.shape[0]

    # Calculate new dimensions
    n_blocks = n // block_size
    remainder = n % block_size

    # New matrix size
    new_n = n_blocks + (1 if remainder > 0 else 0)
    renorm_matrix = np.zeros((new_n, new_n))

    # Average complete blocks
    for i in range(n_blocks):
        for j in range(n_blocks):
            i_start, i_end = i * block_size, (i + 1) * block_size
            j_start, j_end = j * block_size, (j + 1) * block_size
            renorm_matrix[i, j] = np.mean(matrix[i_start:i_end, j_start:j_end])

    # Handle remainder if exists
    if remainder > 0:
        remainder_start = n_blocks * block_size

        # Last row (incomplete blocks × complete blocks)
        for j in range(n_blocks):
            j_start, j_end = j * block_size, (j + 1) * block_size
            renorm_matrix[n_blocks, j] = np.mean(
                matrix[remainder_start:n, j_start:j_end]
            )

        # Last column (complete blocks × incomplete blocks)
        for i in range(n_blocks):
            i_start, i_end = i * block_size, (i + 1) * block_size
            renorm_matrix[i, n_blocks] = np.mean(
                matrix[i_start:i_end, remainder_start:n]
            )

        # Corner (incomplete × incomplete)
        renorm_matrix[n_blocks, n_blocks] = np.mean(
            matrix[remainder_start:n, remainder_start:n]
        )

    return renorm_matrix


def validate_block_sizes(matrix_size, block_sizes):
    """
    Validate that block sizes will produce matrices with at least 2 elements

    Args:
        matrix_size: Size of original matrix (n)
        block_sizes: List of block sizes to validate

    Returns:
        valid_sizes: List of valid block sizes
        invalid_sizes: List of invalid block sizes with reasons
    """
    valid_sizes = []
    invalid_sizes = []

    for b in block_sizes:
        if b < 1:
            invalid_sizes.append((b, "Block size must be >= 1"))
            continue

        new_size = matrix_size // b + (1 if matrix_size % b > 0 else 0)

        if new_size < 2:
            invalid_sizes.append(
                (b, f"Resulting matrix size ({new_size}) would be < 2")
            )
        else:
            valid_sizes.append(b)

    return valid_sizes, invalid_sizes
