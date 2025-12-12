"""
I/O utilities for time series analysis
"""

import os
import pandas as pd
import numpy as np


def read_csv_with_comments(filename):
    """
    Read CSV file with flexible comment lines (starting with #)

    Args:
        filename: Path to CSV file

    Returns:
        df: DataFrame with time series data
        metadata: Dictionary with metadata from comment lines
    """
    metadata = {}

    # Read metadata from comment lines
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                if "=" in line:
                    parts = line[1:].strip().split("=", 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        value = parts[1].strip()
                        metadata[key] = value
            else:
                break

    # Read data using pandas
    df = pd.read_csv(filename, comment="#")

    print(f"Data loaded: {df.shape[0]} samples × {df.shape[1]} variables")
    if metadata:
        print(f"Metadata found: {list(metadata.keys())}")

    return df, metadata


# def save_covariance_matrix2(cov_matrix, filename, variable_names=None):
#     """
#     Save covariance matrix to CSV file

#     Args:
#         cov_matrix: Covariance matrix as numpy array
#         filename: Output filename
#         variable_names: Optional list of variable names
#     """
#     n = cov_matrix.shape[0]

#     if variable_names is None:
#         variable_names = [f"var_{i}" for i in range(n)]

#     df_cov = pd.DataFrame(cov_matrix, columns=variable_names, index=variable_names)
#     df_cov.to_csv(filename, float_format="%.5f")
#     print(f"Covariance matrix saved: {filename}")


# def save_correlation_matrix2(cor_matrix, filename, variable_names=None):
#     """
#     Save correlation matrix to CSV file

#     Args:
#         cor_matrix: Correlation matrix as numpy array
#         filename: Output filename
#         variable_names: Optional list of variable names
#     """
#     n = cor_matrix.shape[0]

#     if variable_names is None:
#         variable_names = [f"var_{i}" for i in range(n)]

#     df_cor = pd.DataFrame(cor_matrix, columns=variable_names, index=variable_names)
#     df_cor.to_csv(filename, float_format="%.5f")
#     print(f"Correlation matrix saved: {filename}")


def save_matrix_for_c(cov_matrix, filename):
    """
    Save covariance matrix in format expected by C programs (CSV)

    Args:
        cov_matrix: Covariance matrix as numpy array
        filename: Output filename
    """
    n = cov_matrix.shape[0]

    with open(filename, "w") as f:
        f.write(",".join([f"var_{i}" for i in range(n)]))
        f.write("\n")

        for i in range(n):
            f.write(",".join([f"{cov_matrix[i, j]:.5e}" for j in range(n)]))
            f.write("\n")

    print(f"Matrix saved for C analysis: {filename}")


def save_variable_order(variable_names, ordered_labels, filename):
    """
    Save variable ordering to CSV file

    Args:
        variable_names: Original variable names
        ordered_labels: Clustered variable names
        filename: Output filename
    """
    n = len(variable_names)
    order_df = pd.DataFrame(
        {
            "original_index": list(range(n)),
            "original_name": variable_names,
            "clustered_index": [
                ordered_labels.index(name) if name in ordered_labels else -1
                for name in variable_names
            ],
            "clustered_name": ordered_labels,
        }
    )
    order_df.to_csv(filename, index=False)
    print(f"Variable ordering saved: {filename}")


# ============================================================================
# BINARY FORMAT I/O FUNCTIONS
# ============================================================================


def save_matrix_binary(matrix, filename, dtype=np.float32):
    """
    Save matrix in binary format (float32 by default)

    Format: [n (int32), m (int32), data (float32 or float64 array)]

    Args:
        matrix: numpy array
        filename: output .bin file
        dtype: data type (np.float32 or np.float64)
    """
    matrix_typed = matrix.astype(dtype)
    with open(filename, "wb") as f:
        # Write dimensions first (2 integers)
        n, m = matrix_typed.shape
        np.array([n, m], dtype=np.int32).tofile(f)
        # Write matrix data
        matrix_typed.tofile(f)

    file_size_mb = os.path.getsize(filename) / (1024 * 1024)
    print(f"Binary matrix saved: {filename} ({n}×{m}, {dtype}, {file_size_mb:.2f} MB)")


def read_matrix_binary(filename, dtype=np.float32):
    """
    Read matrix from binary format

    Format: [n (int32), m (int32), data (float32 or float64 array)]

    Args:
        filename: input .bin file
        dtype: data type (np.float32 or np.float64)

    Returns:
        matrix: numpy array
    """
    with open(filename, "rb") as f:
        # Read dimensions
        dims = np.fromfile(f, dtype=np.int32, count=2)
        n, m = dims
        # Read matrix
        matrix = np.fromfile(f, dtype=dtype, count=n * m)
        matrix = matrix.reshape((n, m))

    print(f"Binary matrix loaded: {filename} ({n}×{m})")
    return matrix


def save_covariance_binary(cov_matrix, filename, variable_names=None):
    """
    Save covariance matrix in binary format + optional variable names

    Args:
        cov_matrix: numpy array
        filename: base filename (will create .bin and optionally .names)
        variable_names: optional list of variable names
    """
    # Ensure .bin extension
    if not filename.endswith(".bin"):
        bin_file = filename + ".bin"
    else:
        bin_file = filename

    # Save binary matrix
    save_matrix_binary(cov_matrix, bin_file, dtype=np.float32)

    # Save variable names separately if provided
    if variable_names is not None:
        names_file = bin_file.replace(".bin", ".names")
        with open(names_file, "w") as f:
            f.write("\n".join(variable_names))
        print(f"Variable names saved: {names_file} ({len(variable_names)} names)")


def read_covariance_binary(filename, has_names=True):
    """
    Read covariance matrix from binary format

    Args:
        filename: .bin file
        has_names: whether to load .names file

    Returns:
        matrix: numpy array
        variable_names: list or None
    """
    # Read binary matrix
    matrix = read_matrix_binary(filename, dtype=np.float32)

    # Try to load variable names if requested
    variable_names = None
    if has_names:
        names_file = filename.replace(".bin", ".names")
        if os.path.exists(names_file):
            with open(names_file, "r") as f:
                variable_names = [line.strip() for line in f if line.strip()]
            print(f"Variable names loaded: {len(variable_names)} names")
        else:
            print(f"Warning: .names file not found: {names_file}")
            print(f"Continuing without variable names...")

    return matrix, variable_names


# def save_correlation_binary2(cor_matrix, filename, variable_names=None):
#     """
#     Save correlation matrix in binary format + optional variable names

#     Args:
#         cor_matrix: numpy array
#         filename: base filename (will create .bin and optionally .names)
#         variable_names: optional list of variable names
#     """
#     # Ensure .bin extension
#     if not filename.endswith(".bin"):
#         bin_file = filename + ".bin"
#     else:
#         bin_file = filename

#     # Save binary matrix
#     save_matrix_binary(cor_matrix, bin_file, dtype=np.float32)

#     # Save variable names separately if provided
#     if variable_names is not None:
#         names_file = bin_file.replace(".bin", ".names")
#         with open(names_file, "w") as f:
#             f.write("\n".join(variable_names))
#         print(f"Variable names saved: {names_file} ({len(variable_names)} names)")


# def read_correlation_binary2(filename, has_names=True):
#     """
#     Read correlation matrix from binary format

#     Args:
#         filename: .bin file
#         has_names: whether to load .names file

#     Returns:
#         matrix: numpy array
#         variable_names: list or None
#     """
#     # Same implementation as read_covariance_binary
#     return read_covariance_binary(filename, has_names)


def convert_csv_to_binary(csv_file, output_file=None, has_header=True, has_index=True):
    """
    Convert CSV covariance/correlation matrix to binary format

    Args:
        csv_file: input CSV file
        output_file: output .bin file (if None, uses csv_file name)
        has_header: whether CSV has column headers
        has_index: whether CSV has row index

    Returns:
        output_file: path to created binary file
    """
    if output_file is None:
        output_file = csv_file.replace(".csv", ".bin")

    # Read CSV
    if has_header and has_index:
        df = pd.read_csv(csv_file, index_col=0)
        variable_names = df.columns.tolist()
    elif has_header and not has_index:
        df = pd.read_csv(csv_file)
        variable_names = df.columns.tolist()
    elif not has_header and has_index:
        df = pd.read_csv(csv_file, header=None, index_col=0)
        variable_names = None
    else:  # no header, no index
        df = pd.read_csv(csv_file, header=None)
        variable_names = None

    matrix = df.values

    # Save as binary
    save_matrix_binary(matrix, output_file, dtype=np.float32)

    # Save variable names if available
    if variable_names is not None:
        names_file = output_file.replace(".bin", ".names")
        with open(names_file, "w") as f:
            f.write("\n".join(variable_names))
        print(f"Variable names saved: {names_file}")

    return output_file


def convert_binary_to_csv(bin_file, output_file=None, has_names=True):
    """
    Convert binary matrix to CSV format

    Args:
        bin_file: input .bin file
        output_file: output CSV file (if None, uses bin_file name)
        has_names: whether to load .names file

    Returns:
        output_file: path to created CSV file
    """
    if output_file is None:
        output_file = bin_file.replace(".bin", ".csv")

    # Read binary
    matrix, variable_names = read_covariance_binary(bin_file, has_names)

    # Create DataFrame
    if variable_names is not None:
        df = pd.DataFrame(matrix, columns=variable_names, index=variable_names)
    else:
        n = matrix.shape[0]
        var_names = [f"var_{i}" for i in range(n)]
        df = pd.DataFrame(matrix, columns=var_names, index=var_names)

    # Save CSV
    df.to_csv(output_file, float_format="%.5e")
    print(f"CSV matrix saved: {output_file}")

    return output_file
