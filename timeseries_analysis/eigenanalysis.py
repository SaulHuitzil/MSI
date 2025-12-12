"""
Eigenvalue analysis using C executables
"""

import subprocess
import os
import pandas as pd
import time
from pathlib import Path

# Configuration - buscar ejecutables en el directorio del paquete
_PACKAGE_DIR = Path(__file__).parent
EIGEN_FULL_EXECUTABLE = str(_PACKAGE_DIR / "eigen_from_cov")
EIGEN_PARTIAL_EXECUTABLE = str(_PACKAGE_DIR / "eigen_from_cov_partial")


def run_eigenvalue_analysis_c(
    cov_file, mode="full", K=None, save_eigenvalues=True, save_eigenvectors=True
):
    """
    Run eigenvalue analysis using C executables

    Args:
        cov_file: Input covariance matrix CSV file
        mode: 'full' or 'partial'
        K: Number of eigenvalues (only for partial mode)
        save_eigenvalues: Whether to save eigenvalues
        save_eigenvectors: Whether to save eigenvectors

    Returns:
        Dictionary with results and timing
    """
    base_name = os.path.splitext(cov_file)[0]
    print(f"Hola ma = {base_name}")

    # Generate output filenames
    if mode == "full":
        eigenval_file = f"{base_name}_eigenvalues_full.csv"
        eigenvec_file = f"{base_name}_eigenvectors_full.csv"
        executable = EIGEN_FULL_EXECUTABLE
        cmd = [executable, cov_file, eigenval_file, eigenvec_file]
        print(f"\nComputing all eigenvalues and eigenvectors...")
    else:  # partial
        if K is None:
            raise ValueError("K must be specified for partial analysis")
        eigenval_file = f"{base_name}_eigenvalues_K{K}.csv"
        eigenvec_file = f"{base_name}_eigenvectors_K{K}.csv"
        executable = EIGEN_PARTIAL_EXECUTABLE
        cmd = [executable, cov_file, str(K), eigenval_file, eigenvec_file]
        print(f"\nComputing the {K} largest eigenvalues and eigenvectors...")

    # Check if executable exists
    if not os.path.exists(executable):
        raise FileNotFoundError(f"Executable not found: {executable}")

    # Run C program
    start_time = time.time()
    subprocess.run(cmd, check=True)
    elapsed_time = time.time() - start_time

    print(f"Time elapsed: {elapsed_time:.2f} seconds")

    # Read results
    eigenvalues = pd.read_csv(eigenval_file)
    eigenvectors = pd.read_csv(eigenvec_file)

    print(f"\nResults ({mode.capitalize()} Analysis):")
    print(f"  - Number of eigenvalues: {len(eigenvalues)}")
    print(f"  - Eigenvector dimensions: {eigenvectors.shape}")
    print(f"  - Largest eigenvalue: {eigenvalues['eigenvalue'].iloc[0]:.6e}")
    print(f"  - Smallest eigenvalue: {eigenvalues['eigenvalue'].iloc[-1]:.6e}")

    results = {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "eigenval_file": eigenval_file if save_eigenvalues else None,
        "eigenvec_file": eigenvec_file if save_eigenvectors else None,
        "time": elapsed_time,
        "mode": mode,
        "K": K if mode == "partial" else len(eigenvalues),
    }

    # Delete files if not saving
    if not save_eigenvalues:
        os.remove(eigenval_file)
        print(f"Eigenvalues not saved (temporary file removed)")
    if not save_eigenvectors:
        os.remove(eigenvec_file)
        print(f"Eigenvectors not saved (temporary file removed)")

    return results
