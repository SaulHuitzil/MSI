"""
Command line interface and main workflow
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from .io_utils import (
    read_csv_with_comments,
    save_matrix_for_c,
    save_covariance_binary,
    # save_correlation_binary,
    read_covariance_binary,
    save_matrix_binary,
)
from .covariance import compute_covariance_matrix
from .eigenanalysis import run_eigenvalue_analysis_c
from .statistics import compute_variance_explained, fit_power_law
from .visualization import (
    plot_eigenspectrum,
    plot_eigenvectors,
    plot_covariance_clustered,
    plot_correlation_clustered,
    plot_dendrogram,
)
from .renormalization import renormalize_matrix, validate_block_sizes


mis_colores = [
    "binary",
    "twilight_shifted",
    "RdBu",
    "BuPu",
    "ocean",
]

umbral_JerClust = 4000


def analyze_timeseries(
    input_file=None,
    mode="full",
    K=None,
    import_covariance=None,
    import_correlation=None,
    has_names=True,
    save_covariance=False,
    save_binary=True,
    save_eigenvalues=True,
    save_eigenvectors=True,
    plot_spectrum=True,
    plot_eigvecs=True,
    plot_cov_clustered=True,
    plot_dendro=True,
    n_eigvecs=3,
    save_plots=True,
    output_dir=None,
    block_sizes=None,
):
    """
    Complete workflow for time series eigenvalue analysis with renormalization

    Args:
        input_file: Path to input CSV file (optional if importing matrices)
        mode: 'full' or 'partial' eigenvalue analysis
        K: Number of eigenvalues (for partial mode)
        import_covariance: Path to covariance .bin file (optional)
        import_correlation: Path to correlation .bin file (optional)
        has_names: Whether imported matrices have .names files
        save_covariance: Whether to save covariance matrix
        save_binary: Whether to use binary format for C communication
        save_eigenvalues: Whether to save eigenvalues
        save_eigenvectors: Whether to save eigenvectors
        plot_spectrum: Whether to plot eigenvalue spectrum
        plot_eigvecs: Whether to plot eigenvectors
        plot_cov_clustered: Whether to plot clustered covariance matrix
        plot_dendro: Whether to plot dendrogram
        n_eigvecs: Number of eigenvectors to plot
        save_plots: Whether to save plots
        output_dir: Output directory (if None, use input file directory)
        block_sizes: List of block sizes for renormalization (default: [1])

    Returns:
        Dictionary with all results
    """
    print("=" * 70)
    print("TIME SERIES EIGENVALUE ANALYSIS")
    print("=" * 70)

    # Default block sizes
    if block_sizes is None:
        block_sizes = [1]

    # Determine if we're importing or computing
    importing_matrices = import_covariance is not None

    # Setup base output directory
    if importing_matrices:
        base_name = os.path.splitext(os.path.basename(import_covariance))[0]
        base_name = base_name.replace("_covariance", "").replace("_correlation", "")
    else:
        if input_file is None:
            raise ValueError("Either input_file or import_covariance must be provided")
        base_name = os.path.splitext(os.path.basename(input_file))[0]

    if output_dir is None:
        if importing_matrices:
            parent_dir = os.path.dirname(import_covariance) or "."
        else:
            parent_dir = os.path.dirname(input_file) or "."
        output_dir = os.path.join(parent_dir, f"{base_name}_analysis")

    os.makedirs(output_dir, exist_ok=True)
    print(f"Base output directory: {output_dir}")

    # Metadata dictionary
    metadata = {}

    # ========================================================================
    # LOAD OR COMPUTE COVARIANCE/CORRELATION MATRICES
    # ========================================================================

    if importing_matrices:
        # Import existing covariance matrix
        print(f"\nImporting covariance matrix from: {import_covariance}")
        cov_matrix_original, variable_names = read_covariance_binary(
            import_covariance, has_names=has_names
        )

        # Import correlation matrix if provided, otherwise compute from covariance
        if import_correlation:
            print(f"Importing correlation matrix from: {import_correlation}")
            cor_matrix_original, _ = read_covariance_binary(
                import_correlation, has_names=has_names
            )
        else:
            print("Computing correlation matrix from covariance...")
            std = np.sqrt(np.diag(cov_matrix_original))
            cor_matrix_original = cov_matrix_original / np.outer(std, std)

        # Use imported variable names or create default
        if variable_names is None:
            n = cov_matrix_original.shape[0]
            variable_names = [f"var_{i}" for i in range(n)]
            print(
                f"No variable names found, using default: var_0, var_1, ..., var_{n-1}"
            )

        metadata["source"] = "imported"
        metadata["covariance_file"] = import_covariance

    else:
        # Original workflow: read data and compute matrices
        print(f"\nReading data from: {input_file}")
        df, metadata = read_csv_with_comments(input_file)

        # Compute covariance matrix
        print("\nComputing covariance and correlation matrices...")
        cov_matrix_original, cor_matrix_original = compute_covariance_matrix(df)

        variable_names = df.columns.tolist()
        metadata["source"] = "computed"
        metadata["input_file"] = input_file

    # ========================================================================
    # VALIDATE BLOCK SIZES
    # ========================================================================

    original_size = cov_matrix_original.shape[0]
    valid_sizes, invalid_sizes = validate_block_sizes(original_size, block_sizes)

    if invalid_sizes:
        print("\n" + "!" * 70)
        print("WARNING: Some block sizes are invalid and will be skipped:")
        for size, reason in invalid_sizes:
            print(f"  K={size}: {reason}")
        print("!" * 70)

    if not valid_sizes:
        raise ValueError("No valid block sizes to process!")

    print(f"\nValid block sizes to process: {valid_sizes}")

    # ========================================================================
    # SAVE ORIGINAL MATRICES (only for K=1, if requested)
    # ========================================================================

    if save_covariance and 1 in valid_sizes:
        k01_dir = os.path.join(output_dir, "K01")
        os.makedirs(k01_dir, exist_ok=True)

        if save_binary:
            # Save in binary format
            cov_bin_file = os.path.join(k01_dir, "covariance.bin")
            # cor_bin_file = os.path.join(k01_dir, "correlation.bin")
            save_covariance_binary(cov_matrix_original, cov_bin_file, variable_names)
            # save_correlation_binary(cor_matrix_original, cor_bin_file, variable_names)

    # ========================================================================
    # PROCESS EACH BLOCK SIZE
    # ========================================================================

    all_results = {}

    for block_size in valid_sizes:
        print("\n" + "=" * 70)
        print(f"PROCESSING RENORMALIZATION WITH K={block_size}")
        print("=" * 70)

        # Create output directory for this block size
        k_dir = os.path.join(output_dir, f"K{block_size:02d}")
        os.makedirs(k_dir, exist_ok=True)
        print(f"Output directory: {k_dir}")

        # ====================================================================
        # RENORMALIZE MATRICES
        # ====================================================================

        if block_size == 1:
            cov_matrix = cov_matrix_original
            cor_matrix = cor_matrix_original
            print(f"Using original matrices (no renormalization)")
        else:
            print(f"\nRenormalizing matrices with block size {block_size}...")
            cov_matrix = renormalize_matrix(cov_matrix_original, block_size)

            # Recalculate correlation from renormalized covariance
            std = np.sqrt(np.diag(cov_matrix))
            cor_matrix = cov_matrix / np.outer(std, std)

        print(f"Matrix size: {cov_matrix.shape[0]} × {cov_matrix.shape[1]}")

        # Update variable names for renormalized matrices
        if block_size == 1:
            var_names_renorm = variable_names
        else:
            n_renorm = cov_matrix.shape[0]
            var_names_renorm = [f"block_{i}" for i in range(n_renorm)]

        # ====================================================================
        # PREPARE TEMPORARY FILE FOR C PROGRAM
        # ====================================================================

        if save_binary:
            cov_file_temp = os.path.join(k_dir, "covariance_temp.bin")
            save_matrix_binary(cov_matrix, cov_file_temp, dtype=np.float32)
            cor_file_temp = os.path.join(k_dir, "correlation_temp.bin")
            save_matrix_binary(cor_matrix, cor_file_temp, dtype=np.float32)
        else:
            cov_file_temp = os.path.join(k_dir, "covariance_temp.csv")
            save_matrix_for_c(cov_matrix, cov_file_temp)
            cor_file_temp = os.path.join(k_dir, "correlation_temp.csv")
            save_matrix_for_c(cor_matrix, cor_file_temp)

        # ====================================================================
        # RUN EIGENVALUE ANALYSIS
        # ====================================================================

        print("\n" + "-" * 70)
        if mode == "full":
            print("EIGENVALUE ANALYSIS (computing all eigenvalues)")
        else:
            print(f"EIGENVALUE ANALYSIS (computing K={K} largest eigenvalues)")
        print("-" * 70)

        results = run_eigenvalue_analysis_c(
            cov_file_temp,
            mode=mode,
            K=K,
            save_eigenvalues=save_eigenvalues,
            save_eigenvectors=save_eigenvectors,
        )

        # Remove temporary covariance files
        if os.path.exists(cov_file_temp):
            os.remove(cov_file_temp)
        if os.path.exists(cor_file_temp):
            os.remove(cor_file_temp)

        # ====================================================================
        # VARIANCE ANALYSIS
        # ====================================================================

        eigenvals = results["eigenvalues"]["eigenvalue"].values
        explained_var, idx_90, idx_95, idx_99 = compute_variance_explained(eigenvals)

        print(f"\nVariance Analysis:")
        if mode == "full":
            print(f"  - 90% variance with first {idx_90} eigenvectors")
            print(f"  - 95% variance with first {idx_95} eigenvectors")
            print(f"  - 99% variance with first {idx_99} eigenvectors")
        else:
            print(f"  - Total variance captured: {explained_var[-1]*100:.2f}%")
            if idx_90 <= K:
                print(f"  - 90% variance with first {idx_90} eigenvectors")
            if idx_95 <= K:
                print(f"  - 95% variance with first {idx_95} eigenvectors")
            if idx_99 <= K:
                print(f"  - 99% variance with first {idx_99} eigenvectors")

        # ====================================================================
        # POWER LAW FIT
        # ====================================================================

        n_skip = min(1, len(eigenvals) // 1)
        slope, r2, fit, indices, evals = fit_power_law(
            eigenvals, fit_fraction=0.6, n_skip=n_skip
        )

        print(f"\nPower Law Fit:")
        print(f"  - Exponent: {slope:.4f}")
        print(f"  - R²: {r2:.4f}")

        # ====================================================================
        # GENERATE PLOTS
        # ====================================================================

        figures = {}

        print(f"\nGenerating plots...")

        # Eigenvalue spectrum
        if plot_spectrum:
            suffix = "full" if mode == "full" else f"K{K}"
            output_file = (
                os.path.join(k_dir, f"spectrum_{suffix}.png") if save_plots else None
            )
            fig = plot_eigenspectrum(
                results, metadata, output_file=output_file, n_skip=n_skip
            )
            figures["spectrum"] = fig
            if not save_plots:
                plt.close(fig)

        # Eigenvectors
        if plot_eigvecs:
            suffix = "full" if mode == "full" else f"K{K}"
            output_file = (
                os.path.join(k_dir, f"eigenvectors_{suffix}.png")
                if save_plots
                else None
            )
            fig = plot_eigenvectors(
                results, metadata, n_vecs=n_eigvecs, output_file=output_file
            )
            figures["eigenvectors"] = fig
            if not save_plots:
                plt.close(fig)

        # Clustered covariance matrix (only for small matrices)
        if plot_cov_clustered:
            m = cov_matrix.shape[0]
            if m < umbral_JerClust:
                for xcolo in mis_colores:
                    output_file_cov = (
                        os.path.join(k_dir, f"covariance_clustered_{xcolo}.png")
                        if save_plots
                        else None
                    )
                    fig_cov = plot_covariance_clustered(
                        cov_matrix,
                        variable_names=var_names_renorm,
                        output_file=output_file_cov,
                        tipo_color=xcolo,
                    )
                    if not save_plots:
                        plt.close(fig_cov)

                    output_file_cor = (
                        os.path.join(k_dir, f"correlation_clustered_{xcolo}.png")
                        if save_plots
                        else None
                    )
                    fig_cor = plot_correlation_clustered(
                        cor_matrix,
                        variable_names=var_names_renorm,
                        output_file=output_file_cor,
                        tipo_color=xcolo,
                    )
                    if not save_plots:
                        plt.close(fig_cor)

                figures[f"covariance_clustered_{xcolo}"] = fig_cov
                figures[f"correlation_clustered_{xcolo}"] = fig_cor
            else:
                print(
                    f"  Skipping clustered matrix plots (matrix too large: {m} > {umbral_JerClust})"
                )

        # Dendrogram (only for small matrices)
        if plot_dendro:
            m = cov_matrix.shape[0]
            if m < umbral_JerClust:
                output_file = (
                    os.path.join(k_dir, "dendrogram.png") if save_plots else None
                )
                fig = plot_dendrogram(
                    cov_matrix, variable_names=var_names_renorm, output_file=output_file
                )
                figures["dendrogram"] = fig
                if not save_plots:
                    plt.close(fig)
            else:
                print(
                    f"  Skipping dendrogram (matrix too large: {m} > {umbral_JerClust})"
                )

        # ====================================================================
        # STORE RESULTS FOR THIS BLOCK SIZE
        # ====================================================================

        all_results[f"K{block_size:02d}"] = {
            "results": results,
            "metadata": metadata.copy(),
            "covariance": cov_matrix,
            "correlation": cor_matrix,
            "variable_names": var_names_renorm,
            "variance": {
                "explained": explained_var,
                "idx_90": idx_90,
                "idx_95": idx_95,
                "idx_99": idx_99,
            },
            "power_law": {"slope": slope, "r_squared": r2},
            "figures": figures,
            "block_size": block_size,
            "output_dir": k_dir,
        }

    # ========================================================================
    # SUMMARY
    # ========================================================================

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)

    print(f"\nProcessed {len(valid_sizes)} block size(s): {valid_sizes}")
    print(f"\nResults saved in: {output_dir}")

    for block_size in valid_sizes:
        k_dir = os.path.join(output_dir, f"K{block_size:02d}")
        print(f"\nK={block_size} outputs in: {k_dir}")

        if save_plots:
            if plot_spectrum:
                suffix = "full" if mode == "full" else f"K{K}"
                print(f"  - {os.path.join(k_dir, f'spectrum_{suffix}.png')}")
            if plot_eigvecs:
                suffix = "full" if mode == "full" else f"K{K}"
                print(f"  - {os.path.join(k_dir, f'eigenvectors_{suffix}.png')}")

            m = all_results[f"K{block_size:02d}"]["covariance"].shape[0]
            if plot_cov_clustered and m < umbral_JerClust:
                for xcolo in mis_colores:
                    print(
                        f"  - {os.path.join(k_dir, f'covariance_clustered_{xcolo}.png')}"
                    )
                    print(
                        f"  - {os.path.join(k_dir, f'correlation_clustered_{xcolo}.png')}"
                    )
            if plot_dendro and m < umbral_JerClust:
                print(f"  - {os.path.join(k_dir, 'dendrogram.png')}")

    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="Time Series Eigenvalue Analysis with Renormalization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full analysis with all options
  python timeseries_eigen_analysis.py data.csv --full --save-cov --save-plots
  
  # Partial analysis (first 100 eigenvalues)
  python timeseries_eigen_analysis.py data.csv --K 100 --save-plots
  
  # Analysis with renormalization (block sizes 1, 2, 3, 4)
  python timeseries_eigen_analysis.py data.csv --full --block-sizes "1,2,3,4"
  
  # Import existing covariance matrix with renormalization
  python timeseries_eigen_analysis.py --import-cov covariance.bin --full --block-sizes "1,2,4,8"
  
  # Import without variable names
  python timeseries_eigen_analysis.py --import-cov covariance.bin --no-names --K 50
  
  # Quick analysis without saving
  python timeseries_eigen_analysis.py data.csv --full --no-save-eigen --no-save-plots
  
  # Use CSV format for C communication (instead of binary)
  python timeseries_eigen_analysis.py data.csv --full --use-csv
  
  # Custom number of eigenvectors to plot
  python timeseries_eigen_analysis.py data.csv --full --n-eigvecs 5 --block-sizes "1,2"
        """,
    )

    # Input file (optional if importing)
    parser.add_argument(
        "input_file",
        nargs="?",
        type=str,
        help="Input CSV file with time series data (not needed if --import-cov is used)",
    )

    # Analysis mode
    analysis_group = parser.add_mutually_exclusive_group(required=True)
    analysis_group.add_argument(
        "--full", action="store_true", help="Compute all eigenvalues"
    )
    analysis_group.add_argument(
        "--K",
        type=int,
        help="Number of largest eigenvalues to compute (partial analysis)",
    )

    # Renormalization options
    parser.add_argument(
        "--block-sizes",
        type=str,
        default="1",
        help='Comma-separated list of block sizes for renormalization (default: "1", e.g., "1,2,3,4")',
    )

    # Import options
    parser.add_argument(
        "--import-cov",
        type=str,
        default=None,
        help="Import covariance matrix from binary file (.bin)",
    )
    parser.add_argument(
        "--import-cor",
        type=str,
        default=None,
        help="Import correlation matrix from binary file (.bin) [optional]",
    )
    parser.add_argument(
        "--no-names",
        action="store_true",
        help="Imported matrices do not have variable names file (.names)",
    )

    # Save options
    parser.add_argument(
        "--save-cov",
        action="store_true",
        help="Save covariance and correlation matrices (only for K=1)",
    )
    parser.add_argument(
        "--use-csv",
        action="store_true",
        help="Use CSV format instead of binary for C communication",
    )
    parser.add_argument(
        "--no-save-eigen",
        action="store_true",
        help="Do not save eigenvalues and eigenvectors",
    )
    parser.add_argument(
        "--no-save-plots", action="store_true", help="Do not save plots"
    )

    # Plot options
    parser.add_argument(
        "--no-plot-spectrum",
        action="store_true",
        help="Do not plot eigenvalue spectrum",
    )
    parser.add_argument(
        "--no-plot-eigvecs", action="store_true", help="Do not plot eigenvectors"
    )
    parser.add_argument(
        "--no-plot-cov",
        action="store_true",
        help="Do not plot clustered covariance matrix",
    )
    parser.add_argument(
        "--no-plot-dendro", action="store_true", help="Do not plot dendrogram"
    )
    parser.add_argument(
        "--n-eigvecs",
        type=int,
        default=3,
        help="Number of eigenvectors to plot (default: 3)",
    )

    # Output directory
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (default: same as input file or imported matrix)",
    )

    # Show plots
    parser.add_argument("--show", action="store_true", help="Show plots interactively")

    args = parser.parse_args()

    # ========================================================================
    # VALIDATE INPUTS
    # ========================================================================

    # Check if we're importing or computing
    importing = args.import_cov is not None

    if not importing and not args.input_file:
        parser.error("Either input_file or --import-cov must be provided")

    if not importing:
        # Check if input file exists
        if not os.path.exists(args.input_file):
            print(f"Error: Input file '{args.input_file}' not found.")
            sys.exit(1)
    else:
        # Check if import file exists
        if not os.path.exists(args.import_cov):
            print(f"Error: Covariance file '{args.import_cov}' not found.")
            sys.exit(1)

        # Check correlation file if provided
        if args.import_cor and not os.path.exists(args.import_cor):
            print(f"Error: Correlation file '{args.import_cor}' not found.")
            sys.exit(1)

    # Determine mode
    if args.full:
        mode = "full"
        K = None
    else:
        mode = "partial"
        K = args.K

    # Parse block sizes
    try:
        block_sizes = [int(k.strip()) for k in args.block_sizes.split(",")]
    except ValueError:
        print(f"Error: Invalid block sizes format: '{args.block_sizes}'")
        print("Block sizes must be comma-separated integers (e.g., '1,2,3,4')")
        sys.exit(1)

    # ========================================================================
    # RUN ANALYSIS
    # ========================================================================

    try:
        results = analyze_timeseries(
            input_file=args.input_file,
            mode=mode,
            K=K,
            import_covariance=args.import_cov,
            import_correlation=args.import_cor,
            has_names=not args.no_names,
            save_covariance=args.save_cov,
            save_binary=not args.use_csv,
            save_eigenvalues=not args.no_save_eigen,
            save_eigenvectors=not args.no_save_eigen,
            plot_spectrum=not args.no_plot_spectrum,
            plot_eigvecs=not args.no_plot_eigvecs,
            plot_cov_clustered=not args.no_plot_cov,
            plot_dendro=not args.no_plot_dendro,
            n_eigvecs=args.n_eigvecs,
            save_plots=not args.no_save_plots,
            output_dir=args.output_dir,
            block_sizes=block_sizes,
        )

        if args.show:
            plt.show()

    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
