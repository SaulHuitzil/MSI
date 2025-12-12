"""
Graphical User Interface for Time Series Eigenvalue Analysis
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import sys
import os
from pathlib import Path

from .cli import analyze_timeseries
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


class RedirectText:
    """Redirect stdout to a text widget"""

    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)
        self.text_widget.update_idletasks()

    def flush(self):
        pass


class TimeSeriesAnalysisGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Time Series Eigenvalue Analysis")
        self.root.geometry("900x1000")

        # Variables
        self.input_file = tk.StringVar()
        self.output_dir = tk.StringVar()
        self.mode = tk.StringVar(value="full")
        self.K_value = tk.IntVar(value=1000)

        # NEW: Import options
        self.import_covariance = tk.StringVar()
        self.has_variable_names = tk.BooleanVar(value=True)
        self.use_binary = tk.BooleanVar(value=True)

        # NEW: Renormalization options
        self.enable_renorm = tk.BooleanVar(value=False)
        self.block_sizes_str = tk.StringVar(value="1,2,3,4")

        # Save options
        self.save_covariance = tk.BooleanVar(value=True)
        self.save_eigenvalues = tk.BooleanVar(value=True)
        self.save_eigenvectors = tk.BooleanVar(value=True)
        self.save_plots = tk.BooleanVar(value=True)

        # Plot options
        self.plot_spectrum = tk.BooleanVar(value=True)
        self.plot_eigvecs = tk.BooleanVar(value=True)
        self.plot_cov_clustered = tk.BooleanVar(value=True)
        self.plot_dendro = tk.BooleanVar(value=True)
        self.n_eigvecs = tk.IntVar(value=15)
        self.show_plots = tk.BooleanVar(value=True)

        self.is_running = False

        self.setup_ui()

    def setup_ui(self):
        """Setup the user interface"""

        # Create main frame with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)

        current_row = 0

        # ====== FILE SELECTION ======
        file_frame = ttk.LabelFrame(main_frame, text="Input File", padding="10")
        file_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )
        file_frame.columnconfigure(1, weight=1)

        ttk.Label(file_frame, text="CSV File (optional if importing):").grid(
            row=0, column=0, sticky=tk.W, padx=5
        )
        ttk.Entry(file_frame, textvariable=self.input_file, width=50).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5
        )
        ttk.Button(file_frame, text="Browse...", command=self.browse_input_file).grid(
            row=0, column=2, padx=5
        )

        current_row += 1

        # ====== ANALYSIS MODE ======
        mode_frame = ttk.LabelFrame(main_frame, text="Analysis Mode", padding="10")
        mode_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )

        ttk.Radiobutton(
            mode_frame,
            text="Full Analysis (all eigenvalues)",
            variable=self.mode,
            value="full",
            command=self.update_mode,
        ).grid(row=0, column=0, sticky=tk.W, padx=5)

        partial_frame = ttk.Frame(mode_frame)
        partial_frame.grid(row=1, column=0, sticky=tk.W, padx=5)

        ttk.Radiobutton(
            partial_frame,
            text="Partial Analysis - K largest eigenvalues:",
            variable=self.mode,
            value="partial",
            command=self.update_mode,
        ).grid(row=0, column=0, sticky=tk.W)

        self.K_spinbox = ttk.Spinbox(
            partial_frame, from_=1, to=10000, textvariable=self.K_value, width=10
        )
        self.K_spinbox.grid(row=0, column=1, padx=10)

        current_row += 1

        # ====== RENORMALIZATION OPTIONS ======
        renorm_frame = ttk.LabelFrame(main_frame, text="Renormalization", padding="10")
        renorm_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )

        ttk.Checkbutton(
            renorm_frame,
            text="Enable renormalization (block averaging)",
            variable=self.enable_renorm,
            command=self.toggle_renorm,
        ).grid(row=0, column=0, columnspan=2, sticky=tk.W, padx=5, pady=2)

        ttk.Label(
            renorm_frame, text="Block sizes (comma-separated, e.g., 1,2,3,4):"
        ).grid(row=1, column=0, sticky=tk.W, padx=5, pady=(5, 0))

        self.renorm_entry = ttk.Entry(
            renorm_frame, textvariable=self.block_sizes_str, width=30
        )
        self.renorm_entry.grid(
            row=1, column=1, sticky=(tk.W, tk.E), padx=5, pady=(5, 0)
        )

        ttk.Label(
            renorm_frame,
            text="Each block size creates a separate analysis (K01, K02, etc.)",
            font=("TkDefaultFont", 8, "italic"),
        ).grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=5, pady=(2, 0))

        # Initially disable renorm entry
        self.renorm_entry.config(state="disabled")

        current_row += 1

        # ====== SAVE OPTIONS ======
        save_frame = ttk.LabelFrame(main_frame, text="Save Options", padding="10")
        save_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )

        ttk.Checkbutton(
            save_frame, text="Save covariance matrix", variable=self.save_covariance
        ).grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(
            save_frame, text="Save eigenvalues", variable=self.save_eigenvalues
        ).grid(row=0, column=1, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(
            save_frame, text="Save eigenvectors", variable=self.save_eigenvectors
        ).grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(save_frame, text="Save plots", variable=self.save_plots).grid(
            row=1, column=1, sticky=tk.W, padx=5, pady=2
        )
        ttk.Checkbutton(
            save_frame, text="Use binary format (.bin)", variable=self.use_binary
        ).grid(row=2, column=0, sticky=tk.W, padx=5, pady=2)

        current_row += 1

        # ====== PLOT OPTIONS ======
        plot_frame = ttk.LabelFrame(main_frame, text="Plot Options", padding="10")
        plot_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )

        ttk.Checkbutton(
            plot_frame, text="Plot eigenvalue spectrum", variable=self.plot_spectrum
        ).grid(row=0, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(
            plot_frame, text="Plot eigenvectors", variable=self.plot_eigvecs
        ).grid(row=0, column=1, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(
            plot_frame,
            text="Plot clustered covariance",
            variable=self.plot_cov_clustered,
        ).grid(row=1, column=0, sticky=tk.W, padx=5, pady=2)
        ttk.Checkbutton(
            plot_frame, text="Plot dendrogram", variable=self.plot_dendro
        ).grid(row=1, column=1, sticky=tk.W, padx=5, pady=2)

        eigvec_frame = ttk.Frame(plot_frame)
        eigvec_frame.grid(row=2, column=0, columnspan=2, sticky=tk.W, padx=5, pady=5)

        ttk.Label(eigvec_frame, text="Number of eigenvectors to plot:").grid(
            row=0, column=0, sticky=tk.W
        )
        ttk.Spinbox(
            eigvec_frame, from_=1, to=20, textvariable=self.n_eigvecs, width=10
        ).grid(row=0, column=1, padx=10)

        ttk.Checkbutton(
            plot_frame,
            text="Open output folder when finished",
            variable=self.show_plots,
        ).grid(row=3, column=0, columnspan=2, sticky=tk.W, padx=5, pady=2)

        current_row += 1

        # ====== OUTPUT DIRECTORY (Optional) ======
        output_frame = ttk.LabelFrame(
            main_frame, text="Output Directory (Optional)", padding="10"
        )
        output_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )
        output_frame.columnconfigure(1, weight=1)

        ttk.Label(output_frame, text="Custom output:").grid(
            row=0, column=0, sticky=tk.W, padx=5
        )
        ttk.Entry(output_frame, textvariable=self.output_dir, width=50).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5
        )
        ttk.Button(output_frame, text="Browse...", command=self.browse_output_dir).grid(
            row=0, column=2, padx=5
        )
        ttk.Label(
            output_frame,
            text="(Leave empty to create folder with CSV name)",
            font=("TkDefaultFont", 8, "italic"),
        ).grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5)

        current_row += 1

        # ====== IMPORT COVARIANCE (OPTIONAL) ======
        import_frame = ttk.LabelFrame(
            main_frame, text="Import Covariance Matrix (Optional)", padding="10"
        )
        import_frame.grid(
            row=current_row, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5
        )
        import_frame.columnconfigure(1, weight=1)

        ttk.Label(
            import_frame, text="Import covariance (.bin, optional if CSV provided):"
        ).grid(row=0, column=0, sticky=tk.W, padx=5)
        ttk.Entry(import_frame, textvariable=self.import_covariance, width=50).grid(
            row=0, column=1, sticky=(tk.W, tk.E), padx=5
        )
        ttk.Button(import_frame, text="Browse...", command=self.browse_import_cov).grid(
            row=0, column=2, padx=5
        )

        ttk.Checkbutton(
            import_frame,
            text="Imported matrix has variable names (.names file)",
            variable=self.has_variable_names,
        ).grid(row=1, column=0, columnspan=3, sticky=tk.W, padx=5, pady=2)

        current_row += 1

        # ====== RUN BUTTON ======
        self.run_button = ttk.Button(
            main_frame, text="Run Analysis", command=self.run_analysis
        )
        self.run_button.grid(
            row=current_row, column=0, columnspan=3, pady=10, sticky=(tk.W, tk.E)
        )

        current_row += 1

        # ====== LOG OUTPUT ======
        log_frame = ttk.LabelFrame(main_frame, text="Analysis Log", padding="10")
        log_frame.grid(
            row=current_row,
            column=0,
            columnspan=3,
            sticky=(tk.W, tk.E, tk.N, tk.S),
            pady=5,
        )
        log_frame.columnconfigure(0, weight=1)
        log_frame.rowconfigure(0, weight=1)

        self.log_text = scrolledtext.ScrolledText(log_frame, height=20, wrap=tk.WORD)
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure grid weights for resizing
        main_frame.rowconfigure(current_row, weight=1)

        # Initial state
        self.update_mode()

    def browse_input_file(self):
        """Browse for input CSV file"""
        filename = filedialog.askopenfilename(
            title="Select CSV file",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if filename:
            self.input_file.set(filename)

    def browse_output_dir(self):
        """Browse for output directory"""
        dirname = filedialog.askdirectory(title="Select output directory")
        if dirname:
            self.output_dir.set(dirname)

    def browse_import_cov(self):
        """Browse for covariance binary file"""
        filename = filedialog.askopenfilename(
            title="Select covariance binary file",
            filetypes=[("Binary files", "*.bin"), ("All files", "*.*")],
        )
        if filename:
            self.import_covariance.set(filename)

    def update_mode(self):
        """Update UI based on selected mode"""
        if self.mode.get() == "full":
            self.K_spinbox.config(state="disabled")
        else:
            self.K_spinbox.config(state="normal")

    def toggle_renorm(self):
        """Toggle renormalization entry state"""
        if self.enable_renorm.get():
            self.renorm_entry.config(state="normal")
        else:
            self.renorm_entry.config(state="disabled")

    def clear_log(self):
        """Clear the log text widget"""
        self.log_text.delete(1.0, tk.END)

    def validate_inputs(self):
        """Validate user inputs"""
        # If importing covariance, input_file is optional
        if self.import_covariance.get():
            if not os.path.exists(self.import_covariance.get()):
                messagebox.showerror(
                    "Error",
                    f"Covariance file does not exist:\n{self.import_covariance.get()}",
                )
                return False
            # When importing, we don't need input_file
        else:
            # If NOT importing, we need the input CSV file
            if not self.input_file.get():
                messagebox.showerror(
                    "Error",
                    "Please select an input CSV file or import a covariance matrix.",
                )
                return False

            if not os.path.exists(self.input_file.get()):
                messagebox.showerror(
                    "Error", f"Input file does not exist:\n{self.input_file.get()}"
                )
                return False

        if self.mode.get() == "partial" and self.K_value.get() < 1:
            messagebox.showerror("Error", "K must be at least 1 for partial analysis.")
            return False

        # Validate block sizes if renormalization is enabled
        if self.enable_renorm.get():
            try:
                block_sizes = [
                    int(k.strip()) for k in self.block_sizes_str.get().split(",")
                ]
                if not block_sizes:
                    messagebox.showerror(
                        "Error", "Please specify at least one block size."
                    )
                    return False
                if any(k < 1 for k in block_sizes):
                    messagebox.showerror("Error", "All block sizes must be >= 1.")
                    return False
            except ValueError:
                messagebox.showerror(
                    "Error",
                    "Invalid block sizes format.\nPlease use comma-separated integers (e.g., 1,2,3,4)",
                )
                return False

        return True

    def run_analysis(self):
        """Run the analysis in a separate thread"""
        if self.is_running:
            messagebox.showwarning("Warning", "Analysis is already running!")
            return

        if not self.validate_inputs():
            return

        # Disable run button
        self.run_button.config(state="disabled")
        self.is_running = True

        # Clear log
        self.clear_log()

        # Redirect stdout to log
        sys.stdout = RedirectText(self.log_text)
        sys.stderr = RedirectText(self.log_text)

        # Run in thread to keep GUI responsive
        thread = threading.Thread(target=self.run_analysis_thread)
        thread.daemon = True
        thread.start()

    def run_analysis_thread(self):
        """Thread function to run analysis"""
        try:
            # Prepare parameters
            mode = self.mode.get()
            K = self.K_value.get() if mode == "partial" else None
            output_dir = self.output_dir.get() if self.output_dir.get() else None

            import_cov = (
                self.import_covariance.get() if self.import_covariance.get() else None
            )

            # Input file is optional when importing
            input_file = self.input_file.get() if self.input_file.get() else None

            # Get block sizes
            if self.enable_renorm.get():
                block_sizes = [
                    int(k.strip()) for k in self.block_sizes_str.get().split(",")
                ]
            else:
                block_sizes = [1]  # Default: no renormalization

            # Run analysis (plots are saved but not shown)
            results = analyze_timeseries(
                input_file=input_file,  # Can be None
                mode=mode,
                K=K,
                import_covariance=import_cov,
                import_correlation=None,  # Could add GUI control for this
                has_names=self.has_variable_names.get(),
                save_covariance=self.save_covariance.get(),
                save_binary=self.use_binary.get(),
                save_eigenvalues=self.save_eigenvalues.get(),
                save_eigenvectors=self.save_eigenvectors.get(),
                plot_spectrum=self.plot_spectrum.get(),
                plot_eigvecs=self.plot_eigvecs.get(),
                plot_cov_clustered=self.plot_cov_clustered.get(),
                plot_dendro=self.plot_dendro.get(),
                n_eigvecs=self.n_eigvecs.get(),
                save_plots=self.save_plots.get(),
                output_dir=output_dir,
                block_sizes=block_sizes,
            )

            # Close all figures to free memory
            plt.close("all")

            # Determine output directory for success message
            if output_dir:
                result_dir = output_dir
            elif import_cov:
                # When importing, use the directory of the imported file
                base_name = os.path.splitext(os.path.basename(import_cov))[0]
                base_name = base_name.replace("_covariance", "")
                parent_dir = os.path.dirname(import_cov) or "."
                result_dir = os.path.join(parent_dir, f"{base_name}_analysis")
            else:
                # When computing, use the directory of the input file
                base_name = os.path.splitext(os.path.basename(input_file))[0]
                parent_dir = os.path.dirname(input_file) or "."
                result_dir = os.path.join(parent_dir, f"{base_name}_analysis")

            # If show_plots is enabled, open the output folder
            if self.show_plots.get():
                self.root.after(0, lambda: self.open_output_folder(result_dir))

        except Exception as e:
            error_msg = f"Error during analysis:\n{str(e)}"
            print(error_msg)
            import traceback

            traceback.print_exc()

            self.root.after(0, lambda: messagebox.showerror("Error", error_msg))

        finally:
            # Reset stdout
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

            # Re-enable run button
            self.root.after(0, lambda: self.run_button.config(state="normal"))
            self.is_running = False

    def open_output_folder(self, folder_path):
        """Open the output folder in the system file browser"""
        import platform
        import subprocess

        try:
            system = platform.system()
            if system == "Darwin":  # macOS
                subprocess.run(["open", folder_path])
            elif system == "Windows":
                subprocess.run(["explorer", folder_path])
            else:  # Linux
                subprocess.run(["xdg-open", folder_path])
        except Exception as e:
            print(f"Could not open folder: {e}")


def launch_gui():
    """Launch the GUI application"""
    root = tk.Tk()

    # Set style
    style = ttk.Style()
    style.theme_use("clam")  # Use 'clam', 'alt', 'default', or 'classic'

    app = TimeSeriesAnalysisGUI(root)
    root.mainloop()


if __name__ == "__main__":
    launch_gui()
