#!/usr/bin/env python3
"""
pysca.notebook_utils

Python API utilities for using pySCA 7.0 in Jupyter notebooks.

This module provides convenient Python functions that wrap the command-line
tools (sca-process-msa and sca-core) for programmatic use in notebooks.
It also provides visualization helpers and data loading utilities.

Example:
    >>> from pysca import notebook_utils as nb
    >>> 
    >>> # Process an MSA
    >>> db_path = nb.process_msa(
    ...     "Inputs/alignment.fasta",
    ...     pdb="Inputs/1ABC.pdb",
    ...     chain="A",
    ...     species="Homo sapiens"
    ... )
    >>> 
    >>> # Run SCA core calculations
    >>> nb.run_sca_core(
    ...     db_path,
    ...     do_sector_id=True,
    ...     kpos=6
    ... )
    >>> 
    >>> # Load and visualize results
    >>> db = nb.load_database(db_path)
    >>> nb.plot_eigenvalues(db)
    >>> nb.plot_ic_heatmap(db)
"""

from __future__ import annotations

import gzip
import logging
import pickle
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, List

import numpy as np

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from bokeh.io import output_notebook, show
    from bokeh.plotting import figure
    from bokeh.models import ColorBar, LinearColorMapper, HoverTool
    from bokeh.palettes import Viridis256
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Import main functions from the command-line scripts
try:
    from pysca.scaProcessMSA_py3_big import main as sca_process_msa_main
    from pysca.scaCore_py3 import main as sca_core_main
except ImportError:
    sca_process_msa_main = None
    sca_core_main = None

# Import core tools
from pysca import scaTools as sca


def load_database(db_path: str | Path) -> Dict[str, Any]:
    """
    Load a pySCA database file (.db or .db.gz).
    
    Parameters:
        db_path: Path to database file
        
    Returns:
        Dictionary containing 'sequence', 'sca', and optionally 'sector' keys
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> print(db.keys())
        dict_keys(['sequence', 'sca', 'sector'])
    """
    db_path = Path(db_path)
    if db_path.suffix == ".gz" or db_path.name.endswith(".db.gz"):
        with gzip.open(db_path, "rb") as f:
            return pickle.load(f)
    with open(db_path, "rb") as f:
        return pickle.load(f)


def process_msa(
    alignment: str | Path,
    output_dir: str | Path = "Outputs",
    pdb: Optional[str | Path] = None,
    chain: str = "A",
    species: Optional[str] = None,
    refseq: Optional[str | Path] = None,
    refpos: Optional[str | Path] = None,
    i_ref: Optional[int] = None,
    parameters: Optional[List[float]] = None,
    precluster: Optional[bool] = None,
    cluster_id: float = 0.85,
    keep_sequences: Optional[List[int]] = None,
    keep_sequences_file: Optional[str | Path] = None,
    initial_trim_gap: float = 0.8,
    matlab: bool = False,
    verbose: bool = True,
) -> Path:
    """
    Process a multiple sequence alignment using sca-process-msa.
    
    This is a Python wrapper around the sca-process-msa command-line tool.
    
    Parameters:
        alignment: Path to input alignment (FASTA/Stockholm/Clustal; optionally .gz)
        output_dir: Directory for output files (default: "Outputs")
        pdb: PDB identifier or path for reference structure
        chain: Chain ID in PDB (default: "A")
        species: Species name for reference sequence search
        refseq: Path to reference sequence FASTA file (alternative to PDB)
        refpos: Path to reference positions file (alternative to PDB)
        i_ref: Reference sequence index (0-based, alternative to PDB/refseq)
        parameters: List of 4 floats [pos_gap, seq_gap, min_seqid, max_seqid]
                   (default: [0.25, 0.2, 0.15, 0.85])
        precluster: Enable/disable MMseqs2 preclustering (None = auto for >50k seqs)
        cluster_id: MMseqs2 identity threshold (default: 0.85)
        keep_sequences: List of 1-based sequence indices to retain during preclustering
        keep_sequences_file: File with 1-based sequence indices (one per line)
        initial_trim_gap: Gap fraction threshold for initial position trimming (default: 0.8)
        matlab: Write MATLAB workspace file (default: False)
        verbose: Enable verbose logging (default: True)
        
    Returns:
        Path to the created database file
        
    Example:
        >>> db_path = process_msa(
        ...     "Inputs/alignment.fasta",
        ...     pdb="Inputs/1ABC.pdb",
        ...     chain="A",
        ...     species="Homo sapiens",
        ...     precluster=True
        ... )
        >>> print(f"Database saved to: {db_path}")
    """
    if sca_process_msa_main is None:
        raise ImportError("scaProcessMSA_py3_big module not available")
    
    # Build argument list
    args = [str(alignment)]
    
    if pdb:
        args.extend(["-s", str(pdb)])
        args.extend(["-c", chain])
    if species:
        args.extend(["-f", species])
    if refseq:
        args.extend(["-r", str(refseq)])
    if refpos:
        args.extend(["-o", str(refpos)])
    if i_ref is not None:
        args.extend(["-i", str(i_ref)])
    
    if parameters:
        args.extend(["-p"] + [str(p) for p in parameters])
    
    if precluster is True:
        args.append("--precluster")
    elif precluster is False:
        args.append("--no-precluster")
    
    if cluster_id != 0.85:
        args.extend(["--cluster-id", str(cluster_id)])
    
    if keep_sequences:
        for idx in keep_sequences:
            args.extend(["--keep-sequences", str(idx)])
    
    if keep_sequences_file:
        args.extend(["--keep-sequences-file", str(keep_sequences_file)])
    
    if initial_trim_gap != 0.8:
        args.extend(["--initial-trim-gap", str(initial_trim_gap)])
    
    if matlab:
        args.append("-m")
    
    if not verbose:
        args.append("--quiet")
    
    # Set output directory via --output
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base_name = Path(alignment).stem.replace(".fasta", "").replace(".fa", "")
    output_file = output_dir / f"{base_name}.db.gz"
    args.extend(["--output", str(output_file)])
    
    # Run the main function
    # Temporarily replace sys.argv to pass arguments
    import sys
    original_argv = sys.argv
    try:
        sys.argv = ["sca-process-msa"] + args
        exit_code = sca_process_msa_main()
    finally:
        sys.argv = original_argv
    
    if exit_code != 0:
        raise RuntimeError(f"sca-process-msa failed with exit code {exit_code}")
    
    return output_file


def run_sca_core(
    database: str | Path,
    output_db: Optional[str | Path] = None,
    norm: str = "frob",
    Ntrials: int = 10,
    lbda: float = 0.01,
    do_seqcorr: bool = False,
    do_sector_id: bool = False,
    kpos: Optional[int] = None,
    kica: Optional[int] = None,
    ic_cutoff: float = 0.95,
    kmax_cap: int = 10,
    max_seqcorr_seqs: Optional[int] = None,
    seqcorr_ref: Optional[int] = None,
    seqcorr_mmseqs2: bool = False,
    seqcorr_keep_indices: Optional[List[int]] = None,
    seqcorr_max_cap: int = 50000,
    auto_seqcorr_subsample: bool = True,
    float32: bool = False,
    matlab: bool = False,
    no_db: bool = False,
    verbose: bool = True,
) -> Optional[Path]:
    """
    Run SCA core calculations using sca-core.
    
    This is a Python wrapper around the sca-core command-line tool.
    
    Parameters:
        database: Path to database from sca-process-msa
        output_db: Output database path (default: same as input with _core suffix)
        norm: Norm type for SCA matrix reduction: 'frob' or 'spec' (default: 'frob')
        Ntrials: Number of randomization trials (default: 10)
        lbda: Regularization parameter (default: 0.01)
        do_seqcorr: Compute sequence correlation matrices (simMat/Useq/Uica)
        do_sector_id: Perform independent component identification
        kpos: Number of eigenmodes to keep (None = auto)
        kica: Number of independent components for ICA (None = use kpos)
        ic_cutoff: T-distribution cutoff percentile for IC identification (default: 0.95)
        kmax_cap: Maximum kpos when auto-selecting (default: 10)
        max_seqcorr_seqs: Maximum sequences for seqSim (None = auto)
        seqcorr_ref: Reference sequence index for seqSim (None = use i_ref from DB)
        seqcorr_mmseqs2: Use MMseqs2 for seqSim subsampling
        seqcorr_keep_indices: Additional 1-based sequence indices to retain in seqSim
        seqcorr_max_cap: Maximum sequences cap for seqSim (default: 50000)
        auto_seqcorr_subsample: Enable automatic 1.5×M_eff subsampling (default: True)
        float32: Use float32 for memory efficiency (default: False)
        matlab: Write MATLAB workspace file (default: False)
        no_db: Suppress database writing (default: False)
        verbose: Enable verbose logging (default: True)
        
    Returns:
        Path to output database file (None if no_db=True)
        
    Example:
        >>> db_path = run_sca_core(
        ...     "Outputs/alignment.db.gz",
        ...     do_sector_id=True,
        ...     kpos=6
        ... )
    """
    if sca_core_main is None:
        raise ImportError("scaCore_py3 module not available")
    
    # Build argument list
    args = [str(database)]
    
    if output_db:
        args.extend(["-o", str(output_db)])
    
    if norm != "frob":
        args.extend(["-n", norm])
    
    if Ntrials != 10:
        args.extend(["-t", str(Ntrials)])
    
    if lbda != 0.01:
        args.extend(["-l", str(lbda)])
    
    if do_seqcorr:
        args.append("--do-seqcorr")
    
    if do_sector_id:
        args.append("--do-sector-id")
    
    if kpos is not None:
        args.extend(["--kpos", str(kpos)])
    
    if kica is not None:
        args.extend(["--kica", str(kica)])
    
    if ic_cutoff != 0.95:
        args.extend(["--ic-cutoff", str(ic_cutoff)])
    
    if kmax_cap != 10:
        args.extend(["--kmax-cap", str(kmax_cap)])
    
    if max_seqcorr_seqs is not None:
        args.extend(["--max-seqcorr-seqs", str(max_seqcorr_seqs)])
    
    if seqcorr_ref is not None:
        args.extend(["--seqcorr-ref", str(seqcorr_ref)])
    
    if seqcorr_mmseqs2:
        args.append("--seqcorr-mmseqs2")
    
    if seqcorr_keep_indices:
        for idx in seqcorr_keep_indices:
            args.extend(["--seqcorr-keep-indices", str(idx)])
    
    if seqcorr_max_cap != 50000:
        args.extend(["--seqcorr-max-cap", str(seqcorr_max_cap)])
    
    if not auto_seqcorr_subsample:
        args.append("--no-auto-seqcorr-subsample")
    
    if float32:
        args.append("--float32")
    
    if matlab:
        args.append("-m")
    
    if no_db:
        args.append("--no-db")
    
    if not verbose:
        args.append("--quiet")
    
    # Run the main function
    # Temporarily replace sys.argv to pass arguments
    import sys
    original_argv = sys.argv
    try:
        sys.argv = ["sca-core"] + args
        exit_code = sca_core_main()
    finally:
        sys.argv = original_argv
    
    if exit_code != 0:
        raise RuntimeError(f"sca-core failed with exit code {exit_code}")
    
    if no_db:
        return None
    
    # Determine output path
    if output_db:
        return Path(output_db)
    else:
        # Default: same as input (sca-core updates in place)
        return Path(database)


# ============================================================================
# Visualization Functions
# ============================================================================

def plot_eigenvalues(
    db: Dict[str, Any],
    n_modes: Optional[int] = None,
    show_randomized: bool = True,
    figsize: Tuple[int, int] = (10, 6),
    ax=None,
) -> None:
    """
    Plot eigenvalues of the SCA matrix.
    
    Parameters:
        db: Database dictionary from load_database()
        n_modes: Number of modes to plot (None = all)
        show_randomized: Overlay randomized trial eigenvalues (default: True)
        figsize: Figure size (width, height)
        ax: Matplotlib axes (None = create new figure)
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> plot_eigenvalues(db, n_modes=20)
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for plotting")
    
    if "sca" not in db:
        raise ValueError("Database must contain 'sca' key")
    
    sca_data = db["sca"]
    if "L" not in sca_data:
        raise ValueError("Database must contain 'sca.L' (eigenvalues)")
    
    L = sca_data["L"]
    if n_modes is None:
        n_modes = len(L)
    else:
        n_modes = min(n_modes, len(L))
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    modes = np.arange(1, n_modes + 1)
    ax.plot(modes, L[:n_modes], "o-", label="SCA eigenvalues", linewidth=2, markersize=6)
    
    if show_randomized and "Lrand" in sca_data:
        Lrand = sca_data["Lrand"]
        n_trials = Lrand.shape[0]
        Lrand_mean = np.mean(Lrand[:, :n_modes], axis=0)
        Lrand_std = np.std(Lrand[:, :n_modes], axis=0)
        
        ax.errorbar(
            modes,
            Lrand_mean,
            yerr=Lrand_std,
            fmt="s-",
            label=f"Randomized (mean ± std, N={n_trials})",
            alpha=0.7,
            markersize=4,
        )
    
    ax.set_xlabel("Eigenmode", fontsize=12)
    ax.set_ylabel("Eigenvalue", fontsize=12)
    ax.set_title("SCA Eigenvalue Spectrum", fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if ax is None:
        plt.tight_layout()
        plt.show()


def plot_ic_heatmap(
    db: Dict[str, Any],
    ic_index: Optional[int] = None,
    figsize: Tuple[int, int] = (12, 8),
    cmap: str = "RdBu_r",
    ax=None,
) -> None:
    """
    Plot heatmap of independent component values across positions.
    
    Parameters:
        db: Database dictionary from load_database()
        ic_index: Which IC to plot (None = plot all ICs)
        figsize: Figure size (width, height)
        cmap: Colormap name (default: "RdBu_r")
        ax: Matplotlib axes (None = create new figure)
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> plot_ic_heatmap(db, ic_index=0)
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for plotting")
    
    if "sector" not in db:
        raise ValueError("Database must contain 'sector' key (run with --do-sector-id)")
    
    sector_data = db["sector"]
    if "Vica" not in sector_data:
        raise ValueError("Database must contain 'sector.Vica' (IC vectors)")
    
    Vica = sector_data["Vica"]
    n_ics, n_pos = Vica.shape
    
    if ic_index is not None:
        if ic_index >= n_ics:
            raise ValueError(f"IC index {ic_index} out of range (0-{n_ics-1})")
        Vica = Vica[ic_index:ic_index+1, :]
        n_ics = 1
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(Vica, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_xlabel("Position (ATS)", fontsize=12)
    ax.set_ylabel("Independent Component", fontsize=12)
    ax.set_title("Independent Component Values", fontsize=14)
    
    if ic_index is None:
        ax.set_yticks(range(n_ics))
        ax.set_yticklabels([f"IC {i+1}" for i in range(n_ics)])
    else:
        ax.set_yticks([0])
        ax.set_yticklabels([f"IC {ic_index+1}"])
    
    plt.colorbar(im, ax=ax, label="IC Value")
    
    if ax is None:
        plt.tight_layout()
        plt.show()


def plot_ic_positions(
    db: Dict[str, Any],
    ic_index: int,
    top_n: Optional[int] = None,
    figsize: Tuple[int, int] = (12, 6),
    ax=None,
) -> None:
    """
    Plot significant positions for a specific independent component.
    
    Parameters:
        db: Database dictionary from load_database()
        ic_index: Which IC to plot (0-based)
        top_n: Number of top positions to highlight (None = use all significant)
        figsize: Figure size (width, height)
        ax: Matplotlib axes (None = create new figure)
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> plot_ic_positions(db, ic_index=0, top_n=20)
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for plotting")
    
    if "sector" not in db:
        raise ValueError("Database must contain 'sector' key")
    
    sector_data = db["sector"]
    if "Vica" not in sector_data:
        raise ValueError("Database must contain 'sector.Vica'")
    
    Vica = sector_data["Vica"]
    if ic_index >= Vica.shape[0]:
        raise ValueError(f"IC index {ic_index} out of range")
    
    ic_values = Vica[ic_index, :]
    positions = np.arange(len(ic_values))
    
    # Get significant positions if available
    significant_pos = None
    if "ic_pos" in sector_data and ic_index < len(sector_data["ic_pos"]):
        significant_pos = set(sector_data["ic_pos"][ic_index])
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    
    # Plot all positions
    ax.plot(positions, ic_values, "o-", alpha=0.3, markersize=3, label="All positions")
    
    # Highlight significant positions
    if significant_pos:
        sig_indices = [i for i in positions if i in significant_pos]
        sig_values = [ic_values[i] for i in sig_indices]
        ax.scatter(sig_indices, sig_values, c="red", s=50, zorder=5, label="Significant positions")
    
    # Highlight top N positions by absolute value
    if top_n:
        abs_values = np.abs(ic_values)
        top_indices = np.argsort(abs_values)[-top_n:][::-1]
        top_values = ic_values[top_indices]
        ax.scatter(top_indices, top_values, c="orange", s=30, zorder=4, alpha=0.7, label=f"Top {top_n} positions")
    
    ax.set_xlabel("Position (ATS)", fontsize=12)
    ax.set_ylabel("IC Value", fontsize=12)
    ax.set_title(f"Independent Component {ic_index+1}", fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    
    if ax is None:
        plt.tight_layout()
        plt.show()


def plot_simmat_interactive(
    db: Dict[str, Any],
    output_notebook_mode: bool = True,
    width: int = 800,
    height: int = 800,
    colormap: str = "Viridis256",
) -> None:
    """
    Create interactive Bokeh heatmap of sequence similarity matrix.
    
    Parameters:
        db: Database dictionary from load_database()
        output_notebook_mode: Use output_notebook() for inline display (default: True)
        width: Plot width in pixels
        height: Plot height in pixels
        colormap: Bokeh palette name
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> plot_simmat_interactive(db)
    """
    if not BOKEH_AVAILABLE:
        raise ImportError("bokeh is required for interactive plots. Install with: pip install bokeh")
    
    if "sca" not in db:
        raise ValueError("Database must contain 'sca' key")
    
    sca_data = db["sca"]
    if "simMat" not in sca_data:
        raise ValueError("Database must contain 'sca.simMat' (run with --do-seqcorr)")
    
    # Import visualization function from scaVisualizeSimMat
    from pysca.scaVisualizeSimMat import create_heatmap, extract_simmat
    
    simMat, headers, selected_indices = extract_simmat(db)
    
    if output_notebook_mode:
        output_notebook()
    
    create_heatmap(
        simMat=simMat,
        headers=headers,
        output_path=None,
        title="Sequence Similarity Matrix",
        colormap=colormap,
        width=width,
        height=height,
        show_plot=True,
    )


# ============================================================================
# Plotly Visualization Functions
# ============================================================================

def plot_eigenvalues_plotly(
    db: Dict[str, Any],
    show_randomized: bool = True,
    height: int = 500,
    width: int = 800,
    show_plot: bool = True,
) -> go.Figure:
    """
    Plot histogram of eigenvalues with randomized distribution overlaid.
    
    Uses all eigenvalues and creates a histogram with Npos bins.
    
    Parameters:
        db: Database dictionary from load_database()
        show_randomized: Overlay randomized trial eigenvalue distribution (default: True)
        height: Plot height in pixels
        width: Plot width in pixels
        show_plot: Whether to display the plot (default: True)
        
    Returns:
        plotly.graph_objects.Figure
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> fig = plot_eigenvalues_plotly(db)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for plotting. Install with: pip install plotly")
    
    if "sca" not in db:
        raise ValueError("Database must contain 'sca' key")
    
    sca_data = db["sca"]
    if "L" not in sca_data:
        raise ValueError("Database must contain 'sca.L' (eigenvalues)")
    
    # Get number of positions for histogram bins
    if "sequence" not in db:
        raise ValueError("Database must contain 'sequence' key")
    n_pos = db["sequence"].get("Npos")
    if n_pos is None:
        raise ValueError("Database must contain 'sequence.Npos' (number of positions)")
    n_pos = int(n_pos)
    
    # Use all eigenvalues
    L = sca_data["L"]
    
    fig = go.Figure()
    
    # Create histogram of SCA eigenvalues with Npos bins
    fig.add_trace(go.Histogram(
        x=L,
        nbinsx=n_pos,
        name='SCA eigenvalues',
        marker=dict(color='#1f77b4', line=dict(color='darkblue', width=1)),
        opacity=0.7,
        hovertemplate='Eigenvalue range: %{x}<br>Count: %{y}<extra></extra>'
    ))
    
    # Add randomized eigenvalue distribution as line plot if available
    if show_randomized and "Lrand" in sca_data:
        Lrand = sca_data["Lrand"]
        n_trials = Lrand.shape[0]
        
        # Flatten all randomized eigenvalues across all trials and modes
        Lrand_flat = Lrand.flatten()
        
        # Create histogram of randomized eigenvalues with same bins
        # Use the same bin edges as the SCA histogram
        L_min = min(L.min(), Lrand_flat.min())
        L_max = max(L.max(), Lrand_flat.max())
        bin_edges = np.linspace(L_min, L_max, n_pos + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Compute histogram of randomized eigenvalues
        hist_rand, _ = np.histogram(Lrand_flat, bins=bin_edges)
        
        # Normalize to match the scale (optional - you might want to adjust this)
        # Scale by the ratio of total counts to make them comparable
        scale_factor = len(L) / (n_trials * Lrand.shape[1])
        hist_rand_scaled = hist_rand * scale_factor
        
        # Add as line plot
        fig.add_trace(go.Scatter(
            x=bin_centers,
            y=hist_rand_scaled,
            mode='lines',
            name=f'Randomized (N={n_trials} trials)',
            line=dict(width=2, color='#ff7f0e', dash='dash'),
            hovertemplate='Eigenvalue: %{x:.4f}<br>Count: %{y:.2f}<extra></extra>'
        ))
    
    fig.update_layout(
        title='SCA Eigenvalue Spectrum',
        xaxis_title='Eigenvalue',
        yaxis_title='Count',
        height=height,
        width=width,
        hovermode='x unified',
        template='plotly_white',
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ),
        bargap=0.1  # Gap between bars in histogram
    )
    
    if show_plot:
        # For Cursor/VS Code, use HTML display since mimetype isn't supported
        # This works in both Jupyter Lab (which supports mimetype) and Cursor (which doesn't)
        try:
            from IPython.display import HTML, display
            # Use HTML display which works in both Jupyter and Cursor/VS Code
            html_str = fig.to_html(include_plotlyjs='cdn', div_id=f"plotly-{id(fig)}")
            display(HTML(html_str))
        except ImportError:
            # If IPython not available, try regular show (will open in browser)
            fig.show()
        except Exception:
            # Final fallback: regular show
            fig.show()
    
    return fig


def plot_ic_heatmap_plotly(
    db: Dict[str, Any],
    ic_index: Optional[int] = None,
    height: int = 600,
    width: int = 1000,
    colorscale: str = "RdBu",
    show_plot: bool = True,
) -> go.Figure:
    """
    Plot heatmap of independent component values across positions using Plotly.
    
    Parameters:
        db: Database dictionary from load_database()
        ic_index: Which IC to plot (None = plot all ICs)
        height: Plot height in pixels
        width: Plot width in pixels
        colorscale: Plotly colorscale name (default: "RdBu")
        show_plot: Whether to display the plot (default: True)
        
    Returns:
        plotly.graph_objects.Figure
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> fig = plot_ic_heatmap_plotly(db)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for plotting. Install with: pip install plotly")
    
    if "sector" not in db:
        raise ValueError("Database must contain 'sector' key (run with --do-sector-id)")
    
    sector_data = db["sector"]
    if "Vica" not in sector_data:
        raise ValueError("Database must contain 'sector.Vica' (IC vectors)")
    
    Vica = sector_data["Vica"]
    n_ics, n_pos = Vica.shape
    
    if ic_index is not None:
        if ic_index >= n_ics:
            raise ValueError(f"IC index {ic_index} out of range (0-{n_ics-1})")
        Vica = Vica[ic_index:ic_index+1, :]
        n_ics = 1
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=Vica,
        x=list(range(n_pos)),
        y=[f"IC {i+1}" for i in range(n_ics)] if ic_index is None else [f"IC {ic_index+1}"],
        colorscale=colorscale,
        colorbar=dict(title="IC Value"),
        hovertemplate='Position: %{x}<br>IC: %{y}<br>Value: %{z:.4f}<extra></extra>'
    ))
    
    fig.update_layout(
        title='Independent Component Values',
        xaxis_title='Position (ATS)',
        yaxis_title='Independent Component',
        height=height,
        width=width,
        template='plotly_white'
    )
    
    if show_plot:
        # For Cursor/VS Code, use HTML display since mimetype isn't supported
        # This works in both Jupyter Lab (which supports mimetype) and Cursor (which doesn't)
        try:
            from IPython.display import HTML, display
            # Use HTML display which works in both Jupyter and Cursor/VS Code
            html_str = fig.to_html(include_plotlyjs='cdn', div_id=f"plotly-{id(fig)}")
            display(HTML(html_str))
        except ImportError:
            # If IPython not available, try regular show (will open in browser)
            fig.show()
        except Exception:
            # Final fallback: regular show
            fig.show()
    
    return fig


def plot_ic_positions_plotly(
    db: Dict[str, Any],
    ic_index: int,
    top_n: Optional[int] = None,
    height: int = 500,
    width: int = 1000,
    show_plot: bool = True,
) -> go.Figure:
    """
    Plot significant positions for a specific independent component using Plotly.
    
    Parameters:
        db: Database dictionary from load_database()
        ic_index: Which IC to plot (0-based)
        top_n: Number of top positions to highlight (None = use all significant)
        height: Plot height in pixels
        width: Plot width in pixels
        show_plot: Whether to display the plot (default: True)
        
    Returns:
        plotly.graph_objects.Figure
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> fig = plot_ic_positions_plotly(db, ic_index=0, top_n=20)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for plotting. Install with: pip install plotly")
    
    if "sector" not in db:
        raise ValueError("Database must contain 'sector' key")
    
    sector_data = db["sector"]
    if "Vica" not in sector_data:
        raise ValueError("Database must contain 'sector.Vica'")
    
    Vica = sector_data["Vica"]
    if ic_index >= Vica.shape[0]:
        raise ValueError(f"IC index {ic_index} out of range")
    
    ic_values = Vica[ic_index, :]
    positions = np.arange(len(ic_values))
    
    # Get significant positions if available
    significant_pos = None
    if "ic_pos" in sector_data and ic_index < len(sector_data["ic_pos"]):
        significant_pos = set(sector_data["ic_pos"][ic_index])
    
    fig = go.Figure()
    
    # Plot all positions
    fig.add_trace(go.Scatter(
        x=positions,
        y=ic_values,
        mode='lines+markers',
        name='All positions',
        line=dict(width=1, color='lightgray'),
        marker=dict(size=3, color='lightgray', opacity=0.5),
        hovertemplate='Position: %{x}<br>IC Value: %{y:.4f}<extra></extra>'
    ))
    
    # Highlight significant positions
    if significant_pos:
        sig_indices = [i for i in positions if i in significant_pos]
        sig_values = [ic_values[i] for i in sig_indices]
        fig.add_trace(go.Scatter(
            x=sig_indices,
            y=sig_values,
            mode='markers',
            name='Significant positions',
            marker=dict(size=8, color='red', line=dict(width=1, color='darkred')),
            hovertemplate='Position: %{x}<br>IC Value: %{y:.4f}<extra></extra>'
        ))
    
    # Highlight top N positions by absolute value
    if top_n:
        abs_values = np.abs(ic_values)
        top_indices = np.argsort(abs_values)[-top_n:][::-1]
        top_values = ic_values[top_indices]
        fig.add_trace(go.Scatter(
            x=top_indices,
            y=top_values,
            mode='markers',
            name=f'Top {top_n} positions',
            marker=dict(size=6, color='orange', line=dict(width=1, color='darkorange'), opacity=0.7),
            hovertemplate='Position: %{x}<br>IC Value: %{y:.4f}<extra></extra>'
        ))
    
    fig.update_layout(
        title=f'Independent Component {ic_index+1}',
        xaxis_title='Position (ATS)',
        yaxis_title='IC Value',
        height=height,
        width=width,
        hovermode='x unified',
        template='plotly_white',
        shapes=[dict(
            type='line',
            x0=0,
            x1=len(positions)-1,
            y0=0,
            y1=0,
            line=dict(color='black', width=1, dash='dash'),
            opacity=0.3
        )]
    )
    
    if show_plot:
        # For Cursor/VS Code, use HTML display since mimetype isn't supported
        # This works in both Jupyter Lab (which supports mimetype) and Cursor (which doesn't)
        try:
            from IPython.display import HTML, display
            # Use HTML display which works in both Jupyter and Cursor/VS Code
            html_str = fig.to_html(include_plotlyjs='cdn', div_id=f"plotly-{id(fig)}")
            display(HTML(html_str))
        except ImportError:
            # If IPython not available, try regular show (will open in browser)
            fig.show()
        except Exception:
            # Final fallback: regular show
            fig.show()
    
    return fig


def plot_simmat_plotly(
    db: Dict[str, Any],
    height: int = 800,
    width: int = 800,
    colorscale: str = "Viridis",
    show_plot: bool = True,
) -> go.Figure:
    """
    Create interactive Plotly heatmap of sequence similarity matrix.
    
    Parameters:
        db: Database dictionary from load_database()
        height: Plot height in pixels
        width: Plot width in pixels
        colorscale: Plotly colorscale name (default: "Viridis")
        show_plot: Whether to display the plot (default: True)
        
    Returns:
        plotly.graph_objects.Figure
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> fig = plot_simmat_plotly(db)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for plotting. Install with: pip install plotly")
    
    if "sca" not in db:
        raise ValueError("Database must contain 'sca' key")
    
    sca_data = db["sca"]
    if "simMat" not in sca_data:
        raise ValueError("Database must contain 'sca.simMat' (run with --do-seqcorr)")
    
    simMat = sca_data["simMat"]
    
    # Convert sparse matrix to dense if needed
    from scipy.sparse import issparse
    if issparse(simMat):
        simMat = simMat.toarray()
    
    # Get headers if available
    headers = None
    if "sequence" in db and "hd" in db["sequence"]:
        headers = db["sequence"]["hd"]
        # If simMat was subsampled, we might have selected_indices
        if "selected_indices" in sca_data:
            selected_indices = sca_data["selected_indices"]
            headers = [headers[i] if i < len(headers) else f"Seq {i}" for i in selected_indices]
    
    # Create hover text
    M = simMat.shape[0]
    hovertext = []
    for i in range(M):
        row = []
        for j in range(M):
            seq_i = headers[i] if headers and i < len(headers) else f"Seq {i}"
            seq_j = headers[j] if headers and j < len(headers) else f"Seq {j}"
            row.append(f"Seq {i} ({seq_i})<br>Seq {j} ({seq_j})<br>Similarity: {simMat[i, j]:.4f}")
        hovertext.append(row)
    
    fig = go.Figure(data=go.Heatmap(
        z=simMat,
        colorscale=colorscale,
        colorbar=dict(title="Similarity"),
        hovertext=hovertext,
        hovertemplate='%{hovertext}<extra></extra>'
    ))
    
    fig.update_layout(
        title='Sequence Similarity Matrix',
        xaxis_title='Sequence Index',
        yaxis_title='Sequence Index',
        height=height,
        width=width,
        template='plotly_white'
    )
    
    if show_plot:
        # For Cursor/VS Code, use HTML display since mimetype isn't supported
        # This works in both Jupyter Lab (which supports mimetype) and Cursor (which doesn't)
        try:
            from IPython.display import HTML, display
            # Use HTML display which works in both Jupyter and Cursor/VS Code
            html_str = fig.to_html(include_plotlyjs='cdn', div_id=f"plotly-{id(fig)}")
            display(HTML(html_str))
        except ImportError:
            # If IPython not available, try regular show (will open in browser)
            fig.show()
        except Exception:
            # Final fallback: regular show
            fig.show()
    
    return fig


def plot_ic_correlation_plotly(
    db: Dict[str, Any],
    height: int = 600,
    width: int = 700,
    colorscale: str = "RdBu",
    show_plot: bool = True,
) -> go.Figure:
    """
    Plot Spearman correlation matrix between independent components using Plotly.
    
    Parameters:
        db: Database dictionary from load_database()
        height: Plot height in pixels
        width: Plot width in pixels
        colorscale: Plotly colorscale name (default: "RdBu")
        show_plot: Whether to display the plot (default: True)
        
    Returns:
        plotly.graph_objects.Figure
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> fig = plot_ic_correlation_plotly(db)
        >>> fig.show()
    """
    if not PLOTLY_AVAILABLE:
        raise ImportError("plotly is required for plotting. Install with: pip install plotly")
    
    if "sector" not in db:
        raise ValueError("Database must contain 'sector' key")
    
    sector_data = db["sector"]
    if "ic_corr_spearman" not in sector_data:
        raise ValueError("Database must contain 'sector.ic_corr_spearman'")
    
    ic_corr = sector_data["ic_corr_spearman"]
    n_ics = len(ic_corr)
    
    # Create hover text with correlation values
    hovertext = []
    for i in range(n_ics):
        row = []
        for j in range(n_ics):
            row.append(f"IC {i+1} vs IC {j+1}<br>Correlation: {ic_corr[i, j]:.4f}")
        hovertext.append(row)
    
    fig = go.Figure(data=go.Heatmap(
        z=ic_corr,
        x=[f"IC {i+1}" for i in range(n_ics)],
        y=[f"IC {i+1}" for i in range(n_ics)],
        colorscale=colorscale,
        zmid=0,  # Center colorscale at 0
        zmin=-1,
        zmax=1,
        colorbar=dict(title="Spearman ρ"),
        hovertext=hovertext,
        hovertemplate='%{hovertext}<extra></extra>',
        text=[[f"{ic_corr[i, j]:.3f}" for j in range(n_ics)] for i in range(n_ics)],
        texttemplate="%{text}",
        textfont={"size": 10}
    ))
    
    fig.update_layout(
        title='Spearman Correlation Between Independent Components',
        xaxis_title='Independent Component',
        yaxis_title='Independent Component',
        height=height,
        width=width,
        template='plotly_white'
    )
    
    if show_plot:
        # For Cursor/VS Code, use HTML display since mimetype isn't supported
        # This works in both Jupyter Lab (which supports mimetype) and Cursor (which doesn't)
        try:
            from IPython.display import HTML, display
            # Use HTML display which works in both Jupyter and Cursor/VS Code
            html_str = fig.to_html(include_plotlyjs='cdn', div_id=f"plotly-{id(fig)}")
            display(HTML(html_str))
        except ImportError:
            # If IPython not available, try regular show (will open in browser)
            fig.show()
        except Exception:
            # Final fallback: regular show
            fig.show()
    
    return fig


# ============================================================================
# Data Summary Functions
# ============================================================================

def summarize_database(db: Dict[str, Any]) -> Dict[str, Any]:
    """
    Print a summary of database contents.
    
    Parameters:
        db: Database dictionary from load_database()
        
    Returns:
        Dictionary with summary statistics
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> summary = summarize_database(db)
    """
    summary = {}
    
    if "sequence" in db:
        seq = db["sequence"]
        # Database uses Nseq, effseqs, Npos (not M, M_eff, L)
        # Direct access since we know these keys exist in current databases
        summary["sequences"] = {
            "M": seq.get("Nseq", seq.get("M", "N/A")),
            "M_eff": seq.get("effseqs", seq.get("M_eff", "N/A")),
            "L": seq.get("Npos", seq.get("L", "N/A")),
            "i_ref": seq.get("i_ref", "N/A"),
        }
    
    if "sca" in db:
        sca_data = db["sca"]
        summary["sca"] = {
            "Csca_shape": sca_data.get("Csca", np.array([])).shape if "Csca" in sca_data else "N/A",
            "has_eigenvalues": "L" in sca_data,
            "has_simMat": "simMat" in sca_data,
            "Ntrials": sca_data.get("Ntrials", "N/A"),
        }
    
    if "sector" in db:
        sector_data = db["sector"]
        summary["sector"] = {
            "kpos": sector_data.get("kpos", "N/A"),
            "kpos_auto": sector_data.get("kpos_auto", "N/A"),
            "n_ics": len(sector_data.get("ic_list", [])) if "ic_list" in sector_data else "N/A",
            "has_Vica": "Vica" in sector_data,
        }
    
    return summary


def print_summary(db: Dict[str, Any]) -> None:
    """
    Print a formatted summary of database contents.
    
    Parameters:
        db: Database dictionary from load_database()
        
    Example:
        >>> db = load_database("Outputs/alignment.db.gz")
        >>> print_summary(db)
    """
    summary = summarize_database(db)
    
    print("=" * 60)
    print("Database Summary")
    print("=" * 60)
    
    if "sequences" in summary:
        print("\nSequence Data:")
        for key, value in summary["sequences"].items():
            print(f"  {key:12s}: {value}")
    
    if "sca" in summary:
        print("\nSCA Data:")
        for key, value in summary["sca"].items():
            print(f"  {key:12s}: {value}")
    
    if "sector" in summary:
        print("\nSector/IC Data:")
        for key, value in summary["sector"].items():
            print(f"  {key:12s}: {value}")
    
    print("=" * 60)

