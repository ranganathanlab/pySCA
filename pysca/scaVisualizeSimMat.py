#!/usr/bin/env python3
"""
scaVisualizeSimMat.py

Interactive visualization of the sequence similarity matrix (simMat) using Bokeh.

This script loads a database file from sca-core and creates an interactive
heatmap visualization of the simMat matrix. The visualization supports:
- Zoom and pan
- Hover tooltips showing sequence indices and similarity values
- Color mapping with customizable colormap
- Export to HTML or PNG

Input: Database file (.db or .db.gz) containing simMat in db['sca']['simMat']
Output: Interactive HTML visualization

Requires:
  - bokeh (install with: pip install bokeh)
  - numpy, scipy (for handling sparse matrices if present)
"""

from __future__ import annotations

import argparse
import gzip
import pickle
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from scipy.sparse import issparse

try:
    from bokeh.io import output_file, save, show, export_png
    from bokeh.plotting import figure
    from bokeh.models import ColorBar, LinearColorMapper, HoverTool, ColumnDataSource
    from bokeh.layouts import column
    from bokeh.models.widgets import Div
    from bokeh.palettes import Viridis256
    BOKEH_AVAILABLE = True
except ImportError:
    BOKEH_AVAILABLE = False


def load_db(path: Path) -> Dict[str, Any]:
    """Load database file (.db or .db.gz)."""
    if path.suffix == ".gz" or path.name.endswith(".db.gz"):
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    with open(path, "rb") as f:
        return pickle.load(f)


def extract_simmat(db: Dict[str, Any]) -> tuple[np.ndarray, Optional[list], Optional[list]]:
    """
    Extract simMat from database.
    
    Returns:
        simMat: numpy array (M x M)
        headers: sequence headers if available (None if not)
        selected_indices: indices of sequences if subsampled (None if not)
    """
    if "sca" not in db:
        raise KeyError("Database missing 'sca' key. Did you run sca-core with --do-seqcorr?")
    
    sca_data = db["sca"]
    if "simMat" not in sca_data:
        raise KeyError("Database missing 'sca.simMat'. Did you run sca-core with --do-seqcorr?")
    
    simMat = sca_data["simMat"]
    
    # Convert sparse matrix to dense if needed
    if issparse(simMat):
        print(f"Converting sparse matrix to dense array (shape: {simMat.shape})...")
        simMat = simMat.toarray()
    
    # Ensure it's a numpy array
    simMat = np.asarray(simMat)
    
    # Ensure it's 2D
    if simMat.ndim != 2:
        raise ValueError(f"simMat must be 2D, got shape {simMat.shape}")
    
    # Get headers if available
    headers = None
    selected_indices = None
    
    if "sequence" in db:
        seq_data = db["sequence"]
        if "hd" in seq_data:
            all_headers = seq_data["hd"]
            
            # Check if we have selected_indices from subsampling
            # This would be stored if seqSim was called with return_indices=True
            # For now, we'll use all headers if available
            if len(all_headers) == simMat.shape[0]:
                headers = all_headers
            # If simMat was subsampled, we might not have the exact headers
            # In that case, headers will be None
    
    return simMat, headers, selected_indices


def create_heatmap(
    simMat: np.ndarray,
    headers: Optional[list] = None,
    output_path: Optional[Path] = None,
    title: str = "Sequence Similarity Matrix",
    colormap: str = "Viridis256",
    width: int = 800,
    height: int = 800,
    show_plot: bool = True,
) -> None:
    """
    Create interactive Bokeh heatmap of simMat.
    
    Parameters:
        simMat: 2D numpy array (M x M) of similarity values
        headers: Optional list of sequence headers/names
        output_path: Path to save HTML file (optional)
        title: Plot title
        colormap: Bokeh palette name (default: Viridis256)
        width: Plot width in pixels
        height: Plot height in pixels
        show_plot: Whether to open plot in browser
    """
    if not BOKEH_AVAILABLE:
        raise ImportError(
            "Bokeh is required for visualization. Install with: pip install bokeh"
        )
    
    M = simMat.shape[0]
    
    # Get colormap
    if colormap == "Viridis256":
        palette = Viridis256
    else:
        # Try to get from bokeh.palettes
        try:
            from bokeh.palettes import all_palettes
            if colormap in all_palettes:
                # Use the largest available palette
                palette = all_palettes[colormap][max(all_palettes[colormap].keys())]
            else:
                print(f"Warning: Colormap '{colormap}' not found, using Viridis256")
                palette = Viridis256
        except Exception:
            palette = Viridis256
    
    # Find min/max for color mapping
    vmin = float(np.nanmin(simMat))
    vmax = float(np.nanmax(simMat))
    
    # Create color mapper
    color_mapper = LinearColorMapper(
        palette=palette,
        low=vmin,
        high=vmax,
        nan_color="gray",
    )
    
    # Prepare data for Bokeh
    # Bokeh's image glyph works with row-major data
    # We'll use image_rgba for better performance with large matrices
    # For smaller matrices, we can use rect glyphs for better interactivity
    
    if M <= 1000:
        # Use rect glyphs for smaller matrices (better interactivity)
        # Create arrays for rectangle centers
        x = np.arange(M)
        y = np.arange(M)
        xx, yy = np.meshgrid(x, y)
        
        # Flatten for plotting
        x_flat = xx.flatten()
        y_flat = yy.flatten()
        values_flat = simMat.flatten()
        
        # Prepare data source
        if headers:
            # Add header columns to source
            x_headers = [headers[int(xi)] if int(xi) < len(headers) else f"seq{int(xi)}" for xi in x_flat]
            y_headers = [headers[int(yi)] if int(yi) < len(headers) else f"seq{int(yi)}" for yi in y_flat]
            source = ColumnDataSource(data={
                "x": x_flat,
                "y": y_flat,
                "value": values_flat,
                "x_header": x_headers,
                "y_header": y_headers,
            })
        else:
            source = ColumnDataSource(data={
                "x": x_flat,
                "y": y_flat,
                "value": values_flat,
            })
        
        # Create hover tool
        if headers:
            hover = HoverTool(tooltips=[
                ("Seq i", "@y{0}"),
                ("Seq j", "@x{0}"),
                ("Similarity", "@value{0.000}"),
                ("Header i", "@y_header"),
                ("Header j", "@x_header"),
            ])
        else:
            hover = HoverTool(tooltips=[
                ("Seq i", "@y{0}"),
                ("Seq j", "@x{0}"),
                ("Similarity", "@value{0.000}"),
            ])
        
        # Create figure
        p = figure(
            title=title,
            x_range=(0, M),
            y_range=(0, M),
            width=width,
            height=height,
            tools=[hover, "pan", "box_zoom", "wheel_zoom", "reset", "save"],
            toolbar_location="above",
        )
        
        # Add rectangles
        p.rect(
            x="x",
            y="y",
            width=1,
            height=1,
            fill_color={"field": "value", "transform": color_mapper},
            line_color=None,
            source=source,
        )
    else:
        # For larger matrices, use image glyph (more efficient)
        # Convert to RGBA image
        # Normalize to 0-255 range
        norm_data = (simMat - vmin) / (vmax - vmin)
        norm_data = np.clip(norm_data, 0, 1)
        
        # Map to palette indices (0-255)
        indices = (norm_data * (len(palette) - 1)).astype(int)
        
        # Create RGBA array
        image_rgba = np.zeros((M, M, 4), dtype=np.uint8)
        for i in range(M):
            for j in range(M):
                if np.isnan(simMat[i, j]):
                    image_rgba[i, j] = [128, 128, 128, 255]  # Gray for NaN
                else:
                    color_hex = palette[indices[i, j]]
                    # Convert hex to RGB
                    rgb = tuple(int(color_hex[k:k+2], 16) for k in (1, 3, 5))
                    image_rgba[i, j] = [rgb[0], rgb[1], rgb[2], 255]
        
        # Create figure
        p = figure(
            title=title,
            x_range=(0, M),
            y_range=(0, M),
            width=width,
            height=height,
            tools="pan,box_zoom,wheel_zoom,reset,save",
            toolbar_location="above",
        )
        
        # Add image
        p.image_rgba(
            image=[image_rgba],
            x=0,
            y=0,
            dw=M,
            dh=M,
        )
        
        # For large matrices, add a simpler hover (shows coordinates only)
        p.add_tools(HoverTool(
            tooltips=[
                ("Seq i", "$y{0}"),
                ("Seq j", "$x{0}"),
            ],
            mode="mouse",
        ))
    
    # Add colorbar
    color_bar = ColorBar(
        color_mapper=color_mapper,
        label_standoff=12,
        border_line_color=None,
        location=(0, 0),
    )
    p.add_layout(color_bar, "right")
    
    # Add axis labels
    p.xaxis.axis_label = "Sequence Index j"
    p.yaxis.axis_label = "Sequence Index i"
    
    # Create info div
    info_text = f"""
    <div style="font-family: Arial; font-size: 12px; margin: 10px;">
        <b>Matrix Info:</b><br>
        Size: {M} × {M}<br>
        Value range: [{vmin:.4f}, {vmax:.4f}]<br>
        Mean: {np.nanmean(simMat):.4f}<br>
        Std: {np.nanstd(simMat):.4f}
    </div>
    """
    info_div = Div(text=info_text, width=width, height=100)
    
    # Combine plot and info
    layout = column(info_div, p)
    
    # Save or show
    if output_path:
        output_file(str(output_path))
        save(layout)
        print(f"Saved visualization to: {output_path}")
    
    if show_plot:
        show(layout)


def main(argv: Optional[list[str]] = None) -> int:
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Interactive visualization of sequence similarity matrix (simMat) using Bokeh",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (opens in browser)
  sca-visualize-simmat Outputs/alignment.db.gz

  # Save to HTML file
  sca-visualize-simmat Outputs/alignment.db.gz -o simmat_plot.html

  # Custom size and colormap
  sca-visualize-simmat Outputs/alignment.db.gz -o plot.html --width 1200 --height 1200 --colormap Inferno256

  # Save as PNG (requires selenium)
  sca-visualize-simmat Outputs/alignment.db.gz -o plot.png
        """
    )
    
    parser.add_argument(
        "database",
        type=Path,
        help="Database file from sca-core (.db or .db.gz) containing simMat",
    )
    
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output file path (HTML or PNG). If not specified, opens in browser only.",
    )
    
    parser.add_argument(
        "--width",
        type=int,
        default=800,
        help="Plot width in pixels (default: 800)",
    )
    
    parser.add_argument(
        "--height",
        type=int,
        default=800,
        help="Plot height in pixels (default: 800)",
    )
    
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Plot title (default: auto-generated from database name)",
    )
    
    parser.add_argument(
        "--colormap",
        type=str,
        default="Viridis256",
        help="Colormap name from Bokeh palettes (default: Viridis256). "
             "Options: Viridis256, Inferno256, Magma256, Plasma256, etc.",
    )
    
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Don't open plot in browser (only save if --output is specified)",
    )
    
    args = parser.parse_args(argv)
    
    # Check Bokeh availability
    if not BOKEH_AVAILABLE:
        print("Error: Bokeh is required but not installed.")
        print("Install with: pip install bokeh")
        return 1
    
    # Load database
    db_path = args.database
    if not db_path.exists():
        print(f"Error: Database file not found: {db_path}")
        return 1
    
    print(f"Loading database: {db_path}")
    try:
        db = load_db(db_path)
    except Exception as e:
        print(f"Error loading database: {e}")
        return 1
    
    # Extract simMat
    print("Extracting simMat from database...")
    try:
        simMat, headers, selected_indices = extract_simmat(db)
        print(f"simMat shape: {simMat.shape}")
        print(f"Value range: [{np.nanmin(simMat):.4f}, {np.nanmax(simMat):.4f}]")
    except Exception as e:
        print(f"Error extracting simMat: {e}")
        return 1
    
    # Determine output path
    output_path = args.output
    if output_path is None and args.no_show:
        # Generate default output name
        output_path = db_path.with_suffix(".simmat.html")
    
    # Determine title
    title = args.title
    if title is None:
        title = f"Sequence Similarity Matrix: {db_path.stem}"
    
    # Create visualization
    try:
        create_heatmap(
            simMat=simMat,
            headers=headers,
            output_path=output_path,
            title=title,
            colormap=args.colormap,
            width=args.width,
            height=args.height,
            show_plot=not args.no_show,
        )
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

