# scaVisualizeSimMat - Interactive SimMat Visualization

Interactive visualization of the sequence similarity matrix (simMat) using Bokeh.

## Installation

Install Bokeh (required dependency):
```bash
pip install bokeh
```

## Usage

### Basic Usage

Visualize simMat from a database file (opens in browser):
```bash
sca-visualize-simmat Outputs/alignment.db.gz
```

Or run the Python script directly:
```bash
python pysca/scaVisualizeSimMat.py Outputs/alignment.db.gz
```

### Save to HTML

Save the visualization to an HTML file:
```bash
sca-visualize-simmat Outputs/alignment.db.gz -o simmat_plot.html
```

### Custom Options

Customize plot size, colormap, and title:
```bash
sca-visualize-simmat Outputs/alignment.db.gz \
    -o plot.html \
    --width 1200 \
    --height 1200 \
    --colormap Inferno256 \
    --title "My SimMat Visualization"
```

### Available Colormaps

Common colormap options (from Bokeh palettes):
- `Viridis256` (default) - perceptually uniform, colorblind-friendly
- `Inferno256` - dark theme, high contrast
- `Magma256` - purple to yellow
- `Plasma256` - purple to pink to yellow
- `Turbo256` - rainbow-like

## Requirements

- Database file from `sca-core` with `--do-seqcorr` flag
- The database must contain `db['sca']['simMat']`

## Features

- **Interactive heatmap** with zoom, pan, and hover tooltips
- **Hover tooltips** showing:
  - Sequence indices (i, j)
  - Similarity values
  - Sequence headers (if available)
- **Colorbar** showing value scale
- **Matrix statistics** displayed above the plot
- **Optimized rendering**:
  - Small matrices (≤1000 sequences): Rect glyphs for precise hover
  - Large matrices (>1000 sequences): Image glyph for performance

## Example Workflow

1. Run sca-core with sequence correlations:
```bash
sca-core Outputs/alignment.db.gz --do-seqcorr
```

2. Visualize the simMat:
```bash
sca-visualize-simmat Outputs/alignment.db.gz -o simmat.html
```

3. Open the HTML file in your browser for interactive exploration.

## Command-Line Options

```
positional arguments:
  database              Database file from sca-core (.db or .db.gz) containing simMat

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output file path (HTML or PNG). If not specified, opens in browser only.
  --width WIDTH         Plot width in pixels (default: 800)
  --height HEIGHT       Plot height in pixels (default: 800)
  --title TITLE         Plot title (default: auto-generated from database name)
  --colormap COLORMAP   Colormap name from Bokeh palettes (default: Viridis256)
  --no-show             Don't open plot in browser (only save if --output is specified)
```

## Notes

- For large matrices (>1000 sequences), the visualization uses image rendering for performance
- Hover tooltips work best with smaller matrices (≤1000 sequences)
- The script automatically handles sparse matrices by converting them to dense arrays
- NaN values are displayed in gray

## Troubleshooting

**Error: "Database missing 'sca.simMat'"**
- Run `sca-core` with the `--do-seqcorr` flag to compute simMat

**Error: "Bokeh is required but not installed"**
- Install Bokeh: `pip install bokeh`

**Large matrices are slow to render**
- This is expected for very large matrices (>10k sequences)
- Consider subsampling when running `sca-core` with `--max-seqcorr-seqs`


