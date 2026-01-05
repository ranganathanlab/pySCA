# Jupyter Notebook Support for pySCA 7.0

pySCA 7.0 includes comprehensive support for Jupyter notebooks, allowing you to perform SCA analysis interactively with rich visualizations and data exploration.

## Installation

To use pySCA in Jupyter notebooks, install the optional notebook dependencies:

```bash
pip install -e .[notebooks]
```

This installs:
- `jupyter` - Jupyter notebook environment
- `ipykernel` - Jupyter kernel (required for conda environments)
- `ipywidgets` - Interactive widgets
- `plotly` - Interactive visualizations (primary, recommended)
- `bokeh` - Alternative visualization library

Alternatively, install individually:

```bash
pip install jupyter ipykernel ipywidgets plotly
```

### Registering Your Environment as a Jupyter Kernel

If you're using a conda environment (like `pysca3`), you need to register it as a Jupyter kernel:

```bash
# Activate your conda environment
conda activate pysca3

# Install ipykernel (if not already installed)
pip install ipykernel

# Register the environment as a Jupyter kernel
python -m ipykernel install --user --name pysca3 --display-name "Python (pysca3)"
```

After registration, you can select "Python (pysca3)" as the kernel when creating or opening a notebook in Jupyter.

## Quick Start

### 1. Import Notebook Utilities

```python
from pysca import notebook_utils as nb
from pysca import scaTools as sca
import plotly.graph_objects as go
import plotly.express as px
```

### 2. Process an MSA

```python
db_path = nb.process_msa(
    "Inputs/alignment.fasta",
    pdb="Inputs/structure.pdb",
    chain="A",
    species="Homo sapiens"
)
```

### 3. Run SCA Core Calculations

```python
output_db = nb.run_sca_core(
    db_path,
    do_sector_id=True,
    kpos=6
)
```

### 4. Load and Visualize Results

```python
db = nb.load_database(output_db)
nb.plot_eigenvalues(db)
nb.plot_ic_heatmap(db)
```

## Notebook Utilities API

### Data Loading

#### `load_database(db_path)`
Load a pySCA database file (.db or .db.gz).

**Parameters:**
- `db_path` (str | Path): Path to database file

**Returns:**
- Dictionary with 'sequence', 'sca', and optionally 'sector' keys

**Example:**
```python
db = nb.load_database("Outputs/alignment.db.gz")
```

### MSA Processing

#### `process_msa(...)`
Process a multiple sequence alignment. Wrapper around `sca-process-msa`.

**Key Parameters:**
- `alignment`: Path to input alignment
- `pdb`: PDB identifier or path (optional)
- `chain`: Chain ID (default: "A")
- `species`: Species name for reference search (optional)
- `precluster`: Enable/disable MMseqs2 preclustering (None = auto)
- `parameters`: List [pos_gap, seq_gap, min_seqid, max_seqid]

**Example:**
```python
db_path = nb.process_msa(
    "Inputs/alignment.fasta",
    pdb="Inputs/1ABC.pdb",
    chain="A",
    precluster=True
)
```

### SCA Core Calculations

#### `run_sca_core(...)`
Run SCA core calculations. Wrapper around `sca-core`.

**Key Parameters:**
- `database`: Path to database from `process_msa`
- `do_sector_id`: Perform independent component identification
- `kpos`: Number of eigenmodes (None = auto)
- `do_seqcorr`: Compute sequence correlations (memory-intensive)
- `Ntrials`: Number of randomization trials (default: 10)

**Example:**
```python
output_db = nb.run_sca_core(
    "Outputs/alignment.db.gz",
    do_sector_id=True,
    kpos=6
)
```

### Visualization Functions

pySCA provides both matplotlib (static) and Plotly (interactive) visualization functions. Plotly is recommended for notebooks as it provides rich interactivity.

#### Plotly Functions (Recommended)

##### `plot_eigenvalues_plotly(db, n_modes=None, show_randomized=True)`
Plot eigenvalues of the SCA matrix with optional randomized comparison using Plotly.

**Example:**
```python
fig = nb.plot_eigenvalues_plotly(db, n_modes=20)
fig.show()
```

##### `plot_ic_heatmap_plotly(db, ic_index=None)`
Plot heatmap of independent component values across positions using Plotly.

**Example:**
```python
# Plot all ICs
fig = nb.plot_ic_heatmap_plotly(db)

# Plot specific IC
fig = nb.plot_ic_heatmap_plotly(db, ic_index=0)
fig.show()
```

##### `plot_ic_positions_plotly(db, ic_index, top_n=None)`
Plot significant positions for a specific independent component using Plotly.

**Example:**
```python
fig = nb.plot_ic_positions_plotly(db, ic_index=0, top_n=20)
fig.show()
```

##### `plot_simmat_plotly(db)`
Create interactive Plotly heatmap of sequence similarity matrix.

**Requirements:** `plotly` installed and `--do-seqcorr` flag used in `run_sca_core`

**Example:**
```python
fig = nb.plot_simmat_plotly(db)
fig.show()
```

##### `plot_ic_correlation_plotly(db)`
Plot Spearman correlation matrix between independent components using Plotly.

**Example:**
```python
fig = nb.plot_ic_correlation_plotly(db)
fig.show()
```

#### Matplotlib Functions (Static)

##### `plot_eigenvalues(db, n_modes=None, show_randomized=True)`
Plot eigenvalues of the SCA matrix with optional randomized comparison.

**Example:**
```python
nb.plot_eigenvalues(db, n_modes=20)
```

##### `plot_ic_heatmap(db, ic_index=None)`
Plot heatmap of independent component values across positions.

**Example:**
```python
# Plot all ICs
nb.plot_ic_heatmap(db)

# Plot specific IC
nb.plot_ic_heatmap(db, ic_index=0)
```

##### `plot_ic_positions(db, ic_index, top_n=None)`
Plot significant positions for a specific independent component.

**Example:**
```python
nb.plot_ic_positions(db, ic_index=0, top_n=20)
```

#### Bokeh Functions (Alternative Interactive)

##### `plot_simmat_interactive(db)`
Create interactive Bokeh heatmap of sequence similarity matrix.

**Requirements:** `bokeh` installed and `--do-seqcorr` flag used in `run_sca_core`

**Example:**
```python
from bokeh.io import output_notebook
output_notebook()
nb.plot_simmat_interactive(db)
```

### Data Summary

#### `summarize_database(db)`
Get summary statistics from database.

**Returns:** Dictionary with summary information

#### `print_summary(db)`
Print formatted summary of database contents.

**Example:**
```python
db = nb.load_database("Outputs/alignment.db.gz")
nb.print_summary(db)
```

## Example Notebooks

Example notebooks are located in the `notebooks/` directory:

- `SCA_Example_Template.ipynb` - Complete workflow template
- `SCA_DHFR.ipynb` - DHFR family analysis (legacy, to be updated)
- `SCA_betalactamase.ipynb` - Beta-lactamase analysis (legacy, to be updated)
- `SCA_G.ipynb` - G-protein analysis (legacy, to be updated)
- `SCA_S1A.ipynb` - S1A protease analysis (legacy, to be updated)

## Workflow Example

Here's a complete example workflow:

```python
# 1. Setup
from pysca import notebook_utils as nb
import matplotlib.pyplot as plt
%matplotlib inline

# 2. Process MSA
db_path = nb.process_msa(
    "Inputs/alignment.fasta",
    pdb="Inputs/structure.pdb",
    chain="A",
    species="Homo sapiens",
    precluster=None  # Auto-enable for >50k sequences
)

# 3. Inspect processed data
db = nb.load_database(db_path)
nb.print_summary(db)

# 4. Run SCA core with IC identification
output_db = nb.run_sca_core(
    db_path,
    do_sector_id=True,
    kpos=None  # Auto-select
)

# 5. Visualize results
db = nb.load_database(output_db)

# Eigenvalue spectrum (Plotly - recommended)
fig = nb.plot_eigenvalues_plotly(db, n_modes=20)
fig.show()

# IC heatmap (Plotly)
fig = nb.plot_ic_heatmap_plotly(db)
fig.show()

# Individual IC positions (Plotly)
fig = nb.plot_ic_positions_plotly(db, ic_index=0, top_n=20)
fig.show()

# IC correlation matrix (Plotly)
fig = nb.plot_ic_correlation_plotly(db)
fig.show()

# 6. Access data programmatically
sector_data = db["sector"]
print(f"Number of ICs: {len(sector_data['ic_list'])}")
print(f"Significant positions for IC 1: {sector_data['sector_ats'][0]}")
```

## Tips and Best Practices

### Memory Management
- For large MSAs (>50k sequences), preclustering is automatically enabled
- Use `float32=True` in `run_sca_core()` for memory efficiency
- Sequence correlations (`--do-seqcorr`) require O(M²) memory - use subsampling

### IC Selection
- Let `kpos=None` to auto-select based on eigenvalue spectrum
- Inspect eigenvalue plot to manually choose `kpos`
- Use `kmax_cap` to limit maximum number of modes

### Visualization
- **Plotly (recommended)**: Use `plot_eigenvalues_plotly()` to assess how many modes to keep
- **Plotly**: `plot_ic_heatmap_plotly()` shows all ICs at once for comparison
- **Plotly**: `plot_ic_positions_plotly()` highlights significant positions for individual ICs
- **Plotly**: `plot_simmat_plotly()` for interactive sequence similarity matrices
- **Plotly**: `plot_ic_correlation_plotly()` for IC correlation analysis
- **Matplotlib**: Static versions available for publication-quality figures

### Data Access
- All data is stored in the database dictionary
- Access via `db["sequence"]`, `db["sca"]`, `db["sector"]`
- Use `print_summary()` to see what's available

## Integration with Command-Line Tools

The notebook utilities are wrappers around the command-line tools (`sca-process-msa` and `sca-core`). You can:

1. **Use notebooks for exploration** - Interactive analysis and visualization
2. **Use command-line for batch processing** - Automated workflows and scripts
3. **Mix both** - Run command-line tools, then load results in notebooks

The database format is identical, so you can seamlessly switch between interfaces.

## Troubleshooting

### Import Errors
If you get import errors, ensure pySCA is installed:
```bash
pip install -e .
```

### Missing Dependencies
Install notebook dependencies:
```bash
pip install jupyter ipywidgets bokeh
```

### Bokeh Not Showing
For Bokeh plots, you need to call:
```python
from bokeh.io import output_notebook
output_notebook()
```

### Memory Issues
- Use `precluster=True` for large MSAs
- Set `float32=True` in `run_sca_core()`
- Avoid `--do-seqcorr` for very large alignments

## Advanced Usage

### Custom Visualizations

You can create custom visualizations using the data directly:

```python
db = nb.load_database("Outputs/alignment.db.gz")

# Access raw data
Vica = db["sector"]["Vica"]  # IC vectors
L = db["sca"]["L"]  # Eigenvalues
Csca = db["sca"]["Csca"]  # SCA matrix

# Create custom plots
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
# ... your custom plotting code ...
```

### Programmatic Workflow

You can also call the underlying functions directly:

```python
from pysca import scaTools as sca

# Direct function calls
alg, headers = sca.readAlg("Inputs/alignment.fasta")
seqw = sca.seqWeights(alg)
# ... etc.
```

## See Also

- [Usage Instructions](USAGE_INSTRUCTIONS.md) - Command-line usage
- [Installation Guide](INSTALLATION.md) - Installation and dependencies
- [SCA Tools Functions](SCATOOLS_FUNCTIONS.md) - Core function reference

