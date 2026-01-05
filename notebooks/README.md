# Running pySCA Notebooks

You can run notebooks in two ways:
1. **In Cursor AI** (recommended) - See [CURSOR_SETUP.md](CURSOR_SETUP.md)
2. **In Jupyter** - Follow the steps below

## Quick Setup Guide for Jupyter

To run the example notebooks in Jupyter, follow these steps:

### 1. Activate Your Conda Environment

```bash
conda activate pysca3
```

### 2. Install Notebook Dependencies

```bash
# From the pySCA root directory
pip install -e .[notebooks]
```

This installs:
- `jupyter` - Jupyter notebook environment
- `ipykernel` - Jupyter kernel
- `ipywidgets` - Interactive widgets
- `plotly` - Interactive visualizations
- `bokeh` - Alternative visualization library

### 3. Register Your Environment as a Jupyter Kernel

```bash
# Make sure you're in the pysca3 environment
python -m ipykernel install --user --name pysca3 --display-name "Python (pysca3)"
```

### 4. Start Jupyter

```bash
# From the pySCA root directory (or notebooks directory)
jupyter notebook

# Or use JupyterLab
jupyter lab
```

### 5. Open and Run the Notebook

1. Navigate to the `notebooks/` directory in Jupyter
2. Open `SCA_PDZ_Example.ipynb`
3. **Important**: Select "Python (pysca3)" as the kernel from the kernel menu (Kernel → Change Kernel → Python (pysca3))
4. Execute cells one by one using `Shift+Enter`

## Troubleshooting

### Kernel Not Found

If you don't see "Python (pysca3)" in the kernel list:
- Make sure you registered the kernel (step 3)
- Try restarting Jupyter
- Check that you're in the correct conda environment: `conda activate pysca3`

### Import Errors

If you get import errors:
- Make sure pySCA is installed: `pip install -e .`
- Verify installation: `python -c "import pysca; print('OK')"`

### Database Not Found

If the notebook can't find `Outputs/PDZ_PF00595_aln.db.gz`:
- Make sure you've run `sca-process-msa` first to create the database
- Check that the file exists: `ls Outputs/PDZ_PF00595_aln.db.gz`

## Available Notebooks

- **SCA_PDZ_Example.ipynb** - Complete workflow example using PDZ domain data
- **SCA_Example_Template.ipynb** - Template for creating your own analyses

