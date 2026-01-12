# pySCA Notebooks

Example notebooks for pySCA will be available here shortly.

## Coming Soon

New example notebooks demonstrating pySCA 7.0 features will be added here.

## Running Notebooks

To run notebooks in Jupyter:

### 1. Install Notebook Dependencies

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

### 2. Register Your Environment as a Jupyter Kernel

```bash
# Make sure you're in your conda/virtual environment
python -m ipykernel install --user --name pysca3 --display-name "Python (pysca3)"
```

### 3. Start Jupyter

```bash
# From the pySCA root directory
jupyter notebook

# Or use JupyterLab
jupyter lab
```

### 4. Open and Run Notebooks

1. Navigate to the `notebooks/` directory in Jupyter
2. Open a notebook
3. Select the appropriate kernel (e.g., "Python (pysca3)")
4. Execute cells one by one using `Shift+Enter`

## Troubleshooting

### Kernel Not Found

If you don't see your kernel in the kernel list:
- Make sure you registered the kernel (step 2)
- Try restarting Jupyter
- Check that you're in the correct conda/virtual environment

### Import Errors

If you get import errors:
- Make sure pySCA is installed: `pip install -e .`
- Verify installation: `python -c "import pysca; print('OK')"`
