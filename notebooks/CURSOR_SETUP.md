# Running Notebooks in Cursor AI

Cursor (based on VS Code) has built-in Jupyter notebook support. You can run notebooks directly in Cursor without needing a separate Jupyter server.

## Setup Steps

### 1. Install Python Extension (if not already installed)

Cursor should have Python support built-in, but make sure you have:
- Python extension (usually pre-installed)
- Jupyter extension (usually pre-installed)

### 2. Install Notebook Dependencies

In your terminal (within Cursor or external):

```bash
# Activate your conda environment
conda activate pysca3

# Install notebook dependencies
pip install -e .[notebooks]
```

### 3. Select Python Interpreter in Cursor

**Method 1: From the Notebook (Recommended)**

1. Open `notebooks/SCA_PDZ_Example.ipynb` in Cursor
2. Look at the top-right of the notebook - you should see a kernel selector (may show "Select Kernel" or a Python version)
3. Click on it
4. Choose "Select Another Kernel..."
5. Choose "Python Environments..."
6. If you see `pysca3`, select it. If not, continue to Method 2

**Method 2: Manual Path Selection**

If `pysca3` doesn't appear in the list:

1. In the kernel selector, choose "Select Another Kernel..."
2. Choose "Python Environments..."
3. Click "Enter interpreter path..." or "Browse..."
4. Navigate to or paste this path:
   ```
   /opt/miniconda3/envs/pysca3/bin/python
   ```
5. Select it

**Method 3: Using Command Palette**

1. Press `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
2. Type "Python: Select Interpreter"
3. Choose "Enter interpreter path..."
4. Paste: `/opt/miniconda3/envs/pysca3/bin/python`
5. Press Enter

### 4. Run Cells

- Click the "Run" button above each cell
- Or use `Shift+Enter` to run the current cell
- Or use `Cmd+Enter` / `Ctrl+Enter` to run without moving to next cell

## Advantages of Using Cursor

✅ **Integrated experience** - Edit code and run notebooks in one place  
✅ **No separate Jupyter server needed** - Runs directly in Cursor  
✅ **Better code editing** - Full IDE features (autocomplete, linting, etc.)  
✅ **Git integration** - Easy version control for notebooks  
✅ **AI assistance** - Cursor's AI can help with notebook code  

## Troubleshooting

### Kernel Not Found / Environment Not Showing

If `pysca3` doesn't appear in the kernel list:

**Solution 1: Manually Enter Path**
1. In the notebook, click the kernel selector (top-right)
2. Choose "Select Another Kernel..." → "Python Environments..."
3. Click "Enter interpreter path..." or "Browse..."
4. Enter: `/opt/miniconda3/envs/pysca3/bin/python`
5. Select it

**Solution 2: Install ipykernel in the Environment**
Sometimes Cursor needs ipykernel installed to detect the environment:

```bash
# In terminal (with pysca3 activated)
conda activate pysca3
pip install ipykernel
```

Then restart Cursor and try selecting the interpreter again.

**Solution 3: Register the Kernel**
```bash
conda activate pysca3
python -m ipykernel install --user --name pysca3 --display-name "Python (pysca3)"
```

After registration, restart Cursor. The kernel should appear in the list.

### Import Errors

If you get import errors:
- Make sure pySCA is installed: `pip install -e .`
- Verify in terminal: `python -c "import pysca; print('OK')"`
- Restart Cursor after installing packages

### Plotly Not Displaying

Plotly plots should display inline in Cursor. If they don't:
- Make sure `plotly` is installed: `pip install plotly`
- Try restarting Cursor
- Check that the cell output is visible (may need to expand output)

## Quick Start

1. Open `notebooks/SCA_PDZ_Example.ipynb` in Cursor
2. Select `pysca3` as the kernel (top-right of notebook)
3. Run cells with `Shift+Enter`

That's it! No need for a separate Jupyter server.

