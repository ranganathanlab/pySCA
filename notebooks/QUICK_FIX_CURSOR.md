# Quick Fix: Select pysca3 Environment in Cursor

If you don't see `pysca3` in the kernel selector, follow these steps:

## Step-by-Step Fix

### 1. Open the Notebook
- Open `notebooks/SCA_PDZ_Example.ipynb` in Cursor

### 2. Click the Kernel Selector
- Look at the **top-right** of the notebook
- You should see something like "Select Kernel" or a Python version number
- **Click on it**

### 3. Enter the Python Path Manually
- Choose **"Select Another Kernel..."**
- Choose **"Python Environments..."**
- Click **"Enter interpreter path..."** (or "Browse...")
- Paste this exact path:
  ```
  /opt/miniconda3/envs/pysca3/bin/python
  ```
- Press Enter or click "Select"

### 4. Verify It Works
- The kernel selector should now show something like "Python 3.11.x ('pysca3': conda)"
- Try running the first code cell (Setup and Imports)
- If you get import errors, see troubleshooting below

## Alternative: Use Command Palette

1. Press `Cmd+Shift+P` (Mac) or `Ctrl+Shift+P` (Windows/Linux)
2. Type: `Python: Select Interpreter`
3. Choose: `Enter interpreter path...`
4. Paste: `/opt/miniconda3/envs/pysca3/bin/python`
5. Press Enter

## If It Still Doesn't Work

### Install ipykernel
```bash
conda activate pysca3
pip install ipykernel
```

Then restart Cursor and try again.

### Register the Kernel
```bash
conda activate pysca3
python -m ipykernel install --user --name pysca3 --display-name "Python (pysca3)"
```

Restart Cursor after registration.

## Verify Installation

After selecting the kernel, run this in the first code cell to verify:

```python
import sys
print(sys.executable)  # Should show /opt/miniconda3/envs/pysca3/bin/python
```

If it shows a different path, the kernel selection didn't work - try the steps again.


