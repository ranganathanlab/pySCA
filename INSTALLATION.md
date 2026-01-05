# pySCA Installation Guide

Complete installation and verification instructions for pySCA and all dependencies.

---

## Table of Contents

1. [Python Package Dependencies](#python-package-dependencies)
2. [External Command-Line Tools](#external-command-line-tools)
3. [Installing pySCA](#installing-pysca)
4. [Verification](#verification)
5. [Troubleshooting](#troubleshooting)

---

## Python Package Dependencies

pySCA requires the following Python packages, which are automatically installed when you install pySCA:

- **biopython**: Sequence alignment and PDB structure handling
- **numpy**: Numerical computations
- **scipy**: Scientific computing (sparse matrices, statistics)
- **matplotlib**: Plotting and visualization
- **argparse**: Command-line argument parsing (included in Python 3.2+)
- **wheel**: Package building (for development)

### Installing Python Dependencies

If installing pySCA from source, you can install dependencies manually:

```bash
pip install biopython numpy scipy matplotlib
```

Or use the requirements file (if available):

```bash
pip install -r requirements.txt
```

---

## External Command-Line Tools

pySCA requires two external command-line tools that must be installed separately:

### FASTA36 (for ggsearch36)

FASTA36 is **required** for reference sequence searching in `sca-process-msa`. The `ggsearch36` program from the FASTA36 package is used to find the best matching sequence in the alignment when using PDB structures or reference sequence files.

#### Linux

```bash
# Clone and compile FASTA36
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make -j2 -f ../make/Makefile.linux all
sudo cp -r ../bin /usr/local
sudo rm /usr/local/bin/README
cd ../..

# Verify installation
which ggsearch36
ggsearch36 -h
```

#### macOS

```bash
# Install via Homebrew (recommended)
brew install fasta

# Or compile from source
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make -j2 -f ../make/Makefile.osx all
sudo cp -r ../bin /usr/local
sudo rm /usr/local/bin/README
cd ../..

# Verify installation
which ggsearch36
ggsearch36 -h
```

#### Windows

1. Install MSYS2 (see [FASTA36 documentation](https://github.com/wrpearson/fasta36))
2. In MSYS2 terminal:
```bash
pacman -Syu
pacman -S git make gcc
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make CC=/usr/bin/gcc LD=/usr/bin/ld -j2 -f ../make/Makefile.linux all
cp -r ../bin /usr/local/
rm /usr/local/bin/README
```
3. Add `/usr/local/bin` to your Windows PATH (see [Windows PATH setup](https://github.com/wrpearson/fasta36#windows))

#### Verify Installation

After installation, verify that `ggsearch36` is in your PATH:

```bash
which ggsearch36
ggsearch36 -h
```

If the command is not found, ensure the installation directory is in your system PATH.

### MMseqs2

MMseqs2 is **optional** but **highly recommended** for preclustering large MSAs (>50k sequences). It significantly speeds up processing and reduces memory usage.

#### Linux/macOS (via Conda - Recommended)

```bash
conda install -c bioconda mmseqs2
```

#### Linux (via Package Manager)

```bash
# Ubuntu/Debian
sudo apt-get install mmseqs2

# Or download from GitHub releases
wget https://github.com/soedinglab/MMseqs2/releases/latest/download/mmseqs-linux-avx2.tar.gz
tar -xzf mmseqs-linux-avx2.tar.gz
sudo cp mmseqs/bin/mmseqs /usr/local/bin/
```

#### macOS (via Homebrew)

```bash
brew install mmseqs2
```

#### Manual Installation

Download pre-compiled binaries from [MMseqs2 releases](https://github.com/soedinglab/MMseqs2/releases) or compile from source:

```bash
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make -j4
make install
export PATH=$PATH:$(pwd)/bin
```

#### Verify Installation

```bash
which mmseqs
mmseqs version
```

If the command is not found, ensure the installation directory is in your system PATH.

---

## Installing pySCA

### Option 1: Development Installation (Recommended)

Install pySCA in editable mode for development:

```bash
cd pySCA
pip install -e .
```

This installs pySCA and makes the command-line scripts (`sca-process-msa`, `sca-core`) available in your PATH.

### Option 2: Standard Installation

```bash
cd pySCA
pip install .
```

### Option 3: Install from Source (No pip)

If you don't want to use pip, you can add the pySCA directory to your Python path:

```bash
export PYTHONPATH=$PYTHONPATH:/path/to/pySCA
```

Then run scripts directly:

```bash
python pysca/scaProcessMSA_py3_big.py --help
python pysca/scaCore_py3.py --help
```

---

## Verification

### Verify Python Package Installation

```bash
python -c "import pysca; print('pySCA installed successfully')"
python -c "import pysca.scaTools as sca; print('scaTools imported successfully')"
```

### Verify Command-Line Scripts

After installation with `pip install -e .`, verify the scripts are available:

```bash
which sca-process-msa
which sca-core
sca-process-msa --help
sca-core --help
```

### Verify External Tools

```bash
# Verify ggsearch36 (required)
which ggsearch36
ggsearch36 -h

# Verify mmseqs2 (optional but recommended)
which mmseqs
mmseqs version
```

### Quick Test Run

Test that pySCA can import and run basic functions:

```python
import numpy as np
from pysca import scaTools as sca

# Test basic function
alg = ["ACDEFGHIKLMNPQRSTVWY", "ACDEFGHIKLMNPQRSTVWY-"]
msa_num = sca.lett2num(alg)
print(f"Test alignment shape: {msa_num.shape}")
print("✓ Basic functions working")
```

---

## Troubleshooting

### "ggsearch36 not found" or "ggsearch36 failed" errors

**Symptoms:**
- `RuntimeError: ggsearch36 failed: ... Make sure ggsearch36 is installed and in your PATH.`
- `FileNotFoundError: ggsearch36 command not found`
- Error occurs when using `-s` (PDB reference) or `--refseq` options

**Solutions:**

1. **Install FASTA36** (see [FASTA36 Installation](#fasta36-for-ggsearch36) above)

2. **Verify installation:**
   ```bash
   which ggsearch36
   ggsearch36 -h
   ```

3. **If installed but not found:**
   - Check that the installation directory (usually `/usr/local/bin`) is in your PATH
   - Add to PATH if needed:
     ```bash
     export PATH=$PATH:/usr/local/bin
     ```
   - For permanent fix, add to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.)

4. **Alternative:** If you cannot install FASTA36, you can use `--refindex` to specify the reference sequence by index instead:
   ```bash
   sca-process-msa alignment.fasta --refindex 0  # Use first sequence
   ```

### "MMseqs2 not found" errors

**Symptoms:**
- Error when using `--precluster` option
- `FileNotFoundError: mmseqs command not found`

**Solutions:**
- Install MMseqs2: `conda install -c bioconda mmseqs2` or download from https://github.com/soedinglab/MMseqs2
- Ensure MMseqs2 is in PATH: `which mmseqs`
- If MMseqs2 is not available, you can disable preclustering with `--no-precluster` (not recommended for large MSAs)

### Python Import Errors

**Symptoms:**
- `ModuleNotFoundError: No module named 'pysca'`
- `ModuleNotFoundError: No module named 'numpy'`

**Solutions:**

1. **Ensure pySCA is installed:**
   ```bash
   pip install -e .
   ```

2. **Check Python environment:**
   ```bash
   python --version  # Should be Python 3.6+
   which python
   ```

3. **Install missing dependencies:**
   ```bash
   pip install biopython numpy scipy matplotlib
   ```

4. **If using conda:**
   ```bash
   conda install biopython numpy scipy matplotlib
   ```

### Command-Line Scripts Not Found

**Symptoms:**
- `command not found: sca-process-msa`
- `command not found: sca-core`

**Solutions:**

1. **Reinstall in editable mode:**
   ```bash
   pip install -e .
   ```

2. **Check installation location:**
   ```bash
   pip show pySCA
   ```

3. **Verify scripts directory:**
   ```bash
   ls bin/sca-process-msa
   ls bin/sca-core
   ```

4. **Run directly from source:**
   ```bash
   python pysca/scaProcessMSA_py3_big.py --help
   python pysca/scaCore_py3.py --help
   ```

### PATH Issues

If tools are installed but not found:

1. **Check current PATH:**
   ```bash
   echo $PATH
   ```

2. **Add installation directories:**
   ```bash
   export PATH=$PATH:/usr/local/bin
   ```

3. **Make permanent (add to `~/.bashrc` or `~/.zshrc`):**
   ```bash
   echo 'export PATH=$PATH:/usr/local/bin' >> ~/.bashrc
   source ~/.bashrc
   ```

---

## System Requirements

### Minimum Requirements

- **Python:** 3.6 or higher
- **RAM:** 4 GB (8 GB recommended for large MSAs)
- **Disk Space:** 500 MB for installation, additional space for data and outputs

### Recommended for Large MSAs (>100k sequences)

- **Python:** 3.8 or higher
- **RAM:** 16 GB or more
- **Disk Space:** 10+ GB for temporary files and outputs
- **MMseqs2:** Required for efficient processing

---

## Next Steps

After successful installation and verification:

1. Read the [Usage Instructions](USAGE_INSTRUCTIONS.md) for `sca-process-msa` and `sca-core`
2. Review [scaTools Functions](SCATOOLS_FUNCTIONS.md) for available functions
3. Check [Optimizations](OPTIMIZATIONS.md) for performance tips

---

For more information, see:
- [pySCA Website](https://ranganathanlab.gitlab.io/pySCA)
- [GitHub Repository](https://github.com/ranganathanlab/pySCA)


