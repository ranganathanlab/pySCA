# pySCA Usage Instructions

SCA represents a tool to examine the evolutionary pressures on and between amino acids in proteins. The analysis involves two steps: (1) `sca-process-msa`, which pre-processes a multiple sequence alignment to trim partial sequences, positions with excessive gaps, and to compute sequence and position weights, (2) `sca-core`, which carries out the core calculations of site-specific conservation, coevolution, and its decomposition to identify statistically significant collective groups of amino acids (sectors). Here, we provide a guide for using `sca-process-msa` and `sca-core`.

**Before you begin:** Make sure you have installed pySCA and all dependencies. See [INSTALLATION.md](INSTALLATION.md) for installation instructions.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [sca-process-msa - Alignment Preprocessing](#sca-process-msa---alignment-preprocessing)
3. [sca-core - Core SCA Calculations](#sca-core---core-sca-calculations)
4. [Complete Workflow Examples](#complete-workflow-examples)
5. [Command-Line Reference](#command-line-reference)
6. [Tips and Best Practices](#tips-and-best-practices)
7. [Troubleshooting](#troubleshooting)

---

## Quick Start

### Step 1: sca-process-msa

```bash
sca-process-msa \
    alignment.fasta \
    -s 1XYZ --chainID A \
    --precluster \
    --cluster-id 0.85 \
    --parameters 0.2 0.2 0.2 0.8 \
    --initial-trim-gap 0.8 \
    --output alignment_processed \
    --matlab \
    --log processing.log
```

**Note:** Use `--precluster` if MSA has >50k sequences. Use `--matlab` if you want a MATLAB database output.

### Step 2: sca-core

```bash
# Complete SCA analysis with sector identification
sca-core \
    Outputs/alignment_processed.db.gz \
    --do-seqcorr \
    --seqcorr-ref 0 \
    --seqcorr-mmseqs2 \
    --do-sector-id \
    --kpos 0 \
    --ic-cutoff 0.95 \
    --lbda 0.01 \
    --norm frob \
    --Ntrials 10 \
    --float32 \
    --matlab \
    --log sca_core.log
```

**Notes:**
- `--kpos 0`: Auto-selects eigenmodes (or specify eigenvalue cutoff)
- `--ic-cutoff 0.95`: >95% CDF of t-locationscale distribution for each independent component (IC)
- `--lbda 0.01`: Small regularization parameter, typically somewhat larger than 1/Meff

---

## sca-process-msa - Alignment Preprocessing

### Overview

Preprocesses multiple sequence alignments (MSAs) for SCA analysis. Filters sequences and positions, computes sequence weights, and prepares the alignment for downstream SCA calculations.

**Input:** Raw alignment file (FASTA/Stockholm/Clustal, optionally gzipped)  
**Output:** Processed alignment database (`.db.gz`) with sequence weights, filtered alignment, and metadata

### Basic Usage

```bash
sca-process-msa alignment.fasta
```

Alternatively, you can run the Python script directly:
```bash
python pysca/scaProcessMSA_py3_big.py alignment.fasta
```

### Input Requirements

- **Alignment file:** FASTA, Stockholm, or Clustal format (optionally `.gz` compressed)
- **Format detection:** Automatically detects format from file content

### Reference Sequence Options

You must specify a reference sequence using one of these methods:

#### Option 1: PDB Structure (Recommended)
```bash
sca-process-msa alignment.fasta \
    -s 1XYZ \
    -c A
```

Or using the long form:
```bash
sca-process-msa alignment.fasta \
    --pdb 1XYZ \
    --chainID A
```

**Note:** Both `-s`/`--pdb` and `-c`/`--chainID` are equivalent. The short forms (`-s`, `-c`) are more concise, while the long forms (`--pdb`, `--chainID`) are more explicit.

- Uses PDB structure to define reference sequence and alignment-to-structure mapping
- Creates distance matrix for structure-based analysis
- Optionally truncate to PDB positions: `--truncate`

#### Option 2: Reference Sequence File
```bash
sca-process-msa alignment.fasta \
    --refseq reference.fasta \
    --refpos positions.txt  # Optional: specific positions to map
```

#### Option 3: Reference Sequence Index
```bash
sca-process-msa alignment.fasta \
    --refindex 0  # Use first sequence (0-based index)
```

#### Option 4: Automatic Selection
```bash
sca-process-msa alignment.fasta
# Automatically chooses reference based on mean pairwise similarity
```

### Filtering Parameters

Default filtering parameters (can be customized):
```
--parameters max_gap_pos max_gap_seq min_SID max_SID
Default: 0.3 0.2 0.15 0.85
```

- **max_gap_pos:** Maximum fraction of gaps allowed per position (default: 0.3)
- **max_gap_seq:** Maximum fraction of gaps allowed per sequence (default: 0.2)
- **min_SID:** Minimum sequence identity threshold (default: 0.15)
- **max_SID:** Maximum sequence identity threshold (default: 0.85)

**Example:**
```bash
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --parameters 0.3 0.3 0.15 0.85
```

### Large MSA Optimization (MMseqs2 Preclustering)

For very large MSAs (>100k sequences), use MMseqs2 preclustering to reduce alignment size before processing:

```bash
sca-process-msa large_alignment.fasta \
    -s 1XYZ --chainID A \
    --precluster \
    --cluster-id 0.85 \
    --cluster-coverage 0.8
```

**Benefits:**
- Reduces alignment size (e.g., 300k → 50k sequences)
- Uses cluster-size weights (more accurate than computing weights on full alignment)
- Much faster processing
- Lower memory usage

**MMseqs2 Parameters:**
- `--cluster-id`: Sequence identity threshold (default: 0.85)
- `--cluster-coverage`: Coverage threshold (default: 0.8)
- `--cluster-cov-mode`: Coverage mode: 0=bidirectional, 1=query, 2=target (default: 0)

### Initial Position Trimming

Remove highly gapped positions before processing:

```bash
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
     --initial-trim-gap 0.8  # Remove positions with >80% gaps for initial reference sequence search
```

### Output Options

#### Default Output
- **Processed alignment:** `Inputs/<basename>_processed.fasta`
- **Database:** `Outputs/<basename>.db.gz`
- Contains: filtered alignment, sequence weights, ATS mapping, metadata

#### Save Numeric MSA
```bash
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --save-msa-num  # Saves numeric MSA to .npz file
```

#### MATLAB Workspace
```bash
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --matlab  # Also write .mat file
```

#### Custom Output Name
```bash
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --output my_analysis  # Output: my_analysis.db.gz
```

### Logging

```bash
# Verbose logging
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --verbose

# Quiet mode (warnings/errors only)
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --quiet

# Log to file
sca-process-msa alignment.fasta \
    -s 1XYZ --chainID A \
    --log processing.log
```

### Complete Example (Large MSA with PDB)

```bash
sca-process-msa \
    PF00071_alignment.fasta \
    -s 1XYZ --chainID A \
    --precluster \
    --cluster-id 0.85 \
    --parameters 0.2 0.2 0.2 0.8 \
    --initial-trim-gap 0.8 \
    --output PF00071_processed \
    --matlab \
    --log processing.log
```

**What this does:**
1. Loads alignment from `PF00071_alignment.fasta`
2. Uses PDB 1XYZ chain A as reference
3. Preclusters with MMseqs2 (85% identity threshold)
4. Filters sequences/positions
5. Computes sequence weights
6. Outputs database, MATLAB file, and log

---

## sca-core - Core SCA Calculations

### Overview

Performs core SCA calculations: site-specific conservation, SCA correlation matrix, and randomized controls. Optionally computes sequence correlations and sector identification.

**Input:** Processed alignment database from `scaProcessMSA_py3_big.py`  
**Output:** SCA results database (`.db.gz`) with correlation matrices, eigenvectors, and optional sector definitions

### Basic Usage

```bash
sca-core processed_alignment.db.gz
```

Alternatively, you can run the Python script directly:
```bash
python pysca/scaCore_py3.py processed_alignment.db.gz
```

### Input Requirements

- **Database file:** Must contain `db['sequence']` with:
  - `alg`: Filtered alignment
  - `seqw`: Sequence weights
  - `msa_num`: Numeric alignment (optional, computed if missing)
  - `ats`: Alignment-to-structure mapping (if PDB reference was used)

### Core SCA Calculations

These are always performed:

1. **Site-specific conservation ** (Di, Dia)
2. **SCA correlation matrix** (Csca)
3. **Projected alignment** (tX)
4. **Projector** (Proj)
5. **Randomized trials** (Vrand, Lrand)

### SCA Parameters

#### Regularization Parameter

Typically small for most MSAs; something a bit larger than 1/Meff, where Meff is the effective number of sequences in the MSA.
```bash
sca-core input.db.gz --lbda 0.01  # Default: 0.01
```

#### Matrix Norm Type
```bash
sca-core input.db.gz --norm frob   # Default: frob
# or
sca-core input.db.gz --norm spec   # Spectral norm
```

#### Randomization Trials
```bash
sca-core input.db.gz --Ntrials 20  # Default: 10
```

### Sequence Correlations (Optional)

Compute sequence similarity matrix and projections. **WARNING:** Creates O(M²) matrix - can be memory-intensive for large MSAs.

#### Enable for Small MSAs
```bash
sca-core input.db.gz --do-seqcorr
```

#### Automatic Subsampling (Recommended for Large MSAs)
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --seqcorr-ref 0  # Reference sequence to always retain
```

This automatically subsamples to **1.5 × effective sequences** (sum of weights), capped at 50k sequences by default.

#### Manual Subsampling
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --max-seqcorr-seqs 10000 \
    --seqcorr-ref 0
```

#### MMseqs2-Based Subsampling (Best Diversity)
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --max-seqcorr-seqs 10000 \
    --seqcorr-ref 0 \
    --seqcorr-mmseqs2
```

#### Keep Additional Sequences
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --seqcorr-ref 0 \
    --seqcorr-keep-indices 5 10 15  # Always retain these sequences
```

#### Disable Automatic Subsampling
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --no-auto-seqcorr-subsample  # Use all sequences (or --max-seqcorr-seqs if set)
```

#### Memory Cap
```bash
sca-core input.db.gz \
    --do-seqcorr \
    --seqcorr-max-cap 100000  # Increase cap from default 50000
```

### Sector Identification (Integrated)

Perform sector identification after SCA core calculations:

```bash
sca-core input.db.gz --do-sector-id
```

#### Sector ID Parameters

```bash
sca-core input.db.gz \
    --do-sector-id \
    --kpos 6 \              # Number of eigenmodes (0=auto)
    --ic-cutoff 0.95 \  # IC selection cutoff (default: 0.95 of CDF of a t-locationscale distribution)
    --kmax-cap 10           # Safety cap on kpos (default: 10)
```

**Automatic kpos selection:**
- If `--kpos 0` (default): Automatically selects based on Lrand eigenvalues
- Compares Csca eigenvalues to randomized MSA eigenvalues
- Selects eigenmodes where Csca eigenvalue > mean second eigenvalue of randomized MSA trials

### Memory Optimization

Use float32 instead of float64 to reduce memory:

```bash
sca-core input.db.gz --float32
```

**Benefits:**
- 50% memory reduction for large matrices
- Minimal impact on precision for most analyses

### Store Randomized Correlation Matrices

```bash
sca-core input.db.gz \
    --store-crand  # Store Crand (average randomized correlation matrix)
```

**Note:** Increases database size significantly. Only use if needed.

### Output Options

#### Default Output
- **Database:** `Outputs/<basename>.db.gz`
- Contains: SCA results in `db['sca']`, optional sector results in `db['sector']`

#### MATLAB Workspace
```bash
sca-core input.db.gz --matlab
```

#### Custom Output Name
```bash
sca-core input.db.gz \
    --output my_sca_results  # Output: my_sca_results.db.gz
```

### Logging

```bash
# Verbose logging
sca-core input.db.gz --verbose

# Quiet mode
sca-core input.db.gz --quiet

# Log to file
sca-core input.db.gz --log sca_core.log
```

### Complete Example (Full Workflow)

```bash
# Complete SCA analysis with sector identification
sca-core \
    PF00071_processed.db.gz \
    --do-seqcorr \
    --seqcorr-ref 0 \
    --seqcorr-mmseqs2 \
    --do-sector-id \
    --kpos 0 \
    --ic-cutoff 0.95 \
    --lbda 0.03 \
    --norm frob \
    --Ntrials 10 \
    --float32 \
    --matlab \
    --log sca_core.log
```

**What this does:**
1. Loads processed alignment database
2. Computes site-specific conservation and SCA matrix
3. Computes randomized controls (10 trials)
4. Computes sequence correlations with MMseqs2 subsampling
5. Performs sector identification with automatic kpos selection
6. Saves results with float32 precision
7. Outputs MATLAB workspace

---

## Complete Workflow Examples

### Example 1: Standard Workflow (Small-Medium MSA)

```bash
# Step 1: Preprocess alignment
sca-process-msa \
    alignment.fasta \
    -s 1XYZ --chainID A \
    --output my_analysis

# Step 2: Run SCA core calculations
sca-core \
    Outputs/my_analysis.db.gz \
    --do-seqcorr \
    --do-sector-id

# Results in: Outputs/my_analysis.db.gz
# Contains: sequence, sca, and sector data
```

### Example 2: Large MSA Workflow (300k+ sequences)

```bash
# Step 1: Preprocess with MMseqs2 preclustering
sca-process-msa \
    large_alignment.fasta \
    -s 1XYZ --chainID A \
    --precluster \
    --cluster-id 0.85 \
    --output large_analysis

# Step 2: Run SCA core with subsampling
sca-core \
    Outputs/large_analysis.db.gz \
    --do-seqcorr \
    --seqcorr-mmseqs2 \
    --seqcorr-max-cap 50000 \
    --do-sector-id \
    --float32

# Results efficiently processed even for 300k sequences
```

### Example 3: High-Precision Analysis (All Sequences)

```bash
# Step 1: Preprocess (no preclustering for maximum accuracy)
sca-process-msa \
    alignment.fasta \
    -s 1XYZ --chainID A \
    --output high_precision

# Step 2: Run SCA core without subsampling
sca-core \
    Outputs/high_precision.db.gz \
    --do-seqcorr \
    --no-auto-seqcorr-subsample \
    --do-sector-id \
    --kpos 10 \
    --Ntrials 20

# Uses all sequences (may be slow/memory-intensive for large MSAs)
```

### Example 4: Quick Analysis (Minimal Parameters)

```bash
# Step 1: Minimal preprocessing
sca-process-msa alignment.fasta -s 1XYZ --chainID A

# Step 2: Minimal SCA core (no sequence correlations, no sector ID)
sca-core Outputs/alignment.db.gz

# Fastest workflow, core SCA results only
```

---

## Command-Line Reference

### sca-process-msa - Alignment Preprocessing

#### Required
- `alignment`: Input alignment file (FASTA/Stockholm/Clustal, optionally .gz)

#### Reference Sequence (choose one)
- `-s PDBID` or `--pdb PDBID`: PDB identifier or path (**Note:** Use `-s`, not `-pdb`, since `-p` is reserved for `--parameters`)
- `--chainID CHAIN`: Chain ID (default: A)
- `--refseq FILE`: Reference sequence FASTA file
- `--refindex INT`: Reference sequence index (0-based)
- (none): Automatic selection

#### Filtering
- `--parameters FLOAT FLOAT FLOAT FLOAT`: [max_gap_pos, max_gap_seq, min_SID, max_SID] (default: 0.3 0.2 0.15 0.85)
- `--initial-trim-gap FLOAT`: Initial gap trim threshold (default: 0.8)
- `--truncate`: Truncate to PDB positions (if using PDB)

#### Large MSA Optimization
- `--precluster`: Enable MMseqs2 preclustering
- `--cluster-id FLOAT`: MMseqs2 identity threshold (default: 0.85)
- `--cluster-coverage FLOAT`: MMseqs2 coverage threshold (default: 0.8)
- `--cluster-cov-mode INT`: Coverage mode: 0=bidirectional, 1=query, 2=target (default: 0)
- `--keep-mmseqs-tmp`: Keep MMseqs2 temporary files (for debugging)

#### Output
- `--output NAME`: Output base name (default: derived from input filename)
- `--matlab`: Write MATLAB workspace (.mat file)
- `--save-msa-num`: Save numeric MSA to .npz file

#### Logging
- `--log FILE`: Write log to file
- `--verbose`: Verbose console logging
- `--quiet`: Quiet mode (warnings/errors only)

### sca-core - Core SCA Calculations

#### Required
- `database`: Input database file (.db or .db.gz) from scaProcessMSA

#### Core Parameters
- `--lbda FLOAT`: Regularization parameter (default: 0.01)
- `--norm {frob,spec}`: Matrix norm type (default: frob)
- `--Ntrials INT`: Number of randomization trials (default: 10)

#### Sequence Correlations (Optional)
- `--do-seqcorr`: Enable sequence correlation calculations
- `--max-seqcorr-M INT`: Safety limit, skip seqcorr if M exceeds this (default: 5000)
- `--max-seqcorr-seqs INT`: Manual max sequences for subsampling
- `--seqcorr-ref INT`: Reference sequence index (always retained)
- `--seqcorr-mmseqs2`: Use MMseqs2 clustering for subsampling
- `--seqcorr-keep-indices INT [INT ...]`: Additional sequences to always retain
- `--seqcorr-max-cap INT`: Maximum sequences cap (default: 50000)
- `--no-auto-seqcorr-subsample`: Disable automatic 1.5×M_eff subsampling

#### Sector Identification (Optional)
- `--do-sector-id`: Enable sector identification
- `--kpos INT`: Number of eigenmodes (0=auto, default: 0)
- `--ic-cutoff FLOAT`: IC selection cutoff (default: 0.95)
- `--kmax-cap INT`: Safety cap on kpos (default: 10)

#### Memory & Storage
- `--float32`: Use float32 instead of float64 (saves memory)
- `--store-crand`: Store randomized correlation matrices (increases DB size)

#### Output
- `--output NAME`: Output base name (default: derived from input filename)
- `--matlab`: Write MATLAB workspace (.mat file)

#### Logging
- `--log FILE`: Write log to file
- `--verbose`: Verbose console logging
- `--quiet`: Quiet mode (warnings/errors only)

---

## Tips and Best Practices

### For Small MSAs (<10k sequences)
- Don't need preclustering or subsampling
- Can enable all optional calculations (seqcorr, sector ID)
- Default settings work well

### For Medium MSAs (10k-100k sequences)
- Consider MMseqs2 preclustering for faster processing
- Use automatic seqSim subsampling (1.5×M_eff)
- Enable sector ID

### For Large MSAs (>100k sequences)
- **Always use MMseqs2 preclustering** in scaProcessMSA
- Use automatic seqSim subsampling in scaCore
- Consider `--float32` to reduce memory
- May need to adjust `--seqcorr-max-cap` based on available memory

### Memory Considerations
- Sequence correlation matrix (simMat) is O(M²) - be careful with large M
- Automatic subsampling (1.5×M_eff) is recommended
- Use `--float32` for 50% memory savings
- Monitor memory usage with system tools

### Performance Tips
- Use MMseqs2 preclustering for very large MSAs (saves hours)
- Use `--quiet` mode for batch processing
- Save logs for debugging with `--log`
- Use `--save-msa-num` only if you need numeric MSA separately

---

## Troubleshooting

### "Out of memory" errors
- Use `--precluster` in scaProcessMSA
- Use automatic subsampling in scaCore (default behavior)
- Reduce `--seqcorr-max-cap`
- Use `--float32`
- Reduce `--Ntrials`

### "ggsearch36 not found" or "ggsearch36 failed" errors

These errors occur when `sca-process-msa` tries to find a reference sequence in the alignment but `ggsearch36` is not installed or not in your PATH.

**Symptoms:**
- `RuntimeError: ggsearch36 failed: ... Make sure ggsearch36 is installed and in your PATH.`
- `FileNotFoundError: ggsearch36 command not found`
- Error occurs when using `-s` (PDB reference) or `--refseq` options

**Solutions:**

1. **Install FASTA36** (see [INSTALLATION.md](INSTALLATION.md))

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

These errors occur when using `--precluster` option but MMseqs2 is not installed.

**Solutions:**
- Install MMseqs2: `conda install -c bioconda mmseqs2` or download from https://github.com/soedinglab/MMseqs2
- Ensure MMseqs2 is in PATH: `which mmseqs`
- See [INSTALLATION.md](INSTALLATION.md) for detailed installation instructions
- If MMseqs2 is not available, you can disable preclustering with `--no-precluster` (not recommended for large MSAs)

### Slow processing
- Enable MMseqs2 preclustering for large MSAs
- Reduce randomization trials: `--Ntrials 5`
- Use `--quiet` to reduce I/O overhead

### Sector ID issues
- Check that Lrand is present in database (from randomization trials)
- Try explicit `--kpos` instead of automatic selection
- Adjust `--ic-cutoff` if too few/many sectors identified

---

## Output Database Structure

### After scaProcessMSA_py3_big.py
```python
db = {
    'sequence': {
        # Final processed alignment
        'alg': list[str],        # Filtered alignment (final processed)
        'hd': list[str],         # Headers (final processed)
        'seqw': np.ndarray,      # Sequence weights
        'Nseq': int,             # Number of sequences
        'Npos': int,             # Number of positions
        'ats': list,             # Alignment-to-structure mapping
        'effseqs': float,        # Effective number of sequences
        
        # Reference sequence information
        'i_ref': int,            # Reference sequence index in final processed alignment (alg)
        'i_ref_original': int,   # Reference sequence index in original input MSA (alg_original)
        'i_ref_was_moved': bool, # Flag: True if reference was moved to index 0 during filtering
        'ref_header_id': str,    # Reference sequence header ID (for finding in sca-core)
        
        # Original input MSA (before filtering, after preclustering if used)
        'alg_original': list[str],  # Original input MSA sequences (or preclustered if preclustering was done)
        'hd_original': list[str],   # Original input MSA headers (or preclustered if preclustering was done)
        
        # Additional metadata
        'preclustered': bool,        # True if MMseqs2 preclustering was used
        'cluster_id': float,         # MMseqs2 cluster identity threshold (if preclustered)
        'trim_parameters': list,     # Filtering parameters [pos_gap, seq_gap, min_seqid, max_seqid]
        'truncate_flag': bool,       # True if truncation was applied
        ...
    }
}
```

**Important notes:**
- `i_ref`: Position of reference sequence in the **final processed alignment** (`alg`). This is the index you use to access the reference in the filtered alignment.
- `i_ref_original`: Position of reference sequence in the **original input MSA** (`alg_original`). Use this with `alg_original`/`hd_original` to access the reference sequence in the original alignment.
- `i_ref_was_moved`: Set to `True` if the reference sequence was filtered out and then added back at index 0. In this case, `i_ref` will be 0, but `i_ref_original` still points to the original position.
- `alg_original`/`hd_original`: Store the original input MSA (or preclustered MSA if preclustering was done). This is the alignment that `i_ref_original` refers to.

### After scaCore_py3.py
```python
db = {
    'sequence': {...},  # From scaProcessMSA
    'sca': {
        'Di': np.ndarray,        # Position conservation
        'Dia': np.ndarray,       # Position-aa conservation
        'Csca': np.ndarray,      # SCA correlation matrix
        'tX': np.ndarray,        # Projected alignment
        'Proj': np.ndarray,      # Projector
        'Vrand': np.ndarray,     # Randomized eigenvectors
        'Lrand': np.ndarray,     # Randomized eigenvalues
        'simMat': np.ndarray,    # Sequence similarity (if --do-seqcorr)
        'Useq': np.ndarray,      # Sequence projections (if --do-seqcorr)
        'Uica': np.ndarray,      # Sequence ICA (if --do-seqcorr)
        ...
    },
    'sector': {                  # If --do-sector-id
        'kpos': int,             # Number of eigenmodes
        'ic_cutoff': float,      # IC selection cutoff parameter
        'ic_pos': list,          # IC positions (0-based indices)
        'ic_ats': list,          # IC ATS labels
        'Vica': np.ndarray,      # ICA-rotated eigenvectors
        ...
    }
}
```

---

For more details, see:
- [INSTALLATION.md](INSTALLATION.md) - Installation and setup
- [SCATOOLS_FUNCTIONS.md](SCATOOLS_FUNCTIONS.md) - Function reference
- [OPTIMIZATIONS.md](OPTIMIZATIONS.md) - Performance optimizations

Or run:
```bash
sca-process-msa --help
sca-core --help
```

**Note:** These commands (`sca-process-msa` and `sca-core`) are wrapper scripts that execute the underlying Python 3 scripts. After installation with `pip install -e .`, they will be available in your PATH. If the commands are not found, you can run the Python scripts directly:
```bash
python pysca/scaProcessMSA_py3_big.py --help
python pysca/scaCore_py3.py --help
```

