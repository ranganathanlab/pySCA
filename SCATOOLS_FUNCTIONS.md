# scaTools.py Function Reference

Comprehensive guide to the core functions in `scaTools.py`, the central module for Statistical Coupling Analysis (SCA) calculations.

---

## Table of Contents

1. [Alignment I/O Functions](#alignment-io-functions)
2. [Sequence Processing Functions](#sequence-processing-functions)
3. [Statistical Functions](#statistical-functions)
4. [SCA Core Functions](#sca-core-functions)
5. [Independent Component Analysis](#independent-component-analysis)
6. [Utility Functions](#utility-functions)

---

## Alignment I/O Functions

### `readAlg(filename)`

Read a multiple sequence alignment from file.

**Supported formats:** FASTA, Stockholm, Clustal (auto-detected)

**Parameters:**
- `filename` (str): Path to alignment file (optionally gzipped)

**Returns:**
- `headers` (list[str]): Sequence headers
- `sequences` (list[str]): Sequences (uppercased, gaps normalized)

**Example:**
```python
headers, sequences = sca.readAlg("alignment.fasta")
```

**Notes:**
- Automatically detects format from file content
- Uses BioPython AlignIO for Stockholm/Clustal when available
- Falls back to minimal parsers if BioPython fails
- Converts Stockholm '.' gaps to '-' for consistency

---

## Sequence Processing Functions

### `lett2num(msa_lett, code="ACDEFGHIKLMNPQRSTVWY")`

Convert letter-based alignment to numeric representation.

**Parameters:**
- `msa_lett` (list[str] or np.ndarray): Letter-based alignment
- `code` (str): Amino acid code (default: standard 20 AAs)

**Returns:**
- `msa_num` (np.ndarray): Numeric alignment (0-19 for AAs, 20 for gap)

**Example:**
```python
msa_num = sca.lett2num(["ACDEFG", "ACD-EF"])
# Returns: array([[0, 1, 2, 3, 4, 5], [0, 1, 2, 20, 3, 4]])
```

**Optimization:** Uses vectorized NumPy lookup table for fast conversion.

---

### `alg2bin(alg, N_aa=20)`

Convert numeric alignment to binary one-hot representation.

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `N_aa` (int): Number of amino acids (default: 20)

**Returns:**
- `Abin` (sparse matrix): Binary matrix (M × (N_aa × L))

**Example:**
```python
Abin = sca.alg2bin(msa_num)
# Each position becomes N_aa columns (one-hot encoding)
```

**Optimization:** Uses advanced NumPy indexing for efficient sparse matrix construction.

---

### `alg2binss(alg, N_aa=20)`

Convert numeric alignment to sparse binary representation (optimized version).

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `N_aa` (int): Number of amino acids (default: 20)

**Returns:**
- `Abin` (sparse CSR matrix): Sparse binary matrix (M × (N_aa × L))

**Notes:**
- More memory-efficient than `alg2bin` for large alignments
- Directly constructs sparse matrix without dense intermediates
- Recommended for alignments with >10k sequences

---

### `filterSeq(alg0, sref=0.5, max_fracgaps=0.2, min_seqid=0.2, max_seqid=0.8)`

Filter sequences by gap fraction and sequence identity.

**Parameters:**
- `alg0` (list[str]): Alignment sequences
- `sref` (int or float): Reference sequence index (for identity calculation)
- `max_fracgaps` (float): Maximum gap fraction per sequence (default: 0.2)
- `min_seqid` (float): Minimum sequence identity to reference (default: 0.2)
- `max_seqid` (float): Maximum sequence identity to reference (default: 0.8)

**Returns:**
- `alg_filtered` (list[str]): Filtered alignment
- `seqw` (np.ndarray): Sequence weights
- `seqkeep` (list[int]): Indices of kept sequences

**Example:**
```python
alg_filtered, seqw, seqkeep = sca.filterSeq(
    alg, sref=0, max_fracgaps=0.2, min_seqid=0.2, max_seqid=0.8
)
```

**Optimization:** Vectorized gap counting and sequence identity calculations using NumPy.

---

### `filterPos(alg, seqw=[1], max_fracgaps=0.2)`

Filter positions by gap fraction.

**Parameters:**
- `alg` (list[str]): Alignment sequences
- `seqw` (array-like): Sequence weights
- `max_fracgaps` (float): Maximum gap fraction per position (default: 0.2)

**Returns:**
- `alg_filtered` (list[str]): Filtered alignment
- `poskeep` (list[int]): Indices of kept positions

**Example:**
```python
alg_filtered, poskeep = sca.filterPos(alg, seqw, max_fracgaps=0.2)
```

**Optimization:** Vectorized gap counting using NumPy.

---

### `seqWeights(alg, max_seqid=0.8, gaps=1, block_size=1024)`

Compute sequence weights based on sequence identity.

**Parameters:**
- `alg` (list[str]): Alignment sequences
- `max_seqid` (float): Maximum sequence identity threshold (default: 0.8)
- `gaps` (int): Whether to include gaps in identity calculation (default: 1)
- `block_size` (int): Block size for large MSAs (default: 1024)

**Returns:**
- `seqw` (np.ndarray): Sequence weights (1 × M)

**Example:**
```python
seqw = sca.seqWeights(alg, max_seqid=0.8)
```

**Optimization:** Blockwise similarity counting for large MSAs to reduce memory usage.

---

### `MSAsearch(hd, algn, seq, species=None)`

Find the best matching sequence in alignment for a reference sequence.

**Parameters:**
- `hd` (list[str]): Alignment headers
- `algn` (list[str]): Alignment sequences
- `seq` (str): Reference sequence to search for
- `species` (str, optional): Species filter to narrow search

**Returns:**
- `index` (int): Index of best matching sequence

**Example:**
```python
i_ref = sca.MSAsearch(headers, sequences, pdb_sequence, species="Homo sapiens")
```

**Notes:**
- Uses `ggsearch36` (FASTA36) for sequence alignment
- Returns detailed alignment statistics in log
- Automatically filters by species if provided

---

### `chooseRefSeq(alg)`

Automatically choose a reference sequence based on mean pairwise similarity.

**Parameters:**
- `alg` (list[str]): Alignment sequences

**Returns:**
- `index` (int): Index of chosen reference sequence

**Example:**
```python
i_ref = sca.chooseRefSeq(alg)
```

**Optimization:** Uses blockwise dot products for large alignments to control memory.

---

### `makeATS(sequences, refpos, refseq, iref=0, truncate=False)`

Create Alignment-To-Structure (ATS) mapping.

**Parameters:**
- `sequences` (list[str]): Alignment sequences
- `refpos` (list): Reference positions (e.g., PDB residue numbers)
- `refseq` (str): Reference sequence
- `iref` (int): Reference sequence index (default: 0)
- `truncate` (bool): Truncate to reference positions (default: False)

**Returns:**
- `sequences_ats` (list[str]): Sequences with ATS mapping
- `ats` (list): ATS position labels

**Example:**
```python
sequences_ats, ats = sca.makeATS(sequences, pdb_positions, pdb_sequence, iref=0)
```

---

## Statistical Functions

### `freq(alg, Naa=20, seqw=None, lbda=0, freq0=None)`

Compute weighted single-site and pairwise amino-acid frequencies.

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `Naa` (int): Number of amino acids (default: 20)
- `seqw` (array-like): Sequence weights (M,) or (1, M)
- `lbda` (float): Regularization parameter (default: 0)
- `freq0` (array-like, optional): Background frequencies

**Returns:**
- `freq1` (np.ndarray): Single-site frequencies (N_aa × L)
- `freq2` (np.ndarray): Pairwise frequencies (N_aa × L × N_aa × L)
- `freq0` (np.ndarray): Background frequencies used

**Example:**
```python
freq1, freq2, freq0 = sca.freq(msa_num, seqw=seqw, lbda=0.03)
```

**Optimization:** Keeps computations sparse for large MSAs, only converts to dense at the end.

---

### `posWeights(alg, seqw=1, lbda=0, N_aa=20, freq0=None, tolerance=1e-12)`

Compute positional weights (Di, Dia) from alignment.

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `seqw` (array-like): Sequence weights
- `lbda` (float): Regularization parameter
- `N_aa` (int): Number of amino acids (default: 20)
- `freq0` (array-like, optional): Background frequencies
- `tolerance` (float): Numerical tolerance (default: 1e-12)

**Returns:**
- `Wia` (np.ndarray): Position-AA weights (N_aa × L)
- `Dia` (np.ndarray): Position-AA conservation (N_aa × L, flattened)
- `Di` (np.ndarray): Position conservation (L,)

**Example:**
```python
Wia, Dia, Di = sca.posWeights(msa_num, seqw=seqw, lbda=0.03)
```

**Notes:**
- `Dia` is a 1D array of shape `(L * N_aa,)` containing position-AA conservation values
- `Di` is a 1D array of shape `(L,)` containing position conservation values

---

## SCA Core Functions

### `scaMat(alg, seqw=1, norm="frob", lbda=0, freq0=np.ones(20) / 21)`

Compute the SCA correlation matrix.

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `seqw` (array-like): Sequence weights (default: uniform)
- `norm` (str): Matrix norm type - 'frob' or 'spec' (default: 'frob')
- `lbda` (float): Regularization parameter (default: 0)
- `freq0` (array-like): Background frequencies

**Returns:**
- `Csca` (np.ndarray): SCA correlation matrix (L × L)
- `tX` (np.ndarray): Projected alignment (M × L)
- `Proj` (np.ndarray): Projector matrix

**Example:**
```python
Csca, tX, Proj = sca.scaMat(msa_num, seqw=seqw, norm='frob', lbda=0.03)
```

**Optimization:** Keeps sparse operations for projection computation, handles both sparse and dense results.

---

### `seqSim(alg, max_seqs=None, i_ref=None, use_mmseqs2=False, cluster_id=0.85, seqw=None, return_indices=False, max_seqs_cap=50000, keep_indices=None, auto_subsample=True)`

Compute sequence similarity matrix with optional subsampling.

**Parameters:**
- `alg` (np.ndarray): Numeric alignment (M × L)
- `max_seqs` (int, optional): Maximum sequences to use (None = auto or all)
- `i_ref` (int, optional): Reference sequence index (always retained)
- `use_mmseqs2` (bool): Use MMseqs2 for subsampling (default: False)
- `cluster_id` (float): MMseqs2 identity threshold (default: 0.85)
- `seqw` (array-like): Sequence weights for weighted selection
- `return_indices` (bool): Return selected indices (default: False)
- `max_seqs_cap` (int): Maximum sequences cap (default: 50000)
- `keep_indices` (list[int]): Additional sequences to always retain
- `auto_subsample` (bool): Auto-subsample to 1.5×M_eff (default: True)

**Returns:**
- `simMat` (np.ndarray): Sequence similarity matrix (M' × M')
- `selected_indices` (np.ndarray, optional): Indices of selected sequences

**Example:**
```python
# Automatic subsampling
simMat = sca.seqSim(msa_num, seqw=seqw, i_ref=0)

# Manual subsampling
simMat, indices = sca.seqSim(msa_num, max_seqs=10000, i_ref=0, return_indices=True)

# MMseqs2-based subsampling
simMat = sca.seqSim(msa_num, seqw=seqw, i_ref=0, use_mmseqs2=True)
```

**Optimization:** Intelligent subsampling for large MSAs, blockwise computation for memory efficiency.

---

### `seqProj(msa_num, seqw, kseq=15, kica=6)`

Compute sequence projections based on eigenvectors of similarity matrix.

**Parameters:**
- `msa_num` (np.ndarray): Numeric alignment (M × L)
- `seqw` (array-like): Sequence weights (1 × M)
- `kseq` (int): Number of eigenvectors (default: 15)
- `kica` (int): Number of ICA components (default: 6)

**Returns:**
- `Useq` (list): Sequence projections (3 variants: raw, weighted, weighted+positional)
- `Uica` (list): ICA-rotated projections (3 variants)

**Example:**
```python
Useq, Uica = sca.seqProj(msa_num, seqw, kseq=30, kica=15)
```

---

### `eigenVect(M)`

Compute eigenvectors and eigenvalues of symmetric matrix.

**Parameters:**
- `M` (np.ndarray): Symmetric matrix

**Returns:**
- `eigenVectors` (np.ndarray): Eigenvectors (sorted by eigenvalue)
- `eigenValues` (np.ndarray): Eigenvalues (sorted descending)

**Example:**
```python
V, L = sca.eigenVect(Csca)
```

---

### `randomize(msa_num, Ntrials, seqw, norm="frob", lbda=0, Naa=20, kmax=50, tolerance=1e-12)`

Generate randomized alignments and compute SCA matrices for statistical comparison.

**Parameters:**
- `msa_num` (np.ndarray): Numeric alignment (M × L)
- `Ntrials` (int): Number of randomization trials
- `seqw` (array-like): Sequence weights
- `norm` (str): Matrix norm type (default: 'frob')
- `lbda` (float): Regularization parameter (default: 0)
- `Naa` (int): Number of amino acids (default: 20)
- `kmax` (int): Maximum eigenmodes to compute (default: 50)
- `tolerance` (float): Numerical tolerance (default: 1e-12)

**Returns:**
- `Vrand` (np.ndarray): Randomized eigenvectors (Ntrials × L × kmax)
- `Lrand` (np.ndarray): Randomized eigenvalues (Ntrials × L)
- `Crand` (np.ndarray, optional): Average randomized correlation matrix

**Example:**
```python
Vrand, Lrand, Crand = sca.randomize(msa_num, Ntrials=10, seqw=seqw, norm='frob', lbda=0.03)
```

---

## Independent Component Analysis

### `rotICA(V, kmax=6, learnrate=0.1, iterations=100000)`

Perform Independent Component Analysis (ICA) rotation on eigenvectors.

**Parameters:**
- `V` (np.ndarray): Eigenvectors (L × kpos)
- `kmax` (int): Number of independent components (default: 6)
- `learnrate` (float): Learning rate (default: 0.1)
- `iterations` (int): Maximum iterations (default: 100000)

**Returns:**
- `Vica` (np.ndarray): ICA-rotated eigenvectors (L × kmax)
- `W` (np.ndarray): Unmixing matrix

**Example:**
```python
Vica, W = sca.rotICA(V, kmax=6)
```

---

### `icList(Vpica, kpos, Csca, p_cut=0.95)`

Identify significant positions for each independent component using t-distribution fits.

**Parameters:**
- `Vpica` (np.ndarray): ICA-rotated eigenvectors (L × kpos)
- `kpos` (int): Number of eigenmodes/components
- `Csca` (np.ndarray): SCA correlation matrix (L × L)
- `p_cut` (float): Cutoff on CDF (default: 0.95)

**Returns:**
- `ic_list` (list[Unit]): List of Unit objects, each containing significant positions
- `icsize` (list[int]): Size of each IC
- `sortedpos` (list[list[int]]): Sorted positions per IC
- `cutoff` (list[float]): Cutoff values per IC
- `scaled_pdf` (list[np.ndarray]): Scaled PDF fits per IC
- `all_fits` (list[tuple]): T-distribution fit parameters (df, loc, scale) per IC

**Example:**
```python
ic_list, icsize, sortedpos, cutoff, scaled_pdf, all_fits = sca.icList(Vica, kpos=6, Csca, p_cut=0.95)
```

**Notes:**
- Fits t-distribution to histogram of each IC
- Assigns positions above cutoff to ICs
- Resolves conflicts when positions appear in multiple ICs

---

### `t(Vica, ic_list)`

Extract position indices from IC list (wrapper function).

**Parameters:**
- `Vica` (np.ndarray): ICA-rotated eigenvectors
- `ic_list` (list[Unit]): List of Unit objects from `icList`

**Returns:**
- `sectors` (list[list[int]]): List of position index lists

**Example:**
```python
sectors = sca.t(Vica, ic_list)
```

---

### `chooseKpos(Lsca, Lrand)`

Automatically choose number of significant eigenmodes based on randomized controls.

**Parameters:**
- `Lsca` (np.ndarray): Eigenvalues of SCA matrix
- `Lrand` (np.ndarray): Eigenvalues from randomized trials (Ntrials × L)

**Returns:**
- `kpos` (int): Number of significant eigenmodes

**Example:**
```python
kpos = sca.chooseKpos(Lsca, Lrand)
```

---

## Utility Functions

### `pdbSeq(pdbid, chain="A", path2pdb=settings.path2structures, calcDist=1)`

Extract sequence and distance matrix from PDB structure.

**Parameters:**
- `pdbid` (str): PDB identifier or path
- `chain` (str): Chain ID (default: 'A')
- `path2pdb` (str): Path to PDB structures directory
- `calcDist` (int): Whether to calculate distance matrix (default: 1)

**Returns:**
- `seq` (str): Sequence from PDB
- `ats` (list): ATS position labels
- `distmat` (np.ndarray): Distance matrix (if calcDist=1)

**Example:**
```python
seq, ats, distmat = sca.pdbSeq("1XYZ", chain="A")
```

---

### `clean_al(alg, code="ACDEFGHIKLMNPQRSTVWY", gap="-")`

Clean alignment by removing non-standard characters.

**Parameters:**
- `alg` (list[str]): Alignment sequences
- `code` (str): Allowed amino acid code
- `gap` (str): Gap character

**Returns:**
- `alg_clean` (list[str]): Cleaned alignment

**Example:**
```python
alg_clean = sca.clean_al(alg)
```

**Optimization:** Uses list comprehensions and `str.join` for faster string building.

---

### `svdss(X, k=6)`

Singular value decomposition for sparse matrices (top k components).

**Parameters:**
- `X` (sparse matrix): Input matrix
- `k` (int): Number of components (default: 6)

**Returns:**
- `u` (np.ndarray): Left singular vectors
- `s` (np.ndarray): Singular values
- `v` (np.ndarray): Right singular vectors

**Example:**
```python
u, s, v = sca.svdss(X, k=10)
```

---

## Function Categories Summary

### Core SCA Workflow

1. **Read alignment:** `readAlg()`
2. **Convert to numeric:** `lett2num()`
3. **Filter sequences:** `filterSeq()`
4. **Filter positions:** `filterPos()`
5. **Compute weights:** `seqWeights()` or use MMseqs2 cluster sizes
6. **Compute frequencies:** `freq()`
7. **Compute positional weights:** `posWeights()`
8. **Compute SCA matrix:** `scaMat()`
9. **Eigendecomposition:** `eigenVect()`
10. **ICA rotation:** `rotICA()`
11. **Identify ICs:** `icList()`

### Sequence Analysis

- `seqSim()`: Sequence similarity matrix
- `seqProj()`: Sequence projections
- `chooseRefSeq()`: Automatic reference selection
- `MSAsearch()`: Find best matching sequence

### Statistical Analysis

- `randomize()`: Generate randomized controls
- `chooseKpos()`: Auto-select significant eigenmodes

---

For detailed function signatures and examples, see the docstrings in `pysca/scaTools.py`.


