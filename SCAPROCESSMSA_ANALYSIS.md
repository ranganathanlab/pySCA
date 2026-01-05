# scaProcessMSA_py3_big.py Analysis and Optimization Plan

## Workflow Overview

The script performs the following steps:

1. **MMseqs2 Preclustering** (optional) - Reduces alignment size
2. **Load Alignment** - Reads FASTA/Stockholm/Clustal
3. **Filter Non-Standard AAs** - Removes invalid sequences
4. **Initial Position Trim** - Removes highly gapped positions
5. **Reference Sequence Determination** - PDB/refseq/auto-select
6. **Alignment-to-Structure Mapping** - makeATS()
7. **Sequence Filtering** - By gap fraction and identity
8. **Position Filtering** - By gap fraction
9. **Weight Computation** - MMseqs2 weights or seqWeights()
10. **Output Writing** - FASTA, pickle DB, optional MATLAB

---

## Identified Bottlenecks

### 1. `filter_nonstandard()` - Line 89-96
**Current**: Iterates through sequences, converts to uppercase, checks set membership
**Issue**: O(N×L) operations with Python loops
**Optimization**: Vectorized NumPy operations

### 2. `write_fasta()` - Line 99-103
**Current**: Writes sequences one by one
**Issue**: Multiple file writes, could be batched
**Optimization**: Batch writing or use BioPython SeqIO

### 3. Distance Matrix Remapping - Line 257-265
**Current**: Uses list comprehension and dictionary lookup
**Issue**: Could be more efficient with NumPy operations
**Optimization**: Vectorized NumPy indexing

### 4. Multiple `readAlg()` calls
**Current**: Called multiple times (lines 215, 269)
**Issue**: Redundant file I/O
**Optimization**: Cache results when possible

### 5. `load_weights_tsv()` - Line 77-86
**Current**: Simple file reading
**Issue**: Could use more efficient parsing
**Optimization**: Already efficient, but could add validation

### 6. Memory Usage
**Current**: Creates multiple intermediate lists
**Issue**: For 300k sequences, this uses significant memory
**Optimization**: Use generators where possible, avoid unnecessary copies

---

## Optimization Strategy


