# Multiple Sequence Alignment Format Support

## Summary

The `readAlg()` function has been extended to support **FASTA**, **Stockholm**, and **Clustal** alignment formats with automatic format detection.

---

## Usage

### Automatic Format Detection (Recommended)

```python
from pysca import scaTools as sca

# Auto-detect format
headers, sequences = sca.readAlg('alignment.fasta')
headers, sequences = sca.readAlg('alignment.stockholm')
headers, sequences = sca.readAlg('alignment.clustal')
```

### Explicit Format Specification

```python
# Explicitly specify format (bypasses auto-detection)
headers, sequences = sca.readAlg('alignment.txt', format='fasta')
headers, sequences = sca.readAlg('alignment.txt', format='stockholm')
headers, sequences = sca.readAlg('alignment.txt', format='clustal')
```

---

## Supported Formats

### 1. FASTA Format

**Features:**
- Sequence headers start with `>`
- Sequence data follows header line
- Multi-line sequences supported
- Most common format

**Example:**
```
>seq1
ACDEFGHIKLMNPQRSTVWY
>seq2
ACDEFGHIKLMNPQRSTVWY
```

### 2. Stockholm Format

**Features:**
- Header line: `# STOCKHOLM 1.0`
- Sequence lines: `sequence_name    sequence_data` (whitespace separated)
- Annotation lines start with `#`
- Ends with `//`
- Used by Pfam database

**Example:**
```
# STOCKHOLM 1.0
seq1    ACDEFGHIKLMNPQRSTVWY
seq2    ACDEFGHIKLMNPQRSTVWY
#=GC RF    xxxxxxxxxxxxxxxxxxxx
//
```

### 3. Clustal Format

**Features:**
- Header: `CLUSTAL` or `CLUSTALW` or `CLUSTAL W`
- Sequence blocks separated by blank lines
- Each line: `sequence_name    sequence_data`
- Conservation line at bottom of each block (with `*`, `:`, `.` symbols)
- Used by ClustalW/ClustalX alignment tools

**Example:**
```
CLUSTAL W (1.81)

seq1    ACDEFGHIKLMNPQRSTVWY
seq2    ACDEFGHIKLMNPQRSTVWY
        *:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*

seq1    ACDEFGHIKLMNPQRSTVWY
seq2    ACDEFGHIKLMNPQRSTVWY
        *:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*:.*
```

---

## Format Detection Algorithm

The auto-detection function (`_detect_alignment_format()`) examines the first 20 lines of the file:

1. **Stockholm detection:**
   - Looks for `# STOCKHOLM` header
   - Or detects `//` terminator with annotation lines (`#`)

2. **Clustal detection:**
   - Looks for `CLUSTAL` header (case-insensitive)
   - Or detects conservation lines (spaces + `*`, `:`, `.` symbols)

3. **FASTA (default):**
   - If no other format detected, assumes FASTA

---

## Implementation Details

### New Functions

1. **`readAlg(filename, format=None)`** - Main function (extended)
   - Auto-detects format if `format=None`
   - Routes to appropriate reader function

2. **`_detect_alignment_format(filename)`** - Format detection
   - Examines file header to determine format
   - Returns: `'fasta'`, `'stockholm'`, or `'clustal'`

3. **`_read_fasta(filename)`** - FASTA reader
   - Original FASTA reading logic (refactored)

4. **`_read_stockholm(filename)`** - Stockholm reader
   - Handles annotation lines
   - Combines multi-line sequences
   - Stops at `//` terminator

5. **`_read_clustal(filename)`** - Clustal reader
   - Skips conservation lines
   - Handles multi-block alignments
   - Combines sequence fragments

---

## Features

### ✅ Robust Parsing
- Handles multi-line sequences (all formats)
- Skips annotation/conservation lines appropriately
- Removes whitespace from sequence data
- Case-insensitive format detection

### ✅ Error Handling
- Validates sequence data exists before adding
- Handles edge cases (empty lines, malformed entries)
- Clear error messages for unsupported formats

### ✅ Backward Compatibility
- Existing code using `readAlg()` continues to work
- FASTA format remains default
- Same return format: `(headers, sequences)`

---

## Examples

### Example 1: Reading Pfam Stockholm Alignment

```python
from pysca import scaTools as sca

# Pfam alignments are typically in Stockholm format
headers, sequences = sca.readAlg('PF00071_full.stockholm')

print(f"Found {len(headers)} sequences")
print(f"Alignment length: {len(sequences[0])}")
```

### Example 2: Reading Clustal Alignment

```python
# Clustal format from alignment tools
headers, sequences = sca.readAlg('alignment.aln')  # Auto-detects Clustal

# Or explicitly specify
headers, sequences = sca.readAlg('alignment.aln', format='clustal')
```

### Example 3: Mixed Format Workflow

```python
# Read different formats in same workflow
fasta_headers, fasta_seqs = sca.readAlg('data1.fasta')
stockholm_headers, stockholm_seqs = sca.readAlg('data2.stockholm')
clustal_headers, clustal_seqs = sca.readAlg('data3.clustal')

# All return same format: (headers, sequences)
# Can be used interchangeably in downstream analysis
```

---

## Testing Recommendations

1. **Test with real alignments** from each format
2. **Verify sequence order** is preserved
3. **Check multi-line sequences** are correctly combined
4. **Test edge cases**: empty files, malformed entries, mixed formats
5. **Compare with BioPython** SeqIO results for validation

---

## Notes

- All sequences are converted to **uppercase** (consistent with original behavior)
- Sequence data whitespace is **removed** automatically
- Format detection is **heuristic** - may fail on unusual files
- For ambiguous cases, use explicit `format` parameter
- Stockholm annotation lines are **skipped** (only sequence data extracted)
- Clustal conservation lines are **skipped** (only sequence data extracted)

---

## Future Enhancements (Optional)

- Support for additional formats (PHYLIP, MSF, etc.)
- Write functions for each format
- Format conversion utilities
- More robust format detection
- Support for sequence quality scores (Stockholm)


