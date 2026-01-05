# seqSim() Subsampling Implementation

## Overview

Enhanced `seqSim()` function with intelligent subsampling capabilities to handle very large MSAs without creating prohibitively large similarity matrices. Supports both MMseqs2 clustering-based and weighted random subsampling, with guaranteed retention of reference sequences.

---

## Problem Statement

For very large MSAs (e.g., 300,000 sequences), computing the full similarity matrix is:
- **Memory intensive**: 300k × 300k = 90 billion elements = **360GB** (impossible)
- **Computationally expensive**: O(N²) operations
- **Often unnecessary**: Full similarity pattern can be approximated with a diverse subsample

---

## Solution: Intelligent Subsampling

The enhanced `seqSim()` function provides two subsampling strategies:

1. **MMseqs2 Clustering** (recommended for large MSAs)
   - Uses sequence clustering to select diverse representatives
   - Better preserves sequence similarity patterns
   - Automatically falls back if MMseqs2 unavailable

2. **Weighted Random Selection** (default fallback)
   - Uses sequence weights for biased random sampling
   - Faster than MMseqs2 (no external tool needed)
   - Good for smaller subsamples or when weights are meaningful

**Key Feature**: Reference sequence (`i_ref`) is **always retained** regardless of subsampling method.

---

## Function Signature

```python
def seqSim(alg, max_seqs=None, i_ref=None, use_mmseqs2=False, 
           cluster_id=0.85, seqw=None, return_indices=False):
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alg` | np.ndarray | required | Numeric alignment (M × L) from `lett2num()` |
| `max_seqs` | int | None | Maximum sequences to include. If None, uses all sequences. |
| `i_ref` | int | None | Index of reference sequence to always retain |
| `use_mmseqs2` | bool | False | Use MMseqs2 clustering for subsampling |
| `cluster_id` | float | 0.85 | MMseqs2 identity threshold (if `use_mmseqs2=True`) |
| `seqw` | np.ndarray | None | Sequence weights for weighted random selection |
| `return_indices` | bool | False | Return selected sequence indices |

### Returns

- **Standard**: `simMat` (np.ndarray or scipy.sparse matrix)
- **With `return_indices=True`**: `(simMat, selected_indices)`

---

## Usage Examples

### 1. Standard Usage (All Sequences)

```python
# Works as before (backward compatible)
simMat = sca.seqSim(alg)
```

### 2. Simple Subsample (Weighted Random)

```python
# Subsample to 10,000 sequences, always keep reference sequence 0
simMat = sca.seqSim(alg, max_seqs=10000, i_ref=0)
```

### 3. MMseqs2-Based Subsample (Recommended)

```python
# Use MMseqs2 clustering for diverse subsampling
simMat, indices = sca.seqSim(
    alg, 
    max_seqs=10000, 
    i_ref=0, 
    use_mmseqs2=True,
    return_indices=True
)
print(f"Selected {len(indices)} sequences")
```

### 4. With Sequence Weights

```python
# Use existing sequence weights for weighted random selection
seqw = sca.seqWeights(alg)  # Compute weights on full alignment
simMat = sca.seqSim(alg, max_seqs=10000, i_ref=0, seqw=seqw)
```

---

## Implementation Details

### MMseqs2 Subsampling (`_subsample_with_mmseqs2`)

**Workflow:**
1. Convert numeric alignment to ungapped FASTA sequences
2. Run MMseqs2 clustering:
   - `mmseqs createdb` - Create database
   - `mmseqs cluster` - Cluster sequences at specified identity
   - `mmseqs createtsv` - Get cluster representatives
3. Parse cluster representatives (one per cluster)
4. Ensure reference sequence is included
5. If more representatives than `max_seqs`, randomly sample down
6. Return selected sequence indices

**Advantages:**
- **Diverse selection**: Representatives from different clusters
- **Preserves patterns**: Captures sequence similarity structure
- **Efficient**: MMseqs2 is highly optimized for clustering
- **Automatic fallback**: Falls back to weighted random if MMseqs2 fails

**Temporary Files:**
- Creates temporary directory for MMseqs2 operations
- Automatically cleans up on completion
- Uses `tempfile.mkdtemp()` for safe temporary directory creation

### Weighted Random Subsampling (`_subsample_weighted_random`)

**Workflow:**
1. Normalize sequence weights
2. Use `randSel()` for weighted random selection
3. Always include reference sequence in `keepSeq` list
4. Return selected sequence indices

**Advantages:**
- **Fast**: No external tools needed
- **Weighted**: Respects sequence importance
- **Simple**: Straightforward implementation

---

## Performance Comparison

### Without Subsampling (300k sequences)

| Metric | Value |
|--------|-------|
| Matrix size | 300k × 300k = 90B elements |
| Memory | ~360GB (dense) |
| Computation | O(N²) = 90 billion operations |
| **Status** | ❌ **Infeasible** |

### With MMseqs2 Subsampling (300k → 10k sequences)

| Metric | Value |
|--------|-------|
| Subsampling time | ~5-10 minutes (MMseqs2 clustering) |
| Matrix size | 10k × 10k = 100M elements |
| Memory | ~400MB (dense) |
| Computation | O(N'²) where N' = 10k |
| **Status** | ✅ **Feasible** |

**Speedup**: ~100-1000x faster (depending on hardware)
**Memory reduction**: ~900x less memory

---

## Backward Compatibility

✅ **Fully backward compatible**:
- All existing code using `seqSim(alg)` continues to work
- Default behavior (no subsampling) unchanged
- No breaking changes to function signature or return values

**Migration path for large MSAs:**
```python
# Old code (may fail for large MSAs)
simMat = sca.seqSim(alg)

# New code (handles large MSAs)
simMat = sca.seqSim(alg, max_seqs=10000, i_ref=0, use_mmseqs2=True)
```

---

## Integration with Existing Code

### scaCore Integration

The `scaCore` script can be updated to use subsampling:

```python
# In scaCore or similar scripts
if msa_num.shape[0] > 50000:
    # Use subsampling for very large alignments
    simMat = sca.seqSim(msa_num, max_seqs=50000, i_ref=i_ref, use_mmseqs2=True)
else:
    # Standard approach for smaller alignments
    simMat = sca.seqSim(msa_num)
```

### chooseRefSeq Integration

The `chooseRefSeq()` function already uses subsampling for >1000 sequences, but could benefit from MMseqs2:

```python
# Current approach (weighted random)
if len(alg) > 1000:
    seqw = seqWeights(alg)
    keep_seq = randSel(seqw, 1000)

# Could use MMseqs2 for better diversity
# (future enhancement)
```

---

## Recommendations

### When to Use Subsampling

**✅ Recommended:**
- MSAs with >50,000 sequences
- Limited memory resources
- Exploratory analysis where full matrix not needed
- When reference sequence must be preserved

**⚠️ Consider full matrix:**
- Small MSAs (<10,000 sequences)
- When exact full similarity matrix is required
- When downstream analysis needs all sequences

### Parameter Selection

**`max_seqs`:**
- **10,000-50,000**: Good balance for most analyses
- **5,000-10,000**: For memory-constrained systems
- **50,000+**: For detailed analysis with sufficient memory

**`cluster_id` (MMseqs2):**
- **0.80-0.85**: Balanced diversity (recommended)
- **0.90**: Less aggressive (more representatives)
- **0.75**: More aggressive (fewer representatives)

**`use_mmseqs2`:**
- **True**: Recommended for large MSAs (>100k sequences)
- **False**: Faster, simpler, good for smaller subsamples

---

## Error Handling

The implementation includes robust error handling:

1. **MMseqs2 not available**: Falls back to weighted random selection with warning
2. **MMseqs2 failure**: Catches exceptions and falls back gracefully
3. **Invalid parameters**: Validates inputs (e.g., i_ref within bounds)
4. **Empty alignments**: Handles edge cases safely

---

## Testing Recommendations

1. **Functional Testing**:
   - Verify reference sequence is always included
   - Check that subsampling reduces matrix size correctly
   - Validate similarity matrix structure (symmetry, diagonal)

2. **Performance Testing**:
   - Benchmark subsampling time vs full matrix computation
   - Measure memory usage with/without subsampling
   - Test with progressively larger MSAs

3. **Validation**:
   - Compare subsampled similarity patterns to full matrix (on smaller test cases)
   - Verify MMseqs2 fallback works correctly
   - Test edge cases (i_ref at boundaries, very small max_seqs)

---

## Future Enhancements

1. **Automatic max_seqs selection**: Based on available memory
2. **Adaptive clustering**: Adjust cluster_id to get desired number of representatives
3. **Hybrid approach**: Combine MMseqs2 clustering with additional random sampling for fine-tuning
4. **Parallel processing**: Speed up MMseqs2 operations for very large MSAs
5. **Caching**: Cache MMseqs2 clustering results for repeated analyses

---

## Summary

✅ **Intelligent subsampling** implemented with MMseqs2 support
✅ **Reference sequence** always retained
✅ **Backward compatible** with existing code
✅ **Robust error handling** with automatic fallback
✅ **Performance optimized** for very large MSAs
✅ **Flexible API** supporting multiple use cases

The enhanced `seqSim()` function now makes similarity matrix computation **feasible for very large MSAs** while preserving the essential sequence similarity patterns needed for analysis.


