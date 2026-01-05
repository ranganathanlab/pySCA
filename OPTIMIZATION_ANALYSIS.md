# scaTools.py Optimization Analysis

## Executive Summary

This document provides a comprehensive analysis of `scaTools.py` (2,324 lines) identifying performance bottlenecks, memory inefficiencies, and code quality issues. The analysis focuses on functions that are computationally intensive and frequently called in the SCA workflow.

---

## Critical Performance Issues

### 1. **scaMat() - O(N_pos² × N_aa²) Double Loop** ⚠️ CRITICAL

**Location:** Lines 1341-1349

**Issue:** Nested loop computing SVD for every position pair is extremely slow.

```python
for i in range(N_pos):
    for j in range(i, N_pos):
        u, s, vt = np.linalg.svd(
            tildeC[N_aa * i : N_aa * (i + 1), N_aa * j : N_aa * (j + 1)]
        )
```

**Impact:** For typical alignments (200-500 positions), this performs 20,000-125,000 SVD operations.

**Optimization:**
- Use vectorized operations where possible
- Consider block-wise SVD computation
- Cache intermediate results
- Use symmetry: `Cspec[j,i] = Cspec[i,j]` (already done, but could be optimized)

### 2. **freq() - Sparse to Dense Conversion** ⚠️ HIGH

**Location:** Lines 1000-1001

**Issue:** Converts sparse matrix to dense unnecessarily:

```python
al2d = alg2binss(alg, Naa)  # Creates sparse matrix
freq1 = seqwn.dot(np.array(al2d.todense()))[0]  # Converts to dense
freq2 = np.array(
    al2d.T.dot(scipy.sparse.diags(seqwn[0], 0)).dot(al2d).todense()
)
```

**Impact:** For large alignments (1000+ sequences, 500+ positions), this creates 500K × 1000 dense matrices.

**Optimization:**
- Keep computations in sparse format
- Use `scipy.sparse` operations: `al2d.T.dot(scipy.sparse.diags(...).dot(al2d))`
- Only convert to dense when absolutely necessary

### 3. **alg2binss() and alg2bin() - Memory Inefficiency** ⚠️ HIGH

**Location:** Lines 731-766

**Issue:** Creates full dense tensor before reshaping:

```python
Abin_tensor = np.zeros((N_aa, N_pos, N_seq))  # Full dense array
for ia in range(N_aa):
    Abin_tensor[ia, :, :] = (alg == ia + 1).T  # Populates entire tensor
Abin = sparsify(Abin_tensor.reshape(...))
```

**Impact:** Allocates N_aa × N_pos × N_seq memory temporarily.

**Optimization:**
- Use `scipy.sparse.coo_matrix` or `csr_matrix` constructor directly
- Build sparse matrix incrementally
- Use NumPy advanced indexing instead of loops

### 4. **scaMat() - Dense Matrix Conversion** ⚠️ HIGH

**Location:** Line 1354

**Issue:** Converts sparse to dense unnecessarily:

```python
al2d = np.array(alg2binss(alg).todense())  # Converts entire sparse matrix
```

**Impact:** For large alignments, this can consume GB of memory.

**Optimization:**
- Keep sparse throughout projection computation
- Use sparse matrix slicing: `al2d[:, N_aa*i : N_aa*(i+1)]`

### 5. **filterSeq() - Inefficient String Operations** ⚠️ MEDIUM

**Location:** Lines 827-840

**Issue:** Uses `.count("-")` and list comprehension with repeated string comparisons:

```python
seqkeep0 = [
    s for s in range(Nseq) if alg0[s].count("-") / Npos < max_fracgaps
]
seqkeep = [
    s for s in seqkeep0
    if sum([alg0[s][i] == alg0[sref][i] for i in range(Npos)]) / Npos > min_seqid
]
```

**Impact:** O(Nseq × Npos) string operations.

**Optimization:**
- Convert to numpy array once, use vectorized operations
- Use `np.sum(alg_num == gap_code, axis=1) / Npos < max_fracgaps`

### 6. **lett2num() - Double Loop** ⚠️ MEDIUM

**Location:** Lines 724-728

**Issue:** Python-level loops instead of vectorized operations:

```python
for s, seq in enumerate(msa_lett):
    for i, lett in enumerate(seq):
        if lett in code:
            msa_num[s, i] = lett2index[lett]
```

**Impact:** O(Nseq × Npos) Python loops.

**Optimization:**
- Use NumPy's `np.searchsorted` or vectorized dictionary lookup
- Consider `np.fromiter` or `np.vectorize` (though may not be faster)

### 7. **numConnected() - Inefficient List Comprehensions** ⚠️ MEDIUM

**Location:** Lines 1539-1542

**Issue:** Nested list comprehensions with repeated max() operations:

```python
Vothers = [0 for i in range(Npos)]
for kk in range(kmax):
    if kk != k:
        Vothers = [max(Vothers[i], Vp[i, kk]) for i in range(Npos)]
```

**Impact:** O(kmax × Npos) with Python loops.

**Optimization:**
- Use NumPy: `Vothers = np.maximum.reduce([Vp[:, kk] for kk in range(kmax) if kk != k])`

### 8. **posWeights() - List Comprehension for Filtering** ⚠️ MEDIUM

**Location:** Line 1217

**Issue:** Creates list of indices instead of boolean mask:

```python
iok = [i for i in range(N_pos * N_aa) if (freq1[i] > 0 and freq1[i] < 1)]
Wia[iok] = ...
```

**Impact:** Memory overhead for large alignments.

**Optimization:**
- Use boolean indexing: `mask = (freq1 > 0) & (freq1 < 1)`
- Slightly faster and more memory efficient

### 9. **chooseRefSeq() - Redundant Computation** ⚠️ MEDIUM

**Location:** Lines 600-610

**Issue:** Computes all pairwise similarities, then extracts upper triangle:

```python
simMat = seqSim(numAlgNew)  # Computes full matrix
listS = [
    simMat[i, j]
    for i in range(simMat.shape[0])
    for j in range(i + 1, simMat.shape[1])
]
```

**Impact:** Unnecessary computation of symmetric matrix values.

**Optimization:**
- Use `np.triu_indices` or `simMat[np.triu_indices_from(simMat, k=1)]`
- More efficient and cleaner code

### 10. **randomize() - Missing Variable** ⚠️ BUG

**Location:** Line 1910

**Issue:** Uses undefined variable `Nseq`:

```python
if isinstance(seqw, int) and seqw == 1:
    seqw = np.ones((1, Nseq))  # Nseq is not defined here!
```

**Should be:**
```python
Nseq = msa_num.shape[0]
```

---

## Memory Optimization Opportunities

### 1. **Avoid Unnecessary Dense Conversions**

**Functions affected:** `freq()`, `scaMat()`, `projAlg()`, `seqSim()`

**Strategy:** Keep data in sparse format until final output.

### 2. **Use Views Instead of Copies**

**Location:** Multiple places using `copy.copy()` or array slicing

**Example - Line 1803:**
```python
Mtr = copy.copy(M)  # Creates full copy
```

**Optimization:**
- Use `np.copy()` only when necessary
- Use `np.triu()` and `np.tril()` views where possible

### 3. **Pre-allocate Arrays**

**Location:** Functions creating arrays in loops

**Example - Line 1858-1865 in randAlg():**
```python
col = np.array([], dtype=int)  # Grows array repeatedly
for aa, M in enumerate(Maa):
    col = np.append(col, np.tile(aa, M))  # Very slow!
```

**Optimization:**
```python
col = np.empty(Mseq, dtype=int)
idx = 0
for aa, M in enumerate(Maa):
    col[idx:idx+M] = aa
    idx += M
```

---

## Code Quality Issues

### 1. **Inconsistent Error Handling**

**Location:** Multiple functions catch `BaseException` (too broad)

**Example - Lines 198, 519, 563:**
```python
except BaseException as e:  # Too broad, catches KeyboardInterrupt, SystemExit
```

**Should be:**
```python
except Exception as e:  # or more specific exceptions
```

### 2. **File Handling Not Using Context Managers**

**Location:** Lines 104, 192, 2208, 2293, 2304, 2316

**Example - Line 104:**
```python
filelines = open(filename, "r").readlines()  # File not explicitly closed
```

**Optimization:**
```python
with open(filename, "r") as f:
    filelines = f.readlines()
```

### 3. **Inefficient String Operations**

**Location:** `readAlg()`, `clean_al()`, `filterSeq()`

**Example - Line 446:**
```python
seq_clean = ""  # String concatenation in loop
for aa in seq:
    if code.find(aa) != -1:  # .find() is inefficient
        seq_clean += aa  # Creates new string each time
```

**Optimization:**
```python
# Convert to list, join at end, or use translate/maketrans
valid_aa = set(code)
seq_list = [aa if aa in valid_aa else gap for aa in seq]
seq_clean = "".join(seq_list)
```

### 4. **Type Checking Issues**

**Location:** Multiple functions checking `isinstance(seqw, int) and seqw == 1`

**Issue:** This pattern is used to check for default value, but is fragile.

**Better approach:**
```python
if seqw is None or (isinstance(seqw, (int, np.integer)) and seqw == 1):
    seqw = np.ones((1, N_seq))
```

### 5. **Magic Numbers and Hardcoded Values**

**Location:** Throughout code (e.g., `N_aa=20`, gap fraction thresholds)

**Recommendation:** Move to constants at module level or function parameters.

### 6. **Code Duplication**

**Locations:**
- `AnnotPfam()` and `AnnotPfamDB()` have similar structure (lines 146-227, 229-307)
- Multiple functions creating sparse diagonal matrices similarly

---

## Numerical Stability Issues

### 1. **Division by Zero Risk**

**Location:** Line 1060 in `svdss()`

**Issue:**
```python
v = X.T.dot(u).dot(np.diag(1 / s))  # Division by zero if s contains zeros
```

**Optimization:**
```python
s_safe = np.where(s > 1e-12, s, 1.0)  # Prevent division by zero
v = X.T.dot(u).dot(np.diag(1 / s_safe))
```

### 2. **Log of Zero**

**Location:** Lines 1221-1229 in `posWeights()`

**Issue:** Computing `np.log()` on values that might be zero (though `iok` filters these).

**Already handled, but could be more explicit with `np.log1p()` or `np.finfo` checks.**

### 3. **Tolerance Checks**

**Location:** Multiple places use hardcoded tolerances (e.g., `1e-12`, `1e-15`)

**Recommendation:** Use `np.finfo(np.float64).eps` or make tolerances parameters.

---

## Specific Function Optimizations

### `scaMat()` - Highest Priority

**Current bottleneck:** Nested loop with SVD (lines 1341-1349)

**Proposed optimization:**
1. Use block-wise computation
2. Exploit symmetry more efficiently
3. Consider iterative methods for large matrices
4. Cache tildeC blocks if memory allows

### `freq()` - High Priority

**Current issue:** Sparse-to-dense conversions

**Proposed optimization:**
1. Keep `al2d` sparse
2. Use `scipy.sparse` operations for `freq2` computation
3. Only convert `freq1` result (which is small: N_aa*N_pos)

### `alg2binss()` - High Priority

**Current issue:** Creates dense tensor first

**Proposed optimization:**
```python
def alg2binss_optimized(alg, N_aa=20):
    """Optimized version building sparse matrix directly."""
    N_seq, N_pos = alg.shape
    rows = []
    cols = []
    data = []
    
    for ia in range(N_aa):
        matches = (alg == ia + 1)
        row_indices, pos_indices = np.where(matches)
        rows.extend(row_indices)
        cols.extend(ia * N_pos + pos_indices)
        data.extend([1] * len(row_indices))
    
    Abin = scipy.sparse.csr_matrix(
        (data, (rows, cols)), 
        shape=(N_seq, N_aa * N_pos)
    )
    return Abin
```

### `filterSeq()` - Medium Priority

**Current issue:** String operations in Python loops

**Proposed optimization:**
Convert to numeric first, use vectorized NumPy operations.

### `lett2num()` - Medium Priority

**Current issue:** Double Python loop

**Proposed optimization:**
Use NumPy vectorized operations or Cython/Numba if needed.

---

## Testing Recommendations

1. **Create benchmark suite** for critical functions
2. **Profile with `cProfile`** to identify actual bottlenecks
3. **Memory profiling** with `memory_profiler`
4. **Unit tests** to ensure optimizations don't change results
5. **Numerical accuracy tests** comparing old vs new implementations

---

## Implementation Priority

### Phase 1: Critical Fixes (Immediate)
1. Fix `randomize()` bug (undefined `Nseq`)
2. Fix file handling (use context managers)
3. Fix `randAlg()` array appending issue

### Phase 2: High Impact Optimizations (Week 1)
1. Optimize `alg2binss()` - sparse matrix construction
2. Optimize `freq()` - keep sparse
3. Optimize `scaMat()` - sparse conversions

### Phase 3: Medium Impact (Week 2)
1. Optimize `filterSeq()` - vectorized operations
2. Optimize `lett2num()` - vectorization
3. Optimize `numConnected()` - NumPy operations

### Phase 4: Code Quality (Week 3)
1. Error handling improvements
2. Code deduplication
3. Documentation updates
4. Type hints (optional)

---

## Expected Performance Gains

Based on typical alignments (500 sequences × 300 positions):

- **scaMat()**: 5-10x speedup (sparse operations + optimized loops)
- **freq()**: 3-5x speedup (sparse operations)
- **alg2binss()**: 2-3x speedup (direct sparse construction)
- **filterSeq()**: 10-20x speedup (vectorized operations)
- **Overall workflow**: 2-4x speedup for end-to-end analysis

**Memory savings:** 50-70% reduction for large alignments (1000+ sequences)

---

## Notes

- Some optimizations may require trade-offs (e.g., memory vs. speed)
- All optimizations should maintain numerical compatibility with existing results
- Consider using `numba` or `cython` for critical loops if further speedup needed
- Profile before and after to measure actual improvements


