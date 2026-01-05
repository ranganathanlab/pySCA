# seqSim() Subsampling Optimization Summary

## Overview

Implemented intelligent automatic subsampling strategy for `seqSim()` that optimizes memory usage while preserving important sequences.

---

## Key Features

### 1. ✅ Automatic Subsampling Based on Effective Sequences

**Default behavior (when `seqw` is provided and `max_seqs` is None):**
- Automatically calculates `max_seqs = 1.5 × sum(seqw)`
- `sum(seqw)` represents the effective number of sequences (M_eff)
- Ensures sufficient sampling while avoiding excessive memory usage

**Example:**
```python
# If sum(seqw) = 5000, then max_seqs = 1.5 × 5000 = 7500
simMat = seqSim(alg, seqw=seqw, i_ref=0)
```

### 2. ✅ Memory Cap Protection

**Parameter: `max_seqs_cap` (default: 50000)**
- Even if 1.5 × sum(seqw) exceeds the cap, limits to `max_seqs_cap`
- Prevents memory issues for very large MSAs
- Issues warning when capping is applied

**Example:**
```python
# If 1.5 × sum(seqw) = 100000, caps to 50000
simMat = seqSim(alg, seqw=seqw, max_seqs_cap=50000)
```

### 3. ✅ Reference Sequence Always Retained

**Parameter: `i_ref`**
- Reference sequence is always included in subsample
- Works with both automatic and manual subsampling

**Example:**
```python
simMat = seqSim(alg, seqw=seqw, i_ref=0)  # Sequence 0 always included
```

### 4. ✅ User-Specified Sequences Always Retained

**New Parameter: `keep_indices`**
- List of additional sequence indices to always retain
- Useful for sequences with special annotations to study
- Combined with `i_ref` - all specified sequences are guaranteed

**Example:**
```python
# Keep reference (0) + sequences 5, 10, 15 (e.g., annotated sequences)
simMat = seqSim(alg, seqw=seqw, i_ref=0, keep_indices=[5, 10, 15])
```

---

## Updated Function Signature

```python
def seqSim(alg, max_seqs=None, i_ref=None, use_mmseqs2=False, 
           cluster_id=0.85, seqw=None, return_indices=False,
           max_seqs_cap=50000, keep_indices=None):
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `alg` | np.ndarray | required | Numeric alignment (M × L) |
| `max_seqs` | int | None | Manual max sequences. If None and seqw provided, auto-calculates as 1.5 × sum(seqw) |
| `i_ref` | int | None | Reference sequence index (always retained) |
| `seqw` | np.ndarray | None | Sequence weights. Used for: (1) auto max_seqs, (2) weighted selection |
| `keep_indices` | list[int] | None | Additional sequence indices to always retain |
| `max_seqs_cap` | int | 50000 | Maximum sequences cap (prevents memory issues) |
| `use_mmseqs2` | bool | False | Use MMseqs2 clustering for subsampling |
| `cluster_id` | float | 0.85 | MMseqs2 identity threshold |
| `return_indices` | bool | False | Return selected indices |

---

## Usage Examples

### Automatic Subsampling (Recommended)

```python
# Automatically subsample to 1.5 × effective sequences
# Keeps reference sequence, respects memory cap
simMat = seqSim(alg, seqw=seqw, i_ref=0)

# With user-specified sequences to keep
simMat = seqSim(alg, seqw=seqw, i_ref=0, keep_indices=[5, 10, 15])

# Custom memory cap
simMat = seqSim(alg, seqw=seqw, i_ref=0, max_seqs_cap=100000)
```

### Manual Subsampling

```python
# Explicitly specify max_seqs (overrides automatic calculation)
simMat = seqSim(alg, max_seqs=10000, i_ref=0)

# Without seqw (no automatic calculation)
simMat = seqSim(alg)  # Uses all sequences
```

### With MMseqs2

```python
# Automatic subsampling with MMseqs2 clustering
simMat = seqSim(alg, seqw=seqw, i_ref=0, use_mmseqs2=True)

# Keep specific sequences with MMseqs2
simMat = seqSim(alg, seqw=seqw, i_ref=0, keep_indices=[5, 10], use_mmseqs2=True)
```

---

## Integration with scaCore_py3.py

### New Command-Line Options

```bash
--seqcorr-keep-indices INT [INT ...]
    Additional sequence indices to always retain in seqSim subsampling.
    Example: --seqcorr-keep-indices 5 10 15

--seqcorr-max-cap INT
    Maximum sequences cap for seqSim subsampling (default 50000).
    Example: --seqcorr-max-cap 100000
```

### Automatic Behavior in scaCore

When `--do-seqcorr` is enabled:
- Automatically uses `seqw` from database
- Automatically calculates `max_seqs = 1.5 × sum(seqw)` (capped)
- Automatically uses `i_ref` from database if available
- User can specify additional sequences with `--seqcorr-keep-indices`

**Example:**
```bash
python scaCore_py3.py input.db \
    --do-seqcorr \
    --seqcorr-keep-indices 5 10 15 \
    --seqcorr-max-cap 100000
```

---

## Algorithm Details

### Subsampling Logic

1. **Collect required sequences:**
   - Add `i_ref` if specified
   - Add all `keep_indices` if specified
   - Validate indices are within bounds

2. **Calculate max_seqs:**
   - If `max_seqs` is None and `seqw` provided:
     - `max_seqs = round(1.5 × sum(seqw))`
     - Apply cap: `max_seqs = min(max_seqs, max_seqs_cap)`
     - Ensure: `max_seqs >= len(must_keep)`
   - If `max_seqs` is None and `seqw` is None:
     - No subsampling (use all sequences)

3. **Perform subsampling:**
   - MMseqs2 clustering (if `use_mmseqs2=True`)
   - OR weighted random selection (if `use_mmseqs2=False`)
   - Always includes all `must_keep` sequences

4. **Verify and report:**
   - Verify all required sequences are included
   - Report subsampling statistics

---

## Benefits

### 1. **Intelligent Defaults**
- Automatically adapts to alignment characteristics (via effective sequences)
- No need to manually guess appropriate `max_seqs`

### 2. **Memory Safety**
- Memory cap prevents excessive memory usage
- Automatic calculation balances sampling vs. memory

### 3. **Preservation of Important Sequences**
- Reference sequence always retained
- User-specified sequences always retained
- Critical for analysis continuity

### 4. **Flexibility**
- Can override automatic calculation with manual `max_seqs`
- Can specify additional sequences to keep
- Can adjust memory cap as needed

---

## Performance Impact

### Example: 300k sequences, M_eff = 5000

**Before (no subsampling):**
- Matrix size: 300k × 300k = 360GB ❌ **Infeasible**

**After (automatic subsampling):**
- `max_seqs = 1.5 × 5000 = 7500` (capped at 50000)
- Matrix size: 7500 × 7500 = 450MB ✅ **Feasible**
- **Memory reduction: 800x**

### Example: 100k sequences, M_eff = 8000

**Automatic calculation:**
- `max_seqs = 1.5 × 8000 = 12000` (within cap)
- Matrix size: 12k × 12k = 1.15GB ✅ **Feasible**

---

## Backward Compatibility

✅ **Fully backward compatible:**
- Default behavior when `max_seqs=None` and `seqw=None`: uses all sequences (unchanged)
- All existing code continues to work
- New parameters are optional with sensible defaults

**Migration:**
```python
# Old code (still works)
simMat = seqSim(alg)

# New optimized code (automatic subsampling)
simMat = seqSim(alg, seqw=seqw, i_ref=0)
```

---

## Disabling Automatic Subsampling

**New Parameter: `auto_subsample` (default: True)**
- Set to `False` to disable automatic 1.5 × M_eff calculation
- When `False` and `max_seqs=None`: uses all sequences (original behavior)
- When `False` but `max_seqs` is specified: uses the specified value

**Example:**
```python
# Disable automatic subsampling, use all sequences
simMat = seqSim(alg, seqw=seqw, i_ref=0, auto_subsample=False)

# Disable automatic but use manual limit
simMat = seqSim(alg, seqw=seqw, i_ref=0, max_seqs=10000, auto_subsample=False)
```

**In scaCore_py3.py:**
```bash
# Disable automatic subsampling
python scaCore_py3.py input.db --do-seqcorr --no-auto-seqcorr-subsample

# Or use manual limit
python scaCore_py3.py input.db --do-seqcorr --max-seqcorr-seqs 10000 --no-auto-seqcorr-subsample
```

---

## Summary

✅ **Automatic subsampling** based on effective sequences (1.5 × sum(seqw)) - **can be disabled**
✅ **Memory cap** protection (default 50000, configurable)
✅ **Reference sequence** always retained
✅ **User-specified sequences** always retained (new `keep_indices` parameter)
✅ **Option to disable** automatic subsampling (new `auto_subsample` parameter)
✅ **Fully backward compatible** with existing code
✅ **Intelligent defaults** reduce need for manual parameter tuning

The optimized `seqSim()` now provides intelligent automatic subsampling while ensuring important sequences are always preserved, making it both more efficient and more reliable for large MSAs. Users can disable automatic subsampling if they prefer explicit control or want to use all sequences.

