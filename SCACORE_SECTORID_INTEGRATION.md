# scaCore_py3.py + scaSectorID_py3.py Integration

## Summary

Successfully integrated sector identification functionality from `scaSectorID_py3.py` into `scaCore_py3.py`, allowing users to run both SCA core calculations and sector identification in a single command.

---

## Integration Details

### New Command-Line Options

Added the following options to `scaCore_py3.py`:

```bash
--do-sector-id           # Enable sector identification after SCA core
--kpos INT              # Number of significant eigenmodes (0=auto)
--sector-cutoff FLOAT   # Cutoff for IC selection (default 0.95)
--kmax-cap INT          # Safety cap on kpos (default 50)
```

### Implementation

1. **Added helper function `_sector_as_index_lists()`**:
   - Converts sector output from `sca.t()` to standardized format
   - Handles multiple input formats (tuples, dicts, arrays, lists)
   - Same logic as `scaSectorID_py3.py` for consistency

2. **Integrated sector ID computation**:
   - Runs after SCA core calculations
   - Only executes if `--do-sector-id` flag is set
   - Uses existing Csca and Lrand from SCA core calculations
   - Stores results in `db['sector']` with same format as standalone scaSectorID

3. **Workflow:**
   ```
   Load DB → SCA Core → Sector ID (optional) → Save DB
   ```

---

## Usage Examples

### Run SCA Core Only (Original Behavior)

```bash
# Works exactly as before - backward compatible
python scaCore_py3.py input.db
```

### Run SCA Core + Sector ID in One Step

```bash
# Enable sector identification
python scaCore_py3.py input.db --do-sector-id

# With custom parameters
python scaCore_py3.py input.db \
    --do-sector-id \
    --kpos 6 \
    --sector-cutoff 0.95 \
    --kmax-cap 50

# Auto-select kpos based on Lrand
python scaCore_py3.py input.db --do-sector-id --kpos 0
```

### Complete Workflow (All in One)

```bash
# Process MSA
python scaProcessMSA_py3_big.py alignment.fasta --precluster

# Run SCA Core + Sector ID
python scaCore_py3.py output.db --do-sector-id

# Now you have db['sequence'], db['sca'], and db['sector']!
```

---

## Benefits

### 1. **Convenience**
- Single command instead of two
- No intermediate database writes/reads
- Faster overall execution (no I/O overhead between steps)

### 2. **Efficiency**
- Reuses computed Csca and Lrand directly
- No need to load database twice
- Reduced file I/O

### 3. **Flexibility**
- Optional: `--do-sector-id` flag (default: False)
- Can still run steps separately if needed
- Same parameters as standalone scaSectorID

### 4. **Backward Compatibility**
- Default behavior unchanged (no sector ID)
- Existing scripts continue to work
- Can still use standalone `scaSectorID_py3.py`

---

## Output Format

The integrated version produces the same `db['sector']` structure as standalone `scaSectorID_py3.py`:

```python
db['sector'] = {
    'kpos': int,                    # Number of eigenmodes used
    'cutoff': float,                # Cutoff parameter used
    'eigvals': np.ndarray,          # Top kpos eigenvalues
    'V': np.ndarray,                # Top kpos eigenvectors
    'Vica': np.ndarray,             # ICA-rotated eigenvectors
    'W': np.ndarray,                # ICA rotation matrix
    'ic_list': list,                # IC position sets
    'pd': list,                     # Probability distributions
    'scaled_pd': list,              # Scaled probability distributions
    'sector_pos': list[list[int]],  # Sector positions (0-based indices)
    'sector_ats': list[list[str]],  # Sector ATS labels
}
```

---

## Algorithm Flow

When `--do-sector-id` is enabled:

1. **Eigendecomposition** of Csca:
   ```python
   Vfull, eigvals = sca.eigenVect(Csca)
   ```

2. **Select kpos** (number of significant eigenmodes):
   - If `--kpos` > 0: use specified value
   - If `--kpos` = 0 (auto):
     - If Lrand available: compare eigvals to Lrand quantiles
     - Else: use default (min(6, num_eigvals))

3. **Extract top kpos eigenmodes**:
   ```python
   V = Vfull[:, :kpos]
   ```

4. **ICA rotation**:
   ```python
   Vica, W = sca.rotICA(V, kmax=kpos)
   ```

5. **Identify IC position sets**:
   ```python
   ic_list, pd, scaled_pd = sca.icList(Vica, kpos, Csca, p_cut=cutoff)
   ```

6. **Extract sectors**:
   ```python
   sectors_raw = sca.t(Vica, ic_list)
   sector_pos = _sector_as_index_lists(sectors_raw, Lpos)
   sector_ats = [[ats_list[i] for i in s] for s in sector_pos]
   ```

---

## Testing Recommendations

1. **Functional Testing:**
   - Compare output with standalone scaSectorID
   - Verify sector_pos and sector_ats match
   - Test with various kpos values (auto and manual)
   - Test with/without Lrand

2. **Backward Compatibility:**
   - Verify default behavior (no --do-sector-id) unchanged
   - Test existing scripts still work
   - Ensure standalone scaSectorID still works

3. **Edge Cases:**
   - Very small alignments
   - Large alignments
   - Missing Lrand (should fall back gracefully)
   - Invalid kpos values

---

## Migration Guide

### For Users Currently Using Both Scripts

**Before:**
```bash
python scaCore_py3.py input.db
python scaSectorID_py3.py input.db --kpos 6
```

**After (Option 1 - Integrated):**
```bash
python scaCore_py3.py input.db --do-sector-id --kpos 6
```

**After (Option 2 - Still Separate):**
```bash
# Still works! No changes needed
python scaCore_py3.py input.db
python scaSectorID_py3.py input.db --kpos 6
```

---

## Performance Impact

- **Time**: Sector ID is relatively fast (seconds to minutes), so integration adds minimal overhead
- **Memory**: No additional memory impact (reuses existing Csca)
- **I/O**: Eliminates one database read/write cycle

---

## Summary

✅ **Successfully integrated** sector ID into scaCore
✅ **Fully backward compatible** (optional flag, default off)
✅ **Same output format** as standalone scaSectorID
✅ **More convenient** for users (one command)
✅ **More efficient** (no intermediate I/O)

Users can now run the complete SCA analysis (core + sector ID) in a single command while maintaining the flexibility to run steps separately if needed.


