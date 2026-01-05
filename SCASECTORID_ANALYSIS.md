# scaSectorID_py3.py Analysis and Integration Plan

## Workflow Overview

`scaSectorID_py3.py` performs sector identification after SCA core calculations:

1. **Load Database** - Requires `db['sequence']` and `db['sca']` from previous steps
2. **Eigendecomposition** - Compute eigenvectors and eigenvalues of Csca
3. **Choose kpos** - Number of significant eigenmodes (auto or manual)
4. **ICA Rotation** - Rotate top kpos eigenvectors using ICA
5. **IC Identification** - Identify position sets for each independent component
6. **Sector Extraction** - Extract sectors using `sca.t()`
7. **Store Results** - Save to `db['sector']`

---

## Dependencies

**Required from `db['sca']`:**
- `Csca` - SCA correlation matrix (required)
- `Lrand` - Randomization eigenvalues (optional, for automatic kpos selection)

**Required from `db['sequence']`:**
- `ats` - Alignment-to-structure mapping (for sector_ats output)

**Functions called:**
- `sca.eigenVect(Csca)` - Eigendecomposition
- `sca.rotICA(V, kmax=kpos)` - ICA rotation
- `sca.icList(Vica, kpos, Csca, p_cut=cutoff)` - IC identification
- `sca.t(Vica, ic_list)` - Sector extraction

---

## Integration Strategy

### Option 1: Make it Optional in scaCore_py3.py (Recommended)

**Pros:**
- Users can still run steps separately if needed
- Backward compatible
- Can skip sector ID if not needed
- Clear separation of concerns

**Implementation:**
- Add `--do-sector-id` flag (default: False)
- Add sector ID options: `--kpos`, `--sector-cutoff`, `--kmax-cap`
- Run sector ID calculations after SCA core calculations
- Store results in `db['sector']`

### Option 2: Always Run (Not Recommended)

**Pros:**
- Simpler for users (one command)
- No intermediate files

**Cons:**
- Always computes sectors even if not needed
- Less flexible
- Harder to debug individual steps

---

## Proposed Integration

Add to `scaCore_py3.py`:

1. **New command-line arguments:**
   ```python
   --do-sector-id          # Enable sector identification
   --kpos INT             # Number of significant eigenmodes (0=auto)
   --sector-cutoff FLOAT  # Cutoff for IC selection (default 0.95)
   --kmax-cap INT         # Safety cap on kpos (default 50)
   ```

2. **Integration point:**
   - After SCA core calculations (after randomization)
   - Before database write
   - Only if `--do-sector-id` is set

3. **Code structure:**
   - Reuse functions from scaSectorID_py3.py (eigenVect, rotICA, icList, t)
   - Import helper functions (_as_index_lists, _mat_sanitize) if needed
   - Or inline the logic for simplicity

---

## Key Considerations

1. **Backward Compatibility:**
   - Default behavior unchanged (no sector ID)
   - Existing scripts continue to work
   - Can still run scaSectorID separately

2. **Error Handling:**
   - Check that Csca exists
   - Handle missing Lrand gracefully (fallback kpos selection)
   - Validate kpos, cutoff values

3. **Performance:**
   - Sector ID is relatively fast (seconds to minutes)
   - No major computational bottlenecks
   - Can add progress logging

4. **Output:**
   - Store in `db['sector']` same as scaSectorID
   - Same format for compatibility
   - MATLAB export should include sector if requested

---

## Implementation Plan

1. Add command-line arguments for sector ID
2. Add sector ID computation function or inline logic
3. Integrate after SCA core calculations
4. Store results in `db['sector']`
5. Update documentation
6. Test with existing databases


