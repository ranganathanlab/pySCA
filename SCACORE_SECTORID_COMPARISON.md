# scaSectorID Integration Status in scaCore

## Core Functionality: ✅ **Fully Integrated**

All essential sector identification functionality is present in `scaCore_py3.py` when `--do-sector-id` is enabled:

1. ✅ **Eigendecomposition** of Csca
2. ✅ **Automatic kpos selection** (using Lrand if available)
3. ✅ **ICA rotation** (rotICA)
4. ✅ **IC identification** (icList)
5. ✅ **Sector extraction** (sca.t)
6. ✅ **Sector index normalization** (`_sector_as_index_lists` helper function)
7. ✅ **ATS label mapping**
8. ✅ **Database storage** (all sector data fields)
9. ✅ **Progress logging**

## Minor Differences

### 1. MATLAB Export

**scaSectorID_py3.py:**
```python
mat_payload = {
    "sequence": db.get("sequence", {}),
    "sca": db.get("sca", {}),
    "sector": db.get("sector", {}),
}
savemat(mat_path, _mat_sanitize(mat_payload), ...)
```

**scaCore_py3.py (current):**
```python
savemat(mat_path, db, ...)  # Exports entire db
```

**Impact:** 
- scaSectorID only exports `sequence`, `sca`, `sector` keys
- scaCore exports entire database (may include other keys if present)
- scaSectorID uses `_mat_sanitize` to handle edge cases (None values, etc.)

**Recommendation:** 
- Current approach is fine for most cases
- If MATLAB compatibility issues arise, can add `_mat_sanitize` function

### 2. Database File Handling

**scaSectorID_py3.py:**
- Writes back to input database file (in-place update)

**scaCore_py3.py:**
- Writes to new output file in `Outputs/` directory

**Impact:**
- scaCore's approach is better for workflow management (non-destructive)
- Original database file is preserved

## Conclusion

✅ **scaSectorID functionality is fully contained in scaCore** for core sector identification.

The only differences are:
1. **MATLAB export format** (minor - exports full db vs selective keys)
2. **File handling** (scaCore uses separate output file - actually better)

**Recommendation:** 
- Use `scaCore_py3.py --do-sector-id` for new workflows
- Keep standalone `scaSectorID_py3.py` for:
  - Legacy compatibility
  - When you only want to run sector ID on existing SCA results
  - When you need the MATLAB-specific sanitization


