# Pre-Release Checklist for pySCA 7.0

## ✅ Completed Checks

### 1. Version Consistency
- ✅ `setup.py`: version="7.0"
- ✅ `docs/source/conf.py`: version = '7.0', release = '7.0'
- ✅ `scripts/runAllNBCalcs.sh`: version=7.0
- ✅ `README.md`: Version 7.0 header

### 2. File Cleanup
- ✅ Legacy notebooks removed
- ✅ Internal development docs removed (CURSOR_SETUP, QUICK_FIX_CURSOR, SYNC_BETWEEN_COMPUTERS, SETUP_NEW_REPO)
- ✅ Temporary log files removed
- ✅ Untitled.ipynb removed

### 3. Documentation
- ✅ README.md updated with version 7.0 and Quick Start
- ✅ INSTALLATION.md reviewed (accurate)
- ✅ USAGE_INSTRUCTIONS.md updated with new terminology (ic_cutoff, ic_pos, ic_ats)
- ✅ notebooks/README.md updated

### 4. Code Changes
- ✅ `--sector-cutoff` → `--ic-cutoff` (scaCore, scaSectorID, notebook_utils)
- ✅ `sector_pos` → `ic_pos` (all files)
- ✅ `sector_ats` → `ic_ats` (all files)
- ✅ Database keys updated consistently

### 5. .gitignore
- ✅ Outputs/ and Inputs/ ignored
- ✅ robustica/ ignored
- ✅ *.log files ignored
- ✅ Standard Python ignores in place

## ⚠️ Issues Found

### 1. **CRITICAL: settings.py contains user-specific paths**
   - **File**: `pysca/settings.py`
   - **Issue**: Lines 17-19 contain `/Users/rama/Dropbox/transfer/dbs/PFAM/...`
   - **Action Required**: Replace with placeholder paths or default values
   - **Status**: ⚠️ NEEDS FIXING

### 2. **RELEASE_7.0_CHECKLIST.md is untracked**
   - **File**: `RELEASE_7.0_CHECKLIST.md`
   - **Issue**: Development checklist file
   - **Action Required**: Either add to .gitignore or remove before release
   - **Status**: ⚠️ NEEDS DECISION

### 3. **Temporary files in code**
   - **File**: `pysca/scaTools.py`
   - **Issue**: Uses `tmp_pdb_seq.fasta` and `tmp_algn_seq.fasta` (but these are cleaned up)
   - **Status**: ✅ OK (files are cleaned up after use)

## 📋 Pre-Release Actions Required

### Before Making Public:

1. **Fix settings.py** - Replace user-specific paths:
   ```python
   path2pfamseq = "/path/to/pfamseq.txt"  # User must set this
   path2pfamseqdb = "/path/to/pfamseq.db"  # User must set this
   ```

2. **Handle RELEASE_7.0_CHECKLIST.md**:
   - Option A: Delete it (it's a development artifact)
   - Option B: Add to .gitignore if you want to keep it locally

3. **Verify all changes are committed**:
   ```bash
   git add -A
   git status  # Review all changes
   git commit -m "Prepare pySCA 7.0 for public release"
   ```

4. **Create version tag**:
   ```bash
   git tag -a v7.0 -m "pySCA version 7.0"
   ```

5. **Final verification**:
   - [ ] No user-specific paths in code
   - [ ] No sensitive information (passwords, API keys, etc.)
   - [ ] All version numbers consistent
   - [ ] Documentation is complete and accurate
   - [ ] .gitignore properly configured
   - [ ] License file present (✅ LICENSE exists)

## 📝 Files to Review Before Release

- [ ] `pysca/settings.py` - **MUST FIX user paths**
- [ ] `RELEASE_7.0_CHECKLIST.md` - Decide: keep or remove
- [ ] All markdown files in root (some are development docs - OK to keep)
- [ ] `docs/source/` - Verify documentation is current

## 🎯 Ready for Release?

**Status**: ⚠️ **NOT YET** - settings.py needs to be fixed

After fixing settings.py, the repository will be ready for public release.
