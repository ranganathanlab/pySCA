# Cleanup Development Files for Public Release

## Files to Remove from Public Release

### Development/Release Documentation (Remove):
- `FINAL_RELEASE_SUMMARY.md`
- `PRE_RELEASE_CHECKLIST.md`
- `RELEASE_7.0_CHECKLIST.md`
- `RELEASE_STEPS.md`
- `RELEASE_STRATEGY.md`
- `UPDATE_OLD_REPO_README.md`
- `USER_DISCOVERY_STRATEGY.md`
- `SETUP_NEW_REPO.md` (already deleted locally)
- `SYNC_BETWEEN_COMPUTERS.md` (already deleted locally)

### Internal Analysis/Optimization Docs (Remove):
- `SCACORE_ANALYSIS.md`
- `SCACORE_OPTIMIZATIONS.md`
- `SCACORE_SECTORID_COMPARISON.md`
- `SCACORE_SECTORID_INTEGRATION.md`
- `SCAPROCESSMSA_ANALYSIS.md`
- `SCAPROCESSMSA_OPTIMIZATIONS.md`
- `SCASECTORID_ANALYSIS.md`
- `OPTIMIZATION_ANALYSIS.md`
- `OPTIMIZATIONS_APPLIED_V2.md`
- `OPTIMIZATIONS_APPLIED.md`
- `OPTIMIZATIONS.md`
- `SEQSIM_OPTIMIZATION_SUMMARY.md`
- `SEQSIM_SUBSAMPLING.md`
- `MMSEQS2_INTEGRATION_ANALYSIS.md`
- `LARGE_MSA_OPTIMIZATIONS.md`
- `MODERNIZATION_GUIDE.md`
- `SCAMAT_PSEUDOCODE.md`

### Internal UX/Design Docs (Remove):
- `UX_DESIGN.md`
- `UX_IMPLEMENTATION_GUIDE.md`
- `UX_QUICK_REFERENCE.md`

### Files to Keep (User-Facing):
- `README.md` ✅
- `INSTALLATION.md` ✅
- `USAGE_INSTRUCTIONS.md` ✅
- `LICENSE` ✅
- `setup.py` ✅
- `FORMAT_SUPPORT.md` ✅ (useful for users)
- `JUPYTER_NOTEBOOKS.md` ✅ (might be useful)
- `SCATOOLS_FUNCTIONS.md` ✅ (function reference - useful)
- `SCAVISUALIZESIMAT_USAGE.md` ✅ (usage guide - useful)
- `VERSIONS.rst` ✅ (version history - useful)

## Commands to Remove Development Files

```bash
cd /Users/ramaranganathan/pySCA

# Remove development/release docs
git rm FINAL_RELEASE_SUMMARY.md \
       PRE_RELEASE_CHECKLIST.md \
       RELEASE_7.0_CHECKLIST.md \
       RELEASE_STEPS.md \
       RELEASE_STRATEGY.md \
       UPDATE_OLD_REPO_README.md \
       USER_DISCOVERY_STRATEGY.md

# Remove analysis/optimization docs
git rm SCACORE_ANALYSIS.md \
       SCACORE_OPTIMIZATIONS.md \
       SCACORE_SECTORID_COMPARISON.md \
       SCACORE_SECTORID_INTEGRATION.md \
       SCAPROCESSMSA_ANALYSIS.md \
       SCAPROCESSMSA_OPTIMIZATIONS.md \
       SCASECTORID_ANALYSIS.md \
       OPTIMIZATION_ANALYSIS.md \
       OPTIMIZATIONS_APPLIED_V2.md \
       OPTIMIZATIONS_APPLIED.md \
       OPTIMIZATIONS.md \
       SEQSIM_OPTIMIZATION_SUMMARY.md \
       SEQSIM_SUBSAMPLING.md \
       MMSEQS2_INTEGRATION_ANALYSIS.md \
       LARGE_MSA_OPTIMIZATIONS.md \
       MODERNIZATION_GUIDE.md \
       SCAMAT_PSEUDOCODE.md

# Remove UX/design docs
git rm UX_DESIGN.md \
       UX_IMPLEMENTATION_GUIDE.md \
       UX_QUICK_REFERENCE.md

# Commit the removals
git commit -m "Remove development and analysis documentation from public release"
```

## Alternative: Keep in Development Branch

Before removing, create a branch that keeps these files:

```bash
# Create branch with all development docs
git checkout -b dev-docs

# Commit current state (with all docs)
git add -A
git commit -m "Development branch with all analysis and documentation"

# Go back to master
git checkout master

# Now remove the files
# ... (use git rm commands above) ...

# Commit
git commit -m "Remove development docs for public release"
```

## Summary

**Remove from public release:**
- 8 release/development strategy docs
- 17 analysis/optimization docs  
- 3 UX/design docs
- **Total: ~28 files**

**Keep for users:**
- Core documentation (README, INSTALLATION, USAGE_INSTRUCTIONS)
- Function reference (SCATOOLS_FUNCTIONS)
- Format support docs
- Version history
