# pySCA 7.0 Release Checklist

## Completed Tasks

### ✅ Version Updates
- [x] Updated `setup.py` version from 6.1 to 7.0
- [x] Updated `docs/source/conf.py` version and release from 6.1 to 7.0
- [x] Updated `scripts/runAllNBCalcs.sh` version from 6.1 to 7.0
- [x] Updated `README.md` to reflect version 7.0

### ✅ File Cleanup
- [x] Removed legacy notebooks:
  - `notebooks/SCA_betalactamase.ipynb`
  - `notebooks/SCA_DHFR.ipynb`
  - `notebooks/SCA_G.ipynb`
  - `notebooks/SCA_S1A.ipynb`
  - `notebooks/SCA_PDZ_Example.ipynb`
  - `notebooks/SCA_Example_Template.ipynb`
  - `Untitled.ipynb` (root directory)
- [x] Removed internal development documentation:
  - `notebooks/CURSOR_SETUP.md`
  - `notebooks/QUICK_FIX_CURSOR.md`
  - `SYNC_BETWEEN_COMPUTERS.md`
  - `SETUP_NEW_REPO.md`
  - `ANNOTATEMSA_OPTIMIZATION.md` (development planning doc)
- [x] Removed temporary log files:
  - `sc_PF00595.log`
  - `spm_PF00595.log`

### ✅ Documentation Updates
- [x] Updated `README.md`:
  - Added version 7.0 header
  - Added Quick Start section
  - Updated description to be more informative
  - Removed references to old notebooks
  - Added links to installation and usage instructions
- [x] Updated `notebooks/README.md`:
  - Removed references to deleted notebooks
  - Updated to indicate new notebooks coming soon
  - Kept general notebook setup instructions
- [x] Reviewed `INSTALLATION.md` - No changes needed (already accurate)
- [x] Reviewed `USAGE_INSTRUCTIONS.md` - No changes needed (already accurate)

### ✅ Git Backup
- [x] Created backup branch `dev-7.0-backup` for current development version

## Files Not Included in Release (Already in .gitignore)
- `robustica/` - Local clone for evaluation, excluded from release
- `Outputs/` - User output directory
- `Inputs/` - User input directory
- `*.log` files - Log files

## Remaining Development Documentation

The following analysis/optimization documentation files remain in the repository. These may be useful for developers but are not essential for end users:

- `OPTIMIZATION_ANALYSIS.md`
- `OPTIMIZATIONS_APPLIED_V2.md`
- `OPTIMIZATIONS_APPLIED.md`
- `OPTIMIZATIONS.md`
- `SCACORE_ANALYSIS.md`
- `SCACORE_OPTIMIZATIONS.md`
- `SCACORE_SECTORID_COMPARISON.md`
- `SCACORE_SECTORID_INTEGRATION.md`
- `SCAPROCESSMSA_ANALYSIS.md`
- `SCAPROCESSMSA_OPTIMIZATIONS.md`
- `SCASECTORID_ANALYSIS.md`
- `SEQSIM_OPTIMIZATION_SUMMARY.md`
- `SEQSIM_SUBSAMPLING.md`
- `MMSEQS2_INTEGRATION_ANALYSIS.md`
- `LARGE_MSA_OPTIMIZATIONS.md`
- `MODERNIZATION_GUIDE.md`
- `FORMAT_SUPPORT.md`
- `JUPYTER_NOTEBOOKS.md`
- `SCATOOLS_FUNCTIONS.md`
- `SCAVISUALIZESIMAT_USAGE.md`
- `SCAMAT_PSEUDOCODE.md`
- `UX_DESIGN.md`
- `UX_IMPLEMENTATION_GUIDE.md`
- `UX_QUICK_REFERENCE.md`
- `VERSIONS.rst`

**Recommendation:** Consider moving these to a `docs/development/` subdirectory or keeping them as-is for developer reference.

## Next Steps for Final Review

1. **Review all changes:**
   ```bash
   git status
   git diff
   ```

2. **Test installation:**
   ```bash
   pip install -e .
   ```

3. **Verify scripts work:**
   ```bash
   sca-process-msa --help
   sca-core --help
   ```

4. **Check documentation links:**
   - Verify all links in README.md work
   - Check that INSTALLATION.md and USAGE_INSTRUCTIONS.md are accessible

5. **Create new example notebooks** (as mentioned, coming soon)

6. **Final commit and tag:**
   ```bash
   git add -A
   git commit -m "Prepare pySCA 7.0 for release"
   git tag -a v7.0 -m "pySCA version 7.0"
   ```

## Summary

pySCA 7.0 is now prepared for release with:
- ✅ Version updated to 7.0 throughout
- ✅ Legacy notebooks removed
- ✅ Internal development files cleaned up
- ✅ Documentation updated and consistent
- ✅ README provides clear quick start
- ✅ Installation and usage instructions verified

The repository is ready for your final review and completeness check.
