# Committing and Pushing pySCA 7.0 Release

## Current Status

✅ Files deleted locally (need to be committed):
- Legacy notebooks (6 files)
- Development docs (4 files)
- Untitled.ipynb

✅ Files modified (need to be committed):
- Version updates
- Code changes (ic_cutoff, ic_pos, ic_ats)
- Documentation updates

## Steps to Commit and Push

### 1. Stage All Changes (Including Deletions)

```bash
# From your pySCA directory
cd /Users/ramaranganathan/pySCA

# Stage everything (including deletions)
git add -A

# Verify what will be committed
git status
```

### 2. Commit All Changes

```bash
git commit -m "Release pySCA 7.0

- Update version to 7.0 throughout codebase
- Rename sector-cutoff to ic-cutoff
- Rename sector_pos/sector_ats to ic_pos/ic_ats
- Remove user-specific paths from settings.py
- Remove legacy notebooks (SCA_DHFR, SCA_G, SCA_S1A, SCA_betalactamase, SCA_PDZ_Example, SCA_Example_Template)
- Remove development documentation (SETUP_NEW_REPO, SYNC_BETWEEN_COMPUTERS, CURSOR_SETUP, QUICK_FIX_CURSOR)
- Remove Untitled.ipynb
- Update documentation with Quick Start guide
- Update terminology (site-specific conservation, IC terminology)"
```

### 3. Create Version Tag

```bash
git tag -a v7.0 -m "pySCA version 7.0 - Public release

Major updates:
- Version 7.0 with updated terminology
- Improved documentation and Quick Start guide
- Code cleanup and modernization
- Removed legacy notebooks and development docs"
```

### 4. Push Everything

```bash
# Push commits
git push origin master

# Push the tag
git push origin v7.0
```

### 5. Verify on GitHub

1. Go to `https://github.com/ranganathanlab/pySCA-dev`
2. Check that:
   - ✅ Legacy notebooks are gone
   - ✅ Development docs are gone
   - ✅ Version is 7.0
   - ✅ Tag v7.0 exists

## What Gets Removed

After committing and pushing, these files will be **permanently removed** from the repository:

### Notebooks (6 files):
- `notebooks/SCA_DHFR.ipynb`
- `notebooks/SCA_G.ipynb`
- `notebooks/SCA_S1A.ipynb`
- `notebooks/SCA_betalactamase.ipynb`
- `notebooks/SCA_PDZ_Example.ipynb`
- `notebooks/SCA_Example_Template.ipynb`

### Development Docs (4 files):
- `SETUP_NEW_REPO.md`
- `SYNC_BETWEEN_COMPUTERS.md`
- `notebooks/CURSOR_SETUP.md`
- `notebooks/QUICK_FIX_CURSOR.md`

### Other:
- `Untitled.ipynb`

## Release Strategy Docs

The following files we created are now in `.gitignore` (won't be committed):
- `RELEASE_*.md`
- `PRE_RELEASE*.md`
- `FINAL_RELEASE*.md`
- `UPDATE_*.md`
- `USER_DISCOVERY*.md`

These are development artifacts and won't be part of the public release.

## Quick Command Sequence

```bash
cd /Users/ramaranganathan/pySCA
git add -A
git commit -m "Release pySCA 7.0"
git tag -a v7.0 -m "pySCA version 7.0"
git push origin master
git push origin v7.0
```

After this, the repository will be clean and ready for public release!
