# Preserving Legacy Files for Development

## The Situation

If you commit and push the deletions, the legacy notebooks will be **removed from the repository**. However, you have several options to preserve them for development:

## Option 1: Keep Legacy Files in a Separate Branch (RECOMMENDED)

Create a branch that keeps the legacy files before committing deletions:

```bash
# From your current pySCA directory
cd /Users/ramaranganathan/pySCA

# Create a branch with current state (includes legacy files)
git checkout -b legacy-files

# Commit current state (with legacy files)
git add -A
git commit -m "Snapshot with legacy notebooks and development docs"

# Go back to master
git checkout master

# Now commit deletions on master
git add -A
git commit -m "Release pySCA 7.0 - remove legacy files"
git push origin master

# Push the legacy-files branch (keeps legacy files available)
git push origin legacy-files
```

**Result:**
- `master` branch: Clean, no legacy files (for public release)
- `legacy-files` branch: Contains all legacy notebooks and docs (for development reference)

## Option 2: Use the Backup Branch We Created

We already created `dev-7.0-backup` earlier. Check if it has the legacy files:

```bash
# Check what's in the backup branch
git checkout dev-7.0-backup
ls notebooks/  # Should show legacy notebooks

# If it has them, you're good - they're preserved
# Go back to master
git checkout master
```

## Option 3: Keep Legacy Files Locally (Not in Git)

If you want to keep them for reference but not in the repository:

```bash
# Before committing deletions, copy legacy files to a backup location
mkdir -p ~/pySCA-legacy-backup
cp -r notebooks/SCA_*.ipynb ~/pySCA-legacy-backup/
cp notebooks/SCA_Example_Template.ipynb ~/pySCA-legacy-backup/
cp SETUP_NEW_REPO.md SYNC_BETWEEN_COMPUTERS.md ~/pySCA-legacy-backup/ 2>/dev/null || true
cp notebooks/CURSOR_SETUP.md notebooks/QUICK_FIX_CURSOR.md ~/pySCA-legacy-backup/ 2>/dev/null || true

# Now commit deletions
git add -A
git commit -m "Release pySCA 7.0"
```

## Option 4: Restore from Git History Later

Even after deletion, you can always get them back from git history:

```bash
# Find the commit before deletions
git log --oneline --all

# Checkout files from that commit
git checkout <commit-hash> -- notebooks/SCA_DHFR.ipynb
git checkout <commit-hash> -- notebooks/SCA_G.ipynb
# etc.
```

## Recommended Approach

**Use Option 1 (Separate Branch):**

1. **Create `legacy-files` branch** with current state (includes legacy files)
2. **Commit deletions on `master`** (clean for public release)
3. **Keep `legacy-files` branch** for development reference
4. **Make `master` public** (clean, no legacy files)
5. **Keep `legacy-files` branch private** (or don't push it if you want it local only)

## Quick Commands

```bash
# Preserve legacy files in a branch
git checkout -b legacy-files
git add -A
git commit -m "Snapshot: legacy notebooks and development docs preserved"
git push origin legacy-files  # Optional: push to remote

# Go back to master and commit deletions
git checkout master
git add -A
git commit -m "Release pySCA 7.0 - remove legacy files"
git push origin master

# Now:
# - master: Clean, ready for public release
# - legacy-files: Has all legacy files for reference
```

## Summary

**Question: Will legacy files be in pySCA-dev after committing deletions?**

- **On `master` branch**: ❌ No, they'll be removed
- **On `legacy-files` branch**: ✅ Yes, they'll be preserved
- **In git history**: ✅ Yes, always accessible via `git checkout <commit> -- <file>`

**Recommendation**: Create a `legacy-files` branch before committing deletions. This way you have:
- Clean `master` for public release
- `legacy-files` branch with everything for development reference
