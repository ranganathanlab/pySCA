# Final Release and Development Setup for pySCA 7.0

## Current Status
✅ pySCA-dev is ready for public release
- All development files removed
- Documentation updated
- Version set to 7.0
- Code changes committed

## Step 1: Final Commit and Tag

```bash
cd /Users/ramaranganathan/pySCA

# Make sure everything is committed
git add -A
git status  # Review what will be committed

# Commit all remaining changes
git commit -m "Release pySCA 7.0

- Update version to 7.0 throughout
- Update README.md and USAGE_INSTRUCTIONS.md
- Remove development documentation files
- Update terminology (ic-cutoff, ic_pos, ic_ats)
- Remove user-specific paths from settings.py"

# Create version tag
git tag -a v7.0 -m "pySCA version 7.0 - Public release"

# Push everything
git push origin master
git push origin v7.0
```

## Step 2: Make pySCA-dev Public

### On GitHub:

1. Go to `https://github.com/ranganathanlab/pySCA-dev`
2. Click **Settings** (top right)
3. Scroll down to **Danger Zone**
4. Click **Change visibility**
5. Select **Make public**
6. Confirm by typing the repository name

### Optional: Rename to pySCA

If you want users to see "pySCA" instead of "pySCA-dev":

1. In the same Settings page
2. Scroll to **Repository name**
3. Change from `pySCA-dev` to `pySCA`
4. Click **Rename**

**Note:** If you rename, you'll need to update your local remote:
```bash
git remote set-url origin https://github.com/ranganathanlab/pySCA.git
```

## Step 3: Archive Old pySCA Repository

### On GitHub:

1. Go to `https://github.com/ranganathanlab/pySCA` (the old one)
2. Click **Settings** → **General**
3. Scroll to **Danger Zone**
4. Click **Archive this repository**
5. Confirm

**Before archiving, update the README** (as we discussed earlier):
- Add notice at top pointing to new v7.0 repository
- Then archive

## Step 4: Create GitHub Release

1. Go to the public repository (pySCA-dev or pySCA)
2. Click **Releases** → **Create a new release**
3. **Tag**: Select `v7.0`
4. **Release title**: `pySCA 7.0`
5. **Description**:
   ```
   # pySCA 7.0
   
   Python 3 implementation of Statistical Coupling Analysis (SCA) for studying evolutionary conservation in proteins.
   
   ## What's New in 7.0
   
   - Updated terminology: IC-based naming (ic_cutoff, ic_pos, ic_ats)
   - Improved documentation with Quick Start guide
   - Code cleanup and modernization
   - Better support for large MSAs with MMseqs2
   - Enhanced sector identification workflow
   
   ## Installation
   
   ```bash
   git clone https://github.com/ranganathanlab/pySCA.git
   cd pySCA
   pip install -e .
   ```
   
   See [INSTALLATION.md](INSTALLATION.md) for detailed setup instructions.
   ```
6. Check **Set as the latest release**
7. Click **Publish release**

## Step 5: Set Up Development Repository

You have two options:

### Option A: Keep pySCA-dev for Development (Recommended)

**If you renamed pySCA-dev to pySCA:**

1. **Create new private repository** on GitHub:
   - Name: `pySCA-dev`
   - Visibility: **Private**
   - Don't initialize with README

2. **Add as remote and push:**
   ```bash
   # Add new dev repo as remote
   git remote add dev https://github.com/ranganathanlab/pySCA-dev.git
   
   # Push current state to dev repo
   git push dev master
   
   # Create dev branch for ongoing work
   git checkout -b dev
   git push dev dev
   ```

3. **Update your workflow:**
   - Develop on `dev` branch in `pySCA-dev` (private)
   - Push releases to `pySCA` (public) when ready

### Option B: Use Same Repo with Protected Master

**If you kept pySCA-dev as the public repo:**

1. **Create development branch:**
   ```bash
   git checkout -b dev
   git push origin dev
   ```

2. **Protect master branch** (GitHub Settings → Branches):
   - Require pull request reviews
   - Require status checks
   - Include administrators

3. **Workflow:**
   - Develop on `dev` branch
   - Create pull requests to merge into `master`
   - Master stays stable for releases

## Step 6: Continue Development Workflow

### For Ongoing Development:

```bash
# Work on development branch/repo
git checkout dev
# or work in pySCA-dev repo

# Make changes
# ... edit files ...

# Commit
git add -A
git commit -m "Description of changes"

# Push to dev
git push origin dev
# or git push dev dev if using separate repo
```

### When Ready for Next Release:

```bash
# Merge dev into master (or create PR)
git checkout master
git merge dev

# Update version number
# Edit setup.py, docs/source/conf.py, etc.

# Commit version update
git add -A
git commit -m "Bump version to 7.1"

# Create new tag
git tag -a v7.1 -m "pySCA version 7.1"

# Push to public
git push origin master
git push origin v7.1

# Create GitHub release (as in Step 4)
```

## Summary Checklist

- [ ] Final commit and tag v7.0 created
- [ ] Pushed to remote (master and tag)
- [ ] Made repository public (or renamed to pySCA)
- [ ] Archived old pySCA repository (with updated README)
- [ ] Created GitHub release for v7.0
- [ ] Set up development repository/branch
- [ ] Verified public repository is accessible

## Quick Command Summary

```bash
# Final release
git add -A
git commit -m "Release pySCA 7.0"
git tag -a v7.0 -m "pySCA version 7.0"
git push origin master
git push origin v7.0

# Set up development (Option A - separate repo)
git remote add dev https://github.com/ranganathanlab/pySCA-dev.git
git checkout -b dev
git push dev dev

# Set up development (Option B - same repo)
git checkout -b dev
git push origin dev
```

After completing these steps, you'll have:
- ✅ Public pySCA repository with v7.0 release
- ✅ Archived old pySCA repository (with notice)
- ✅ Development environment set up for ongoing work
