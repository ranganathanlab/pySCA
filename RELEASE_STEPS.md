# Step-by-Step Release Process for pySCA 7.0

## Prerequisites
- All changes are ready (we've verified this)
- You have access to the GitHub/GitLab repository
- You have a backup branch (we created `dev-7.0-backup`)

## Step 1: Final Review and Commit

```bash
# Review all changes
git status
git diff --stat

# Stage all changes (including deletions)
git add -A

# Commit with descriptive message
git commit -m "Release pySCA 7.0

- Update version to 7.0 throughout codebase
- Rename sector-cutoff to ic-cutoff
- Rename sector_pos/sector_ats to ic_pos/ic_ats  
- Remove user-specific paths from settings.py
- Remove legacy notebooks and development docs
- Update documentation with Quick Start guide
- Update terminology (site-specific conservation, IC terminology)"

# Verify commit
git log -1
```

## Step 2: Create Version Tag

```bash
# Create annotated tag for version 7.0
git tag -a v7.0 -m "pySCA version 7.0 - Public release

Major updates:
- Version 7.0 with updated terminology
- Improved documentation and Quick Start guide
- Code cleanup and modernization"

# Verify tag
git tag -l
git show v7.0
```

## Step 3: Push to Public Repository

```bash
# Push commits to master/main branch
git push origin master

# Push the version tag
git push origin v7.0

# Verify on GitHub/GitLab that everything is pushed
```

## Step 4: Make Old Repository Private (if separate repo)

If you have a separate old pySCA repository that should be made private:

### On GitHub:
1. Go to the old repository (e.g., `ranganathanlab/pySCA-old` or similar)
2. Click **Settings** (top right)
3. Scroll down to **Danger Zone**
4. Click **Change visibility**
5. Select **Make private**
6. Confirm by typing the repository name

### On GitLab:
1. Go to the old repository
2. Click **Settings** → **General**
3. Expand **Visibility, project features, permissions**
4. Change **Project visibility** to **Private**
5. Click **Save changes**

## Step 5: Set Up Development Repository

You have several options:

### Option A: Use the backup branch for development (Recommended)

```bash
# Switch to development branch
git checkout dev-7.0-backup

# Create a new branch for ongoing development
git checkout -b dev

# Push development branch (make it private if needed)
git push origin dev

# Continue working on dev branch
# Master stays stable for releases
```

### Option B: Create a separate development repository

```bash
# On GitHub/GitLab, create a new repository:
# - Name: pySCA-dev or pySCA-development
# - Visibility: Private (for development)
# - Don't initialize with README

# Add it as a remote
git remote add dev https://github.com/ranganathanlab/pySCA-dev.git

# Push current state to dev repo
git push dev master

# Set up development branch
git checkout -b dev
git push dev dev
```

### Option C: Use the same repo with protected master branch

```bash
# On GitHub/GitLab:
# 1. Go to repository Settings → Branches
# 2. Add branch protection rule for master/main:
#    - Require pull request reviews
#    - Require status checks
#    - Include administrators
# 3. All development happens on dev/feature branches
# 4. Merge to master only via pull requests
```

## Step 6: Update Public Repository Description

On GitHub/GitLab, update the repository description:
- **Title**: pySCA - Statistical Coupling Analysis for Proteins
- **Description**: Python 3 implementation of Statistical Coupling Analysis (SCA) for studying evolutionary conservation in proteins. Version 7.0
- **Topics**: Add relevant tags (bioinformatics, proteins, evolution, SCA, etc.)

## Step 7: Create Release on GitHub/GitLab

### On GitHub:
1. Go to repository → **Releases** → **Create a new release**
2. **Tag**: Select `v7.0`
3. **Release title**: pySCA 7.0
4. **Description**: 
   ```
   ## pySCA 7.0 Release
   
   Major updates:
   - Updated terminology (IC-based naming)
   - Improved documentation with Quick Start guide
   - Code cleanup and modernization
   - Removed legacy notebooks (new ones coming soon)
   
   See [USAGE_INSTRUCTIONS.md](USAGE_INSTRUCTIONS.md) for details.
   ```
5. Check **Set as the latest release**
6. Click **Publish release**

### On GitLab:
1. Go to repository → **Releases** → **New release**
2. **Tag name**: `v7.0`
3. **Release title**: pySCA 7.0
4. **Release notes**: (same as above)
5. Click **Create release**

## Step 8: Verify Public Release

1. Visit the public repository URL
2. Verify:
   - [ ] Version 7.0 is visible in README
   - [ ] All files are present
   - [ ] No user-specific paths in settings.py
   - [ ] Documentation is complete
   - [ ] Release/tag v7.0 is visible
   - [ ] License file is present

## Step 9: Continue Development Workflow

### For ongoing development:

```bash
# Work on development branch
git checkout dev
git pull origin dev

# Make changes, commit
git add .
git commit -m "Description of changes"

# Push to dev branch
git push origin dev

# When ready for next release:
# 1. Merge dev into master (or create PR)
# 2. Update version number
# 3. Create new tag
# 4. Push to public
```

## Summary Checklist

- [ ] All changes committed
- [ ] Version tag v7.0 created
- [ ] Pushed to public repository
- [ ] Old repository made private (if applicable)
- [ ] Development branch/repo set up
- [ ] Release created on GitHub/GitLab
- [ ] Repository description updated
- [ ] Public release verified

## Quick Command Reference

```bash
# Complete release sequence:
git add -A
git commit -m "Release pySCA 7.0"
git tag -a v7.0 -m "pySCA version 7.0"
git push origin master
git push origin v7.0

# Set up development:
git checkout -b dev
git push origin dev
```
