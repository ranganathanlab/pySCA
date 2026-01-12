# Release Strategy for pySCA 7.0

## Repository Structure

### Current State:
- **pySCA-dev**: Your development repository (currently: `https://github.com/ranganathanlab/pySCA-dev.git`)
- **pySCA** (old): The previous public repository (needs to be made private)

### Recommended Approach:

## Option 1: Rename pySCA-dev to pySCA (Recommended)

**What happens:**
1. `pySCA-dev` → renamed to `pySCA` and made **public**
2. Old `pySCA` → renamed to `pySCA-legacy` or `pySCA-v6` and made **private**
3. Create new `pySCA-dev` repository (private) for ongoing development

**What users see:**
- Public repository: `https://github.com/ranganathanlab/pySCA`
- Release title: **"pySCA 7.0"** (or just "7.0")
- Users download: Source code zip/tar.gz from the release, or clone `git clone https://github.com/ranganathanlab/pySCA.git`
- Install command: `pip install git+https://github.com/ranganathanlab/pySCA.git` or `pip install -e .` after cloning

**Benefits:**
- Clean: Users see "pySCA" not "pySCA-dev"
- Standard: Follows typical open-source naming
- Simple: One public repo for releases

## Option 2: Keep pySCA-dev as Development, Push to Separate pySCA

**What happens:**
1. Keep `pySCA-dev` as **private** development repository
2. Create/use public `pySCA` repository for releases
3. Push releases from dev to public

**What users see:**
- Public repository: `https://github.com/ranganathanlab/pySCA`
- Release title: **"pySCA 7.0"**
- Users download: From public `pySCA` repo (not pySCA-dev)

**Benefits:**
- Separation: Dev work stays private until ready
- Control: You decide when to push to public

## What Users Will Download

### From GitHub Release:
1. **Source code (zip)**: `pySCA-7.0.zip` or `Source code (zip)` button
2. **Source code (tar.gz)**: `pySCA-7.0.tar.gz` or `Source code (tar.gz)` button
3. **Clone repository**: `git clone https://github.com/ranganathanlab/pySCA.git`

### Installation Methods:
```bash
# Method 1: Clone and install
git clone https://github.com/ranganathanlab/pySCA.git
cd pySCA
pip install -e .

# Method 2: Direct pip install (if you set up PyPI later)
pip install pySCA

# Method 3: Download release zip and install
# Download v7.0 release zip, extract, then:
cd pySCA-7.0
pip install -e .
```

## Release Title and Description

**Release Title:** `pySCA 7.0` (or just `7.0`)

**Release Description Example:**
```
# pySCA 7.0

Python 3 implementation of Statistical Coupling Analysis (SCA) for studying evolutionary conservation in proteins.

## What's New in 7.0

- Updated terminology: IC-based naming (ic_cutoff, ic_pos, ic_ats)
- Improved documentation with Quick Start guide
- Code cleanup and modernization
- Removed legacy notebooks (new examples coming soon)

## Installation

```bash
git clone https://github.com/ranganathanlab/pySCA.git
cd pySCA
pip install -e .
```

See [INSTALLATION.md](INSTALLATION.md) for detailed setup instructions.

## Documentation

- [Quick Start Guide](USAGE_INSTRUCTIONS.md#quick-start)
- [Installation Instructions](INSTALLATION.md)
- [Usage Guide](USAGE_INSTRUCTIONS.md)
- [Website](https://ranganathanlab.gitlab.io/pySCA)
```

## Recommended Workflow

**I recommend Option 1** (rename pySCA-dev to pySCA):

1. **Make pySCA-dev public** and rename to `pySCA`
2. **Make old pySCA private** and rename to `pySCA-legacy`
3. **Create new private pySCA-dev** for ongoing development
4. **Push releases** from dev to public pySCA when ready

This way:
- Users see a clean "pySCA" repository
- You have a private dev environment
- Releases are clear and professional

## Summary

**Users will:**
- See repository name: **pySCA** (not pySCA-dev)
- Download release: **pySCA 7.0** (title)
- Get source code: Full repository or release zip/tar.gz
- Install: `pip install -e .` after cloning

**You will:**
- Develop in: Private `pySCA-dev` repository
- Release from: Public `pySCA` repository
- Tag releases: `v7.0`, `v7.1`, etc.
