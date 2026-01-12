# How to Update Old pySCA Repository README

## Step 1: Clone or Access the Old pySCA Repository

You have two options:

### Option A: Clone to a Temporary Location

```bash
# Navigate to a temporary directory (not your current pySCA directory)
cd ~/Desktop  # or any location outside your current pySCA

# Clone the old repository
git clone https://github.com/ranganathanlab/pySCA.git pySCA-old

# Navigate into it
cd pySCA-old
```

### Option B: Add as a Remote to Current Directory

```bash
# From your current pySCA directory, add old repo as a remote
cd /Users/ramaranganathan/pySCA

# Add old repository as a remote (use a different name)
git remote add old-pySCA https://github.com/ranganathanlab/pySCA.git

# Fetch from old repository
git fetch old-pySCA

# Checkout the old repository's master branch to a new branch
git checkout -b old-pySCA-readme old-pySCA/master

# Or clone to a separate directory (recommended to avoid confusion)
```

## Step 2: Update the README

### Recommended README Update:

```bash
# Edit the README.md file
# Add this at the very top (before the existing content)
```

**New README content (add at top):**

```markdown
# ⚠️ pySCA (v6.x - Archived)

**This repository has been archived.**

**For the latest version (v7.0), please visit:**
👉 **[pySCA v7.0](https://github.com/ranganathanlab/pySCA)** 👈

---

This version (v6.x) is maintained for historical reference only. 
All new development and releases are in the [pySCA v7.0 repository](https://github.com/ranganathanlab/pySCA).

---

## What's New in v7.0?

- Updated terminology (IC-based naming)
- Improved documentation with Quick Start guide
- Code cleanup and modernization
- Better support for large MSAs

[View pySCA v7.0 →](https://github.com/ranganathanlab/pySCA)

---

*Below is the original README for this archived version:*

---

[Original README content continues here...]
```

## Step 3: Commit and Push

```bash
# If you cloned to a separate directory:
cd ~/Desktop/pySCA-old  # or wherever you cloned it

# Stage the README change
git add README.md

# Commit
git commit -m "Add notice about pySCA v7.0 and archive status"

# Push to old repository
git push origin master
```

## Step 4: Verify

1. Go to `https://github.com/ranganathanlab/pySCA` (old repo)
2. Check that README shows the notice at the top
3. Verify the link to new repository works

## Alternative: Quick Update via GitHub Web Interface

If you prefer not to clone:

1. Go to `https://github.com/ranganathanlab/pySCA`
2. Click on `README.md`
3. Click the **pencil icon** (Edit this file)
4. Add the notice at the top
5. Scroll down, add commit message: "Add notice about pySCA v7.0"
6. Click **Commit changes**

## Template for README Update

Here's a complete template you can use:

```markdown
# ⚠️ pySCA (v6.x - Archived)

<div align="center">

**This repository has been archived.**

**For the latest version (v7.0), please visit:**

# 👉 [pySCA v7.0](https://github.com/ranganathanlab/pySCA) 👈

</div>

---

This version (v6.x) is maintained for historical reference only. 
All new development and releases are in the [pySCA v7.0 repository](https://github.com/ranganathanlab/pySCA).

## What's New in v7.0?

- ✅ Updated terminology (IC-based naming: `ic_cutoff`, `ic_pos`, `ic_ats`)
- ✅ Improved documentation with Quick Start guide
- ✅ Code cleanup and modernization
- ✅ Better support for large MSAs with MMseqs2
- ✅ Enhanced sector identification workflow

**[View pySCA v7.0 →](https://github.com/ranganathanlab/pySCA)**

---

<details>
<summary><b>Original README (v6.x)</b></summary>

[Paste original README content here, or keep existing content below]

</details>

---

*This repository is archived. For the current version, see [pySCA v7.0](https://github.com/ranganathanlab/pySCA).*
```

## Quick Command Summary

```bash
# Clone old repo
cd ~/Desktop
git clone https://github.com/ranganathanlab/pySCA.git pySCA-old
cd pySCA-old

# Edit README.md (add notice at top)
# ... edit file ...

# Commit and push
git add README.md
git commit -m "Add notice about pySCA v7.0 and archive status"
git push origin master
```

## After Updating README

1. **Archive the repository** (Settings → Archive this repository)
2. **Verify** the notice is visible
3. **Test** that links work correctly
