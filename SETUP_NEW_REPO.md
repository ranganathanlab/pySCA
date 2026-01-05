# Setting up New Development Repository

## Option 1: Create new repo on ranganathanlab (Recommended)

1. Go to https://github.com/organizations/ranganathanlab/repositories/new
2. Repository name: `pySCA-dev` or `pySCA-7.0`
3. Description: "Development version of pySCA 7.0 with modernizations"
4. Visibility: Private (for development) or Public (if ready)
5. Do NOT initialize with README, .gitignore, or license (we already have these)

Then run:
```bash
cd /Users/rama/cAI/pySCA
git remote add dev https://github.com/ranganathanlab/pySCA-dev.git
# OR if using personal account:
git remote add dev https://github.com/ramaranganathan/pySCA-dev.git
```

## Option 2: Use personal account (temporary)

1. Go to https://github.com/new
2. Repository name: `pySCA-dev`
3. Description: "Development version of pySCA 7.0"
4. Visibility: Private or Public
5. Do NOT initialize with README, .gitignore, or license

Then run:
```bash
cd /Users/rama/cAI/pySCA
git remote add dev https://github.com/ramaranganathan/pySCA-dev.git
```

## After creating the repo, push your code:

```bash
# Stage all changes
git add .

# Commit
git commit -m "pySCA 7.0: Modernized version with robustica integration, improved logging, and documentation"

# Push to new dev repo
git push dev master

# Or create a new branch for development:
git checkout -b dev-7.0
git push dev dev-7.0
```

## On your other computer:

```bash
# Clone the new dev repo
git clone https://github.com/ranganathanlab/pySCA-dev.git pySCA
# OR
git clone https://github.com/ramaranganathan/pySCA-dev.git pySCA

cd pySCA
# Set up environment
conda create -n pysca3 python=3.11
conda activate pysca3
pip install -e .
pip install -e ".[notebooks]"
```

## When ready to replace public repo:

1. Merge dev branch into master
2. Update the existing ranganathanlab/pySCA repo
3. Or archive old repo and make dev repo the new public one
