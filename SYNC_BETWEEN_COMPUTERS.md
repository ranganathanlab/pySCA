# Syncing pySCA Between Multiple Computers

## Overview

To work on pySCA from multiple computers, you'll use Git to sync changes through the GitHub repository. Cursor AI sessions are local to each computer, but your code changes are synced via Git push/pull.

## Initial Setup on Your Other Computer

### 0. Install Miniconda (if not already installed)

If you don't have conda installed on your other computer:

**macOS:**
```bash
# Download Miniconda installer
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh

# Run installer
bash Miniconda3-latest-MacOSX-x86_64.sh

# Follow the prompts, then restart your terminal or run:
source ~/.zshrc  # or ~/.bash_profile
```

**Linux:**
```bash
# Download Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run installer
bash Miniconda3-latest-Linux-x86_64.sh

# Follow the prompts, then restart your terminal or run:
source ~/.bashrc
```

**Verify installation:**
```bash
conda --version
# Should show something like: conda 23.x.x
```

### 1. Clone the Repository

```bash
# Navigate to where you want the project (e.g., your home directory)
cd ~

# Clone the repository
git clone https://github.com/ranganathanlab/pySCA-dev.git pySCA

# Or if you prefer a different directory name:
# git clone https://github.com/ranganathanlab/pySCA-dev.git ~/cAI/pySCA
```

### 2. Set Up Python Environment

```bash
cd pySCA

# Create conda environment (same as on your main computer)
conda create -n pysca3 python=3.10 -y
conda activate pysca3

# Install pySCA and Python dependencies
pip install -e ".[notebooks]"
```

### 3. Install External Tools

You'll need to install two external command-line tools. These are **system-level** installations (not Python packages), so they work outside the conda environment.

#### FASTA36 (Required)

FASTA36 provides the `ggsearch36` program, which is required for reference sequence searching in `sca-process-msa`.

**macOS (via Homebrew - Recommended):**
```bash
# Install Homebrew if you don't have it: https://brew.sh
brew install fasta

# Verify installation
which ggsearch36
ggsearch36 -h
```

**macOS (from source):**
```bash
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make -j2 -f ../make/Makefile.osx all
sudo cp -r ../bin /usr/local
sudo rm /usr/local/bin/README
cd ../..

# Verify installation
which ggsearch36
ggsearch36 -h
```

**Linux (from source):**
```bash
git clone https://github.com/wrpearson/fasta36.git
cd fasta36/src
make -j2 -f ../make/Makefile.linux all
sudo cp -r ../bin /usr/local
sudo rm /usr/local/bin/README
cd ../..

# Verify installation
which ggsearch36
ggsearch36 -h
```

**If `ggsearch36` is not found after installation:**
- Make sure `/usr/local/bin` is in your PATH
- Check: `echo $PATH | grep /usr/local/bin`
- If not, add it to your shell config file (`~/.zshrc` or `~/.bashrc`)

#### MMseqs2 (Optional but Recommended)

MMseqs2 is highly recommended for preclustering large MSAs (>50k sequences). It significantly speeds up processing.

**macOS/Linux (via Conda - Recommended):**
```bash
# Make sure you're in the pysca3 environment
conda activate pysca3
conda install -c bioconda mmseqs2

# Verify installation
which mmseqs
mmseqs version
```

**macOS (via Homebrew):**
```bash
brew install mmseqs2

# Verify installation
which mmseqs
mmseqs version
```

**Linux (via Package Manager):**
```bash
# Ubuntu/Debian
sudo apt-get install mmseqs2

# Verify installation
which mmseqs
mmseqs version
```

**Note:** If you install MMseqs2 via conda, it will only be available when the `pysca3` environment is activated. If you install via Homebrew or system package manager, it will be available system-wide.

For more detailed installation instructions and troubleshooting, see `INSTALLATION.md` in the repository.

### 4. Configure Git (if not already done)

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### 5. Set Up Authentication

You'll need to authenticate with GitHub. Use the same Personal Access Token method:

1. Use the same token you created earlier, OR
2. Create a new token at https://github.com/settings/tokens
3. Configure credential helper:
   ```bash
   git config --global credential.helper osxkeychain  # macOS
   # or
   git config --global credential.helper store       # Linux
   ```

## Daily Workflow: Syncing Changes

### On Computer A (this computer):

**Before starting work:**
```bash
cd /Users/rama/cAI/pySCA
git pull dev master  # or: git pull origin master (if you set up tracking)
```

**After making changes:**
```bash
# Stage changes
git add .

# Commit
git commit -m "Description of changes"

# Push to GitHub
git push dev master
```

### On Computer B (other computer):

**Before starting work:**
```bash
cd ~/pySCA  # or wherever you cloned it
git pull origin master  # Pull latest changes from GitHub
```

**After making changes:**
```bash
# Stage changes
git add .

# Commit
git commit -m "Description of changes"

# Push to GitHub
git push origin master
```

**Then back on Computer A:**
```bash
cd /Users/rama/cAI/pySCA
git pull dev master  # Get the changes from Computer B
```

## Setting Up Remote Tracking (Recommended)

To make pushing/pulling easier, set up branch tracking:

### On Computer A (this computer):
```bash
cd /Users/rama/cAI/pySCA
git branch --set-upstream-to=dev/master master
# Now you can just use: git pull and git push
```

### On Computer B (other computer):
```bash
cd ~/pySCA
git branch --set-upstream-to=origin/master master
# Now you can just use: git pull and git push
```

## Important Notes

### Cursor AI Sessions
- **Cursor AI sessions are local** - they don't automatically sync
- **Code changes sync via Git** - commit and push to share
- Each computer has its own Cursor AI session and history

### Best Practices

1. **Always pull before starting work:**
   ```bash
   git pull
   ```

2. **Commit frequently:**
   ```bash
   git add .
   git commit -m "Brief description"
   git push
   ```

3. **Check status before pushing:**
   ```bash
   git status
   git log --oneline -5  # See recent commits
   ```

4. **Handle conflicts if they occur:**
   - If both computers modified the same file, Git will mark conflicts
   - Resolve conflicts manually, then:
     ```bash
     git add .
     git commit -m "Resolved merge conflicts"
     git push
     ```

### Working Directory Structure

You can clone to the same path on both computers for consistency:
- **Computer A:** `/Users/rama/cAI/pySCA`
- **Computer B:** `/Users/rama/cAI/pySCA` (if username is the same)

Or use different paths:
- **Computer A:** `/Users/rama/cAI/pySCA`
- **Computer B:** `~/projects/pySCA`

Git doesn't care about the local path - it syncs through the remote repository.

## Quick Reference

```bash
# Pull latest changes
git pull

# Check what changed
git status
git log --oneline -5

# Commit and push changes
git add .
git commit -m "Your message"
git push

# See remote repositories
git remote -v
```

## Troubleshooting

### "Your branch is ahead of 'origin/master'"
- You have local commits that haven't been pushed
- Run: `git push`

### "Your branch is behind 'origin/master'"
- Remote has changes you don't have locally
- Run: `git pull`

### "Updates were rejected because the remote contains work"
- Someone (or you on another computer) pushed changes
- Run: `git pull` first, then `git push`

### Authentication fails
- Make sure you've set up the Personal Access Token
- Check: `git config --global credential.helper osxkeychain`
- Re-authenticate if needed

