# Syncing pySCA Between Multiple Computers

## Overview

To work on pySCA from multiple computers, you'll use Git to sync changes through the GitHub repository. Cursor AI sessions are local to each computer, but your code changes are synced via Git push/pull.

## Initial Setup on Your Other Computer

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

# Install dependencies
pip install -e ".[notebooks]"

# Install external tools (FASTA36, MMseqs2)
# Follow instructions in INSTALLATION.md
```

### 3. Configure Git (if not already done)

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

### 4. Set Up Authentication

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

