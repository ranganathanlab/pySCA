# pySCA UX Improvements - Quick Reference

This document provides a quick overview of the UX improvements with concrete examples.

---

## 🚀 Quick Start Improvements

### Before
```bash
# Required: 10+ arguments, must understand all parameters
sca-process-msa alignment.fasta -s 1XYZ --chainID A \
  --precluster --cluster-id 0.85 --cluster-coverage 0.8 \
  --parameters 0.2 0.2 0.2 0.8 --initial-trim-gap 0.8

sca-core Outputs/alignment.db.gz \
  --do-seqcorr --seqcorr-mmseqs2 --do-sector-id \
  --kpos 0 --sector-cutoff 0.95 --lbda 0.03 --norm frob --Ntrials 10
```

### After
```bash
# Option 1: Preset (simplest)
sca-run alignment.fasta -s 1XYZ --preset standard

# Option 2: Config file (reusable)
sca-run alignment.fasta --config my_config.yaml

# Option 3: Auto-detect (smart defaults)
sca-run alignment.fasta -s 1XYZ --auto
```

---

## 📊 Progress Indicators

### Before
```
[INFO] Starting preprocessing...
[INFO] Loading alignment...
# ... long silent wait ...
[INFO] Done
```

### After
```
Processing alignment... [████████████░░░░░░░░] 60%
  - Filtering sequences: ✓ done
  - Computing weights: ✓ done
  - MMseqs2 clustering: [████░░░░] 40%
    Estimated time remaining: 2m 15s
```

---

## ⚙️ Preset Configurations

### Available Presets

#### Quick
```bash
sca-run alignment.fasta -s 1XYZ --preset quick
```
- Fast processing
- Lenient filtering
- Always preclusters
- Skips sequence correlations
- Minimal output

**Use when:** You need fast results, exploring data

#### Standard (Recommended)
```bash
sca-run alignment.fasta -s 1XYZ --preset standard
```
- Balanced parameters
- Auto-detects preclustering need
- Includes sequence correlations
- Includes sector identification
- Standard output

**Use when:** Most analyses, default choice

#### Detailed
```bash
sca-run alignment.fasta -s 1XYZ --preset detailed
```
- Strict filtering
- No preclustering (full accuracy)
- More randomization trials
- Complete analysis
- Maximum accuracy

**Use when:** Publication-quality analysis, small MSAs

---

## 📝 Configuration Files

### Create Config File
```bash
# Save current run as config
sca-process-msa alignment.fasta -s 1XYZ --save-config my_config.yaml
```

### Use Config File
```bash
# Use saved config
sca-run alignment.fasta --config my_config.yaml
```

### Config File Format
```yaml
name: "my_analysis"
description: "Custom configuration"

preprocessing:
  reference:
    pdb: "1XYZ"
    chain: "A"
  filtering:
    max_gap_pos: 0.2
    max_gap_seq: 0.2
  preclustering:
    enabled: true
    cluster_id: 0.85

core:
  regularization: 0.03
  trials: 10
  sequence_correlations: true
  sector_identification: true
```

---

## ✅ Input Validation

### Validate Before Running
```bash
# Check inputs without processing
sca-process-msa alignment.fasta -s 1XYZ --validate

# Output:
✓ Validation passed
  - Alignment file: OK (125,432 sequences, 284 positions)
  - PDB structure: OK (1XYZ Chain A)
  - Disk space: OK (15 GB available)
  - Estimated memory: 2.3 GB (within limits)
```

### Auto-Validation
Validation runs automatically and shows warnings:
```bash
sca-process-msa large_alignment.fasta -s 1XYZ

⚠️  Validation Warning:
  - Very large MSA (250,000 sequences)
  - Consider using --precluster for faster processing
  - Expected time without preclustering: ~8 hours
```

---

## 🎯 Enhanced Error Messages

### Before
```
Error: File not found: 1XYZ
```

### After
```
❌ Error: PDB file not found: 1XYZ

Possible solutions:
  1. Verify PDB ID at https://www.rcsb.org/structure/1XYZ
  2. Check internet connection (PDB files downloaded automatically)
  3. Provide local path: --pdb /path/to/1XYZ.pdb
  4. Use --refseq to provide a reference sequence file instead

See: https://ranganathanlab.gitlab.io/pySCA/usage.html#reference-sequence
```

---

## 📈 Summary Output

### After Processing
```
╔═══════════════════════════════════════════════╗
║         SCA Preprocessing Complete            ║
╠═══════════════════════════════════════════════╣
║ Input sequences:        125,432               ║
║ Output sequences:       24,891                ║
║ Positions:              284                   ║
║ Effective sequences:    18,234                ║
║ Output file:            Outputs/alignment.db.gz ║
║ Processing time:        8m 23s                ║
╚═══════════════════════════════════════════════╝

Next steps:
  Run: sca-core Outputs/alignment.db.gz --do-sector-id
```

---

## 🔍 Help System

### Standard Help
```bash
sca-process-msa --help
```

### Examples
```bash
sca-process-msa --help-examples
# Shows common usage examples
```

### MSA Size Recommendations
```bash
sca-process-msa --help-msa-size
# Shows recommendations for different MSA sizes:
#   - Small (<10k): No preclustering, full analysis
#   - Medium (10k-100k): Optional preclustering
#   - Large (>100k): Recommended preclustering
```

### Parameter Guide
```bash
sca-process-msa --help-params
# Detailed explanation of all parameters
```

---

## 🎨 Interactive Mode

### Guided Setup
```bash
sca-process-msa alignment.fasta --interactive

Welcome to pySCA!
Let's set up your analysis...

Reference sequence method:
  1. PDB structure (recommended)
  2. Reference sequence file
  3. Sequence index
  4. Auto-select
Choice [1]: 1

PDB ID: 1XYZ
Chain ID [A]: A

Analysis type:
  1. Quick (fast, minimal)
  2. Standard (balanced, recommended)
  3. Detailed (comprehensive, slow)
Choice [2]: 2

Large MSA detected (125,432 sequences)
Enable preclustering? [Y/n]: Y

Output name [alignment]: my_analysis

Starting analysis...
```

---

## 📦 Unified Workflow

### Single Command for Complete Analysis
```bash
# Complete workflow: preprocessing → core → sector ID
sca-run alignment.fasta -s 1XYZ --preset standard

# Equivalent to:
# sca-process-msa alignment.fasta -s 1XYZ --preset standard
# sca-core Outputs/alignment.db.gz --preset standard
```

### Workflow Status
```bash
# Check status of current/incomplete workflow
sca-status Outputs/alignment.db.gz

Status: Preprocessing complete
  - Database: Outputs/alignment.db.gz ✓
  - SCA core: Not run
  - Sector ID: Not run

To continue: sca-core Outputs/alignment.db.gz --do-sector-id
```

### Resume Workflow
```bash
# Resume from last successful step
sca-run alignment.fasta -s 1XYZ --resume
```

---

## 🎯 Parameter Recommendations

### Auto-Recommendations
```bash
sca-process-msa alignment.fasta --recommend

Analyzing alignment...
  Sequences: 125,432
  Positions: 284
  Alignment quality: Good

Recommendations:
  ✓ Preclustering: Recommended (>100k sequences)
    Suggested: --precluster --cluster-id 0.85
  
  ✓ Sequence correlations: Recommended
    (MSA size after preclustering will be manageable)
  
  ✓ Memory optimization: Not needed
    (Estimated memory: 2.3 GB)

Suggested command:
  sca-run alignment.fasta -s 1XYZ --preset standard
```

---

## 📁 Output Organization

### Structured Output
```
Outputs/
  alignment_20240115_143022/
    ├── alignment.db.gz           # Main database
    ├── alignment_processed.fasta # Processed alignment
    ├── alignment.log             # Detailed log
    ├── config.yaml               # Configuration used
    └── summary.txt               # Human-readable summary
```

### Summary File
```
# summary.txt
============================================================
SCA Analysis Summary
============================================================
Analysis Date: 2024-01-15 14:30:22

## Input
  Alignment file: alignment.fasta
  Reference: PDB 1XYZ Chain A

## Preprocessing Results
  Input sequences: 125,432
  Final sequences: 24,891
  Positions: 284
  Effective sequences: 18,234

## Output Files
  Database: alignment.db.gz
  Processed alignment: alignment_processed.fasta
  Log: alignment.log
============================================================
```

---

## 🔄 Migration Guide

### Current Users

If you have existing workflows, they still work! New features are optional:

```bash
# Old way (still works)
sca-process-msa alignment.fasta -s 1XYZ --precluster ...

# New way (easier)
sca-run alignment.fasta -s 1XYZ --preset standard
```

### Converting Existing Workflows

1. **Save as config:**
   ```bash
   # Run with your current parameters, save as config
   sca-process-msa alignment.fasta -s 1XYZ [your params] --save-config my_config.yaml
   ```

2. **Use config:**
   ```bash
   sca-run alignment.fasta --config my_config.yaml
   ```

---

## 📚 Key Benefits Summary

| Feature | Benefit |
|---------|---------|
| **Presets** | Start quickly without parameter tuning |
| **Progress bars** | See what's happening during long runs |
| **Config files** | Save and reuse parameter sets |
| **Validation** | Catch errors before processing |
| **Better errors** | Get actionable suggestions |
| **Unified workflow** | One command for complete analysis |
| **Auto-detection** | Smart defaults based on MSA size |
| **Summary output** | Clear view of results |

---

## 🚦 Choosing the Right Approach

### New Users
→ Use `--preset standard` or `--interactive`

### Experienced Users
→ Use config files or explicit parameters

### Batch Processing
→ Use config files for consistency

### Exploratory Analysis
→ Use `--preset quick`

### Publication Analysis
→ Use `--preset detailed` or custom config

---

For detailed information, see:
- `UX_DESIGN.md` - Full design document
- `UX_IMPLEMENTATION_GUIDE.md` - Implementation details
- `USAGE_INSTRUCTIONS.md` - Complete usage guide


