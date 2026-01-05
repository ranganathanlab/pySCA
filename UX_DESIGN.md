# pySCA UX Design Document

## Executive Summary

This document outlines the UX design for pySCA, a command-line tool for Statistical Coupling Analysis (SCA) of protein sequences. The design focuses on improving usability, discoverability, and user guidance while maintaining the power and flexibility of the current system.

---

## Current State Analysis

### Existing Strengths

1. **Comprehensive functionality**: Full SCA workflow with preprocessing, core calculations, and sector identification
2. **Flexible parameters**: Extensive customization options for advanced users
3. **Large MSA support**: Optimizations for handling very large alignments (100k+ sequences)
4. **Good logging**: Structured logging with file output support
5. **Multiple input formats**: Support for FASTA, Stockholm, Clustal

### Current Pain Points

1. **Parameter complexity**: Many parameters with unclear defaults or interactions
2. **No progress indicators**: Long-running operations provide no feedback
3. **Steep learning curve**: Users must understand many options before getting started
4. **Workflow fragmentation**: Multiple commands required for complete analysis
5. **Error messages**: Technical errors without guidance on fixes
6. **No presets**: Must configure everything manually even for common use cases
7. **Parameter discovery**: Difficult to know which parameters matter for a given MSA size
8. **No validation**: Errors discovered late in the workflow
9. **Configuration management**: No way to save/reuse parameter sets

---

## Design Principles

1. **Progressive disclosure**: Simple defaults for beginners, full control for experts
2. **Smart defaults**: Auto-detect optimal settings based on MSA characteristics
3. **Clear feedback**: Progress indicators, status messages, and helpful errors
4. **Workflow guidance**: Validate inputs and suggest corrections
5. **Reusability**: Configuration files for common parameter sets
6. **Discoverability**: Better help text, examples, and documentation

---

## Proposed Improvements

### 1. Command-Line Interface Enhancements

#### 1.1 Preset Configurations

Add preset modes for common use cases:

```bash
# Quick analysis (fast, minimal)
sca-process-msa alignment.fasta -s 1XYZ --preset quick

# Standard analysis (balanced)
sca-process-msa alignment.fasta -s 1XYZ --preset standard

# Detailed analysis (comprehensive, slower)
sca-process-msa alignment.fasta -s 1XYZ --preset detailed
```

**Preset definitions:**
- **quick**: Minimal filtering, no preclustering, basic output
- **standard**: Balanced filtering, auto preclustering for large MSAs, standard output
- **detailed**: Strict filtering, no preclustering, full output including sequence correlations

#### 1.2 Auto-Detection and Recommendations

Automatically analyze MSA and suggest optimal parameters:

```bash
sca-process-msa alignment.fasta -s 1XYZ --auto
# Analyzes: sequence count, alignment quality, memory requirements
# Suggests: preclustering, filtering parameters, memory optimizations
```

#### 1.3 Interactive Mode

Guide new users through parameter selection:

```bash
sca-process-msa alignment.fasta --interactive
# Prompts for:
# - Reference sequence method
# - Analysis goals (quick/standard/detailed)
# - Large MSA optimizations
# - Output options
```

#### 1.4 Improved Help System

Enhanced help with examples and recommendations:

```bash
sca-process-msa --help          # Standard help
sca-process-msa --help-examples # Usage examples
sca-process-msa --help-msa-size # Recommendations by MSA size
sca-process-msa --help-params   # Detailed parameter guide
```

#### 1.5 Unified Workflow Command

Single command for complete workflow:

```bash
sca-run alignment.fasta -s 1XYZ --preset standard
# Runs: preprocessing → core calculations → sector ID
# Single command, single output database
```

### 2. Progress Indicators and Feedback

#### 2.1 Progress Bars

Show progress for long-running operations:

```
Processing alignment... [████████████░░░░░░░░] 60%
  - Filtering sequences: done
  - Computing weights: done
  - MMseqs2 clustering: [████░░░░] 40%
```

#### 2.2 Status Messages

Clear status updates at each stage:

```
[INFO] Starting MSA preprocessing...
[INFO] Input: 125,432 sequences, 284 positions
[INFO] Auto-detected: Large MSA → enabling preclustering
[INFO] Preclustering with MMseqs2... (estimated time: 5 min)
[INFO] Preclustering complete: 25,231 representative sequences
[INFO] Filtering sequences...
[INFO] Filtered: 24,891 sequences remain (98.7% retention)
[INFO] Computing sequence weights...
[INFO] Effective sequences (M_eff): 18,234
[SUCCESS] Preprocessing complete
```

#### 2.3 Estimated Time

Provide time estimates for operations:

```
[INFO] Estimated processing time: 12-15 minutes
[INFO] Memory usage: ~2.3 GB (within limits)
```

#### 2.4 Summary Output

Clear summary at completion:

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

### 3. Configuration Files

#### 3.1 Configuration File Format

YAML-based configuration for reusable parameter sets:

```yaml
# sca_config.yaml
name: "standard_analysis"
description: "Standard SCA analysis with balanced parameters"

preprocessing:
  reference:
    pdb: "1XYZ"
    chain: "A"
  
  filtering:
    max_gap_pos: 0.2
    max_gap_seq: 0.2
    min_SID: 0.2
    max_SID: 0.8
    initial_trim_gap: 0.8
  
  preclustering:
    enabled: true
    cluster_id: 0.85
    cluster_coverage: 0.8
    auto_enable_if_large: true  # Auto-enable if >100k sequences

core:
  regularization: 0.03
  norm: "frob"
  trials: 10
  sequence_correlations: true
  sector_identification: true
  sector_cutoff: 0.95
  
  memory:
    float32: false
    auto_optimize: true  # Use float32 for large MSAs

output:
  matlab: false
  save_msa_numeric: false
```

#### 3.2 Using Configuration Files

```bash
# Use a config file
sca-process-msa alignment.fasta --config sca_config.yaml

# Create config from current run
sca-process-msa alignment.fasta -s 1XYZ --save-config my_config.yaml

# List available configs
sca-config list

# Validate config
sca-config validate sca_config.yaml
```

### 4. Input Validation and Error Handling

#### 4.1 Pre-flight Checks

Validate inputs before processing:

```bash
sca-process-msa alignment.fasta -s 1XYZ --validate
# Checks:
# - File format and readability
# - PDB file accessibility
# - Alignment quality (gaps, sequence length)
# - Memory requirements estimate
# - Disk space availability
```

#### 4.2 Better Error Messages

Provide actionable error messages:

```
❌ Error: PDB file not found: 1XYZ
   Possible causes:
     - PDB ID is incorrect
     - No internet connection (PDB files downloaded automatically)
     - PDB file path is incorrect
   
   Solutions:
     - Verify PDB ID at https://www.rcsb.org/
     - Check internet connection
     - Provide local path: --pdb /path/to/1XYZ.pdb
```

#### 4.3 Workflow Validation

Check that workflow steps are compatible:

```bash
# In sca-core, check that preprocessing was done correctly
[WARN] Database missing 'ats' mapping. Sector identification may be limited.
       Consider re-running sca-process-msa with --pdb option.
```

### 5. Parameter Recommendations

#### 5.1 MSA Size-Based Recommendations

Auto-suggest parameters based on MSA size:

```bash
sca-process-msa alignment.fasta --recommend
# Analyzes MSA and suggests:
# - Preclustering: Recommended (125k sequences > 100k threshold)
# - Memory optimization: Not needed (<50k sequences after preclustering)
# - Sequence correlations: Recommended (manageable size)
```

#### 5.2 Parameter Warnings

Warn about potentially problematic parameter combinations:

```
[WARN] Large MSA (125k sequences) without preclustering.
       This may take several hours and use significant memory.
       Consider: --precluster
```

### 6. Output Organization

#### 6.1 Structured Output Directory

Organize outputs in a clear directory structure:

```
Outputs/
  alignment_20240115_143022/
    ├── alignment.db.gz           # Main database
    ├── alignment_processed.fasta # Processed alignment
    ├── alignment.log             # Log file
    ├── config.yaml               # Configuration used
    └── summary.txt               # Human-readable summary
```

#### 6.2 Summary Reports

Generate human-readable summary reports:

```markdown
# SCA Analysis Summary

**Analysis Date:** 2024-01-15 14:30:22
**Input File:** alignment.fasta
**Reference:** PDB 1XYZ Chain A

## Preprocessing Results
- Input sequences: 125,432
- Final sequences: 24,891 (19.9%)
- Positions: 284
- Effective sequences: 18,234

## Processing Parameters
- Filtering: max_gap_pos=0.2, max_gap_seq=0.2
- Sequence identity: 0.2 - 0.8
- Preclustering: Yes (MMseqs2, 85% identity)

## Next Steps
1. Run core SCA calculations:
   sca-core Outputs/alignment_20240115_143022/alignment.db.gz

2. For full analysis including sector ID:
   sca-core Outputs/alignment_20240115_143022/alignment.db.gz --do-sector-id
```

### 7. Workflow Simplification

#### 7.1 Single-Command Workflow

Unified command for complete analysis:

```bash
# Complete workflow in one command
sca-run alignment.fasta -s 1XYZ --preset standard

# Equivalent to:
# 1. sca-process-msa alignment.fasta -s 1XYZ --preset standard
# 2. sca-core Outputs/alignment.db.gz --do-seqcorr --do-sector-id
```

#### 7.2 Workflow Resumption

Resume interrupted workflows:

```bash
# If interrupted, resume from last successful step
sca-run alignment.fasta -s 1XYZ --resume

# Check workflow status
sca-status Outputs/alignment.db.gz
```

### 8. Documentation Integration

#### 8.1 Context-Sensitive Help

Help integrated into commands:

```bash
# Help for specific parameter
sca-process-msa --help-param precluster

# Help with examples for your use case
sca-process-msa alignment.fasta --help-use-case large-msa
```

#### 8.2 Quick Start Guide

Built-in quick start:

```bash
sca-quickstart
# Interactive guide for first-time users
```

---

## Implementation Roadmap

### Phase 1: Foundation (High Priority)
1. Progress indicators using `tqdm`
2. Enhanced error messages with suggestions
3. Pre-flight validation checks
4. Improved summary output

### Phase 2: Usability (Medium Priority)
5. Preset configurations (quick/standard/detailed)
6. Auto-detection and recommendations
7. Configuration file system
8. Better help system

### Phase 3: Workflow (Lower Priority)
9. Unified workflow command (`sca-run`)
10. Interactive mode
11. Workflow resumption
12. Structured output organization

---

## Example: Before and After

### Before (Current)
```bash
$ sca-process-msa alignment.fasta -s 1XYZ --chainID A --precluster --cluster-id 0.85 --parameters 0.2 0.2 0.2 0.8
[INFO] Starting preprocessing...
[INFO] Loading alignment...
[INFO] Running MMseqs2...
# ... long wait with no feedback ...
[INFO] Done

$ sca-core Outputs/alignment.db.gz --do-seqcorr --seqcorr-mmseqs2 --do-sector-id
[INFO] Starting SCA core...
# ... long wait with no feedback ...
[INFO] Done
```

### After (Proposed)
```bash
$ sca-run alignment.fasta -s 1XYZ --preset standard

╔═══════════════════════════════════════════════╗
║        pySCA Analysis - Standard Preset      ║
╠═══════════════════════════════════════════════╣
║ Input: alignment.fasta                        ║
║ Reference: PDB 1XYZ Chain A                   ║
║ MSA: 125,432 sequences × 284 positions        ║
╚═══════════════════════════════════════════════╝

[1/2] Preprocessing alignment...
  Analyzing MSA... [████████████] 100%
  → Large MSA detected: Auto-enabling preclustering
  
  Preclustering... [████████░░░░] 80%
    - Sequences processed: 100,346 / 125,432
    - Representatives: 24,891
    - Estimated time remaining: 1m 15s
  
  Filtering sequences... [████████████] 100%
  Computing weights... [████████████] 100%
  
  ✓ Preprocessing complete (8m 23s)
    Output: 24,891 sequences, M_eff = 18,234

[2/2] Running SCA core calculations...
  Computing positional weights... [████████████] 100%
  Computing correlation matrix... [████████░░░░] 85%
  Running randomization trials... [█████░░░░░░░] 45%
    Trial 5/10
  Computing sequence correlations... [████░░░░░░] 30%
  Sector identification... [░░░░░░░░░░░░] 0%
  
  ✓ SCA core complete (15m 42s)

╔═══════════════════════════════════════════════╗
║            Analysis Complete!                 ║
╠═══════════════════════════════════════════════╣
║ Results: Outputs/alignment_20240115_143022/   ║
║   - alignment.db.gz (main database)           ║
║   - summary.txt (readable summary)            ║
║   - alignment.log (detailed log)              ║
║                                                ║
║ Total time: 24m 5s                            ║
║ Total sectors identified: 3                   ║
╚═══════════════════════════════════════════════╝

Next: Open summary.txt or view results in Jupyter notebook
```

---

## Success Metrics

1. **Usability**
   - Time to first successful run: < 5 minutes for new users
   - Number of command-line arguments needed for common case: < 3
   - User errors caught before processing: > 80%

2. **Feedback**
   - Progress visibility: All operations > 30s show progress
   - Error message quality: All errors include actionable suggestions
   - Time estimates: Available for operations > 1 minute

3. **Adoption**
   - Configuration file usage: > 50% of regular users
   - Preset usage: > 70% of new users
   - Workflow completion rate: > 90% for standard analyses

---

## Appendix: Command Reference Comparison

### Current Commands
```bash
# Preprocessing
sca-process-msa alignment.fasta -s 1XYZ --chainID A \
  --precluster --cluster-id 0.85 \
  --parameters 0.2 0.2 0.2 0.8 \
  --initial-trim-gap 0.8

# Core analysis
sca-core Outputs/alignment.db.gz \
  --do-seqcorr --seqcorr-mmseqs2 \
  --do-sector-id --kpos 0 \
  --sector-cutoff 0.95 \
  --lbda 0.03 --norm frob --Ntrials 10
```

### Proposed Commands
```bash
# Option 1: Preset (simple)
sca-run alignment.fasta -s 1XYZ --preset standard

# Option 2: Config file (reusable)
sca-run alignment.fasta --config my_standard_config.yaml

# Option 3: Explicit (advanced users)
sca-run alignment.fasta -s 1XYZ \
  --precluster \
  --do-seqcorr --do-sector-id
  # ... (auto-detects optimal sub-parameters)
```

---

## Conclusion

This UX design focuses on making pySCA more accessible while preserving its power and flexibility. By providing smart defaults, clear feedback, and workflow guidance, we can significantly improve the user experience for both new and experienced users.

The proposed improvements can be implemented incrementally, starting with high-impact, low-effort changes (progress indicators, better errors) and gradually adding more sophisticated features (presets, configuration files, unified workflows).


