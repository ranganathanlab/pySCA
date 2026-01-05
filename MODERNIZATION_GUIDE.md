# pySCA Modernization Guide

## Overview

This document describes the modernization of `scaTools.py` and associated scripts. The new version (`scaTools_v2.py`) provides:

1. **Better Organization**: Functions grouped by functionality
2. **Modern Python**: Type hints, f-strings, better error handling
3. **All Optimizations**: All performance improvements applied
4. **Backward Compatibility**: Legacy function names maintained as aliases
5. **Better Documentation**: Comprehensive docstrings

## Migration Strategy

### Phase 1: Gradual Migration (Recommended)

1. **Keep both versions**: Use `scaTools_v2.py` alongside original
2. **Test compatibility**: All existing scripts should work with aliases
3. **Update incrementally**: Migrate scripts one at a time

### Phase 2: Full Migration

1. Replace `scaTools.py` with `scaTools_v2.py`
2. Update imports in scripts
3. Test all workflows

## Key Changes

### Function Naming

| Old Name | New Name | Status |
|----------|----------|--------|
| `readAlg` | `read_alignment` | ✅ Alias maintained |
| `lett2num` | `letters_to_numbers` | ✅ Alias maintained |
| `alg2bin` | `alignment_to_binary` | ✅ Alias maintained |
| `seqWeights` | `compute_sequence_weights` | ✅ Alias maintained |
| `filterSeq` | `filter_sequences` | ✅ Alias maintained |
| `freq` | `compute_frequencies` | ✅ Alias maintained |
| `posWeights` | `compute_position_weights` | ✅ Alias maintained |
| `scaMat` | `compute_sca_matrix` | ✅ Alias maintained |
| `randomize` | `randomize_sca` | ✅ Alias maintained |

### Code Organization

**Old Structure:**
- Functions in order of first use
- Mixed concerns
- Inconsistent naming

**New Structure:**
```python
# Constants
# Data Classes
# Alignment I/O
# Alignment Processing
# Statistical Functions
# Eigenvalue Decomposition
# Independent Component Analysis
# SCA Core Functions
# Randomization
# Sector Analysis
# Utility Functions
# Legacy Aliases
```

### Performance Improvements

All optimizations from `OPTIMIZATION_ANALYSIS.md` are included:

1. ✅ Sparse matrix operations (`alg2binss`)
2. ✅ Vectorized calculations (`filterSeq`, `filterPos`)
3. ✅ Efficient frequency computation
4. ✅ Optimized randomization
5. ✅ Better memory management

### New Features

1. **Format Support**: FASTA, Stockholm, Clustal (auto-detect)
2. **Type Hints**: Better IDE support and documentation
3. **Error Handling**: More informative error messages
4. **Constants**: Centralized configuration

## Testing Checklist

- [ ] `scaProcessMSA` works with new version
- [ ] `scaCore` produces identical results
- [ ] `scaSectorID` identifies sectors correctly
- [ ] `annotateMSA` handles all formats
- [ ] All notebook examples run successfully
- [ ] Performance benchmarks match or exceed original

## Next Steps

1. Complete remaining functions (MSAsearch, makeATS, annotation functions, PDB functions)
2. Update bin scripts to use modern practices
3. Add comprehensive unit tests
4. Update documentation


