# scaCore_py3.py Analysis and Optimization Plan

## Workflow Overview

The script performs the core SCA calculations:

1. **Load Database** - Reads processed alignment from scaProcessMSA
2. **Sequence Correlations** (optional):
   - Compute sequence similarity matrix (simMat)
   - Compute sequence projections (Useq, Uica)
3. **Positional Weights** - Compute Di, Dia
4. **SCA Matrix** - Compute Csca, tX, Proj
5. **Randomization Trials** - Compute Vrand, Lrand, Crand (Ntrials iterations)
6. **Save Results** - Write database and optional MATLAB file

---

## Identified Issues and Bottlenecks

### 1. ❌ **CRITICAL BUG: Incorrect seqSim() call** (Line 195)

**Problem:**
```python
simMat = sca.seqSim(msa_num, seqw)  # WRONG!
```

The actual signature is:
```python
def seqSim(alg, max_seqs=None, i_ref=None, use_mmseqs2=False, ...)
```

So `seqw` is being interpreted as `max_seqs`, which will cause incorrect behavior!

**Fix:** Remove `seqw` parameter (it's not used in seqSim anyway for the basic call)

---

### 2. ⚠️ **Missing Subsampling Support for Large MSAs**

**Problem:**
- For large MSAs (>50k sequences), `seqSim()` will create huge matrices
- No way to use the new subsampling features we just added

**Solution:** Add command-line options for subsampling:
- `--max-seqcorr-seqs`: Maximum sequences for seqSim subsampling
- `--seqcorr-ref`: Reference sequence index for subsampling
- `--seqcorr-mmseqs2`: Use MMseqs2 for subsampling

---

### 3. ⚠️ **No Progress Tracking for Randomization**

**Problem:**
- `randomize()` runs Ntrials iterations (default 10, can be much more)
- Each trial does full SCA matrix computation
- No progress feedback, user doesn't know how long it will take

**Solution:** Add progress logging inside the randomize loop or before/after calls

---

### 4. ⚠️ **Memory Usage**

**Problem:**
- Large matrices stored as float64 by default
- Multiple large matrices in memory simultaneously
- No automatic memory management

**Solution:** 
- Already has `--float32` option, but could be default for large MSAs
- Better memory usage logging
- Consider memory-efficient storage

---

### 5. ⚠️ **Redundant Type Conversions**

**Problem:**
- Line 212-214: Converts seqw to array and reshapes, but it's already an array from line 173
- Multiple dtype conversions could be optimized

**Solution:** Streamline type handling

---

### 6. ⚠️ **Database I/O**

**Problem:**
- Reads entire database into memory
- Writes entire database (including large matrices)

**Solution:** Already uses gzip compression, but could optimize further

---

## Optimization Strategy

1. **Fix seqSim() call** - Remove incorrect seqw parameter
2. **Add subsampling options** - Enable seqSim subsampling for large MSAs
3. **Add progress tracking** - Better user feedback for long operations
4. **Optimize memory usage** - Better defaults and logging
5. **Streamline code** - Remove redundant operations


