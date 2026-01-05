# Pseudocode: How `scaMat` Computes `projMat` and `tX`

## Overview

The `scaMat` function computes the SCA positional correlation matrix (`Csca`), the projected alignment (`tX`), and the projector (`projMat`). This document focuses on the computation of `projMat` and `tX`.

## Inputs

- `alg`: M × L numeric alignment (M sequences, L positions)
  - Each element is an integer 0-19 representing amino acids (0 = gap)
- `seqw`: M × 1 vector of sequence weights (default: uniform weights = 1)
- `norm`: Matrix norm type ("frob" or "spec") for correlation matrix
- `lbda`: Regularization parameter (pseudocount)
- `freq0`: 20 × 1 background amino acid frequencies (default: uniform = 1/21)

## Outputs

- `Csca`: L × L SCA positional correlation matrix
- `tX`: M × L projected alignment (dimension-reduced MSA)
- `Proj`: (L × 20) × 1 projector vector (reshaped to `projMat`: L × 20)

## Pseudocode

```
FUNCTION scaMat(alg, seqw, norm, lbda, freq0):
    // ============================================
    // STEP 1: Initialize dimensions
    // ============================================
    M = number of sequences in alg
    L = number of positions in alg
    N_aa = 20  // number of amino acids
    
    // Normalize sequence weights if needed
    IF seqw is integer AND seqw == 1:
        seqw = ones(M)  // uniform weights
    
    // ============================================
    // STEP 2: Compute amino acid frequencies
    // ============================================
    // freq1: (L × 20) × 1 vector of single-position frequencies
    //        freq1[i*20 + a] = frequency of amino acid 'a' at position 'i'
    // freq2: (L × 20) × (L × 20) matrix of pairwise frequencies
    //        freq2[i*20 + a, j*20 + b] = frequency of pair (a at i, b at j)
    freq1, freq2, freq0 = freq(alg, Naa=20, seqw=seqw, lbda=lbda, freq0=freq0)
    
    // ============================================
    // STEP 3: Compute position weights
    // ============================================
    // Wpos: (L × 20) × 1 vector of position-specific weights
    //       Wpos[i*20 + a] = weight for amino acid 'a' at position 'i'
    Wpos = posWeights(alg, seqw, lbda)[0]
    
    // ============================================
    // STEP 4: Compute correlation tensor
    // ============================================
    // tildeC: (L × 20) × (L × 20) correlation tensor
    //         tildeC[i*20 + a, j*20 + b] = weighted correlation between
    //         amino acid 'a' at position 'i' and amino acid 'b' at position 'j'
    tildeC = outer_product(Wpos, Wpos) * (freq2 - outer_product(freq1, freq1))
    
    // ============================================
    // STEP 5: Compute positional correlation matrix (Csca)
    // ============================================
    // This step reduces the (L × 20) × (L × 20) tensor to L × L matrix
    // by taking SVD of each 20 × 20 block and using the spectral or
    // Frobenius norm. This is NOT directly used for projMat/tX, but
    // is computed as part of the SCA matrix.
    
    Cspec = zeros(L, L)
    Cfrob = zeros(L, L)
    
    FOR i = 0 TO L-1:
        FOR j = i TO L-1:
            // Extract 20 × 20 block for positions i and j
            block = tildeC[i*20:(i+1)*20, j*20:(j+1)*20]
            
            // SVD of the block
            u, s, vt = SVD(block)
            
            // Spectral norm (largest singular value)
            Cspec[i, j] = s[0]
            
            // Frobenius norm (sqrt of sum of squared singular values)
            Cfrob[i, j] = sqrt(sum(s^2))
            
            // Make symmetric
            Cspec[j, i] = Cspec[i, j]
            Cfrob[j, i] = Cfrob[i, j]
    
    // Select norm type
    IF norm == "frob":
        Csca = Cfrob
    ELSE:
        Csca = Cspec
    
    // ============================================
    // STEP 6: Compute the projector (Proj)
    // ============================================
    // Proj: (L × 20) × 1 vector
    //       Proj[i*20 + a] = Wpos[i*20 + a] * freq1[i*20 + a]
    //       This is the element-wise product of position weights and frequencies
    Proj = Wpos * freq1  // element-wise multiplication
    
    // ============================================
    // STEP 7: Convert alignment to binary sparse representation
    // ============================================
    // al2d: M × (L × 20) sparse binary matrix
    //       al2d[s, i*20 + a] = 1 if sequence s has amino acid 'a' at position i
    //                          = 0 otherwise
    //       This is a "one-hot" encoding of the alignment
    al2d = alg2binss(alg)
    
    // ============================================
    // STEP 8: Compute projMat and tX
    // ============================================
    // Initialize outputs
    tX = zeros(M, L)           // M × L projected alignment
    projMat = zeros(L, 20)      // L × 20 projector matrix
    
    FOR i = 0 TO L-1:
        // Extract projector for position i (20 × 1 vector)
        Projati = Proj[i*20:(i+1)*20]
        
        // Normalize the projector to unit length
        norm_val = sqrt(sum(Projati^2))
        IF norm_val > 0:
            Projati = Projati / norm_val
        
        // Store normalized projector in projMat
        projMat[i, :] = Projati
        
        // Extract binary representation for position i from al2d
        // al2d_slice: M × 20 sparse matrix
        //             al2d_slice[s, a] = 1 if sequence s has amino acid 'a' at position i
        al2d_slice = al2d[:, i*20:(i+1)*20]
        
        // Project: tX[s, i] = sum over amino acids 'a' of (al2d_slice[s, a] * Projati[a])
        //         This is the dot product of the binary sequence vector with the projector
        dot_result = al2d_slice * Projati  // matrix-vector product
        
        // Convert sparse result to dense if needed
        IF dot_result is sparse:
            tX[:, i] = dot_result.toarray().flatten()
        ELSE:
            tX[:, i] = dot_result.flatten()
    
    // ============================================
    // RETURN
    // ============================================
    RETURN Csca, tX, Proj
    // Note: Proj is (L × 20) × 1, but can be reshaped to L × 20 as projMat
END FUNCTION
```

## Key Mathematical Operations

### 1. Projector Construction
```
Proj[i*20 + a] = Wpos[i*20 + a] × freq1[i*20 + a]
```
- `Wpos`: Position-specific weights (accounts for sequence diversity, gaps, etc.)
- `freq1`: Amino acid frequencies at each position
- The product gives a weighted frequency that emphasizes important positions/amino acids

### 2. Projection Operation
```
tX[s, i] = Σ_a (al2d[s, i*20 + a] × Projati[a])
```
- For each sequence `s` and position `i`:
  - `al2d[s, i*20 + a]` is 1 if sequence `s` has amino acid `a` at position `i`, else 0
  - Multiply by the normalized projector `Projati[a]`
  - Sum over all 20 amino acids
- This reduces the 20-dimensional amino acid space at each position to a single scalar value

### 3. Dimension Reduction
- **Input**: M × L alignment with 20 amino acid states per position → M × (L × 20) binary matrix
- **Output**: M × L matrix with one scalar value per position
- **Reduction**: 20 → 1 dimension per position

## Interpretation

- **`projMat`**: L × 20 matrix where each row is a normalized projector for one position
  - Each row sums to unit length (after normalization)
  - Encodes the "direction" in amino acid space that captures the most information at that position
  
- **`tX`**: M × L projected alignment
  - Each element `tX[s, i]` is a scalar projection of sequence `s` at position `i`
  - Represents how well sequence `s` matches the "preferred" amino acid pattern at position `i`
  - Positive values: sequence matches the pattern
  - Negative values: sequence deviates from the pattern
  - Magnitude: strength of match/deviation

## Relationship to SCA Matrix

The projector `Proj` is derived from the same frequency and weight calculations used to build the correlation tensor `tildeC`, which is then reduced to the SCA matrix `Csca`. The projection `tX` is used in subsequent ICA analysis to identify independent components in sequence space.

