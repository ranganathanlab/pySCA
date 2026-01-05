#!/usr/bin/env python3
"""
pySCA - Statistical Coupling Analysis Toolbox (Modernized Version)

A comprehensive Python implementation of Statistical Coupling Analysis (SCA)
for analyzing protein sequence alignments to identify functionally important
amino acid positions and coevolution patterns.

Authors:
    Olivier Rivoire (olivier.rivoire@ujf-grenoble.fr)
    Kimberly Reynolds (kimberly.reynolds@utsouthwestern.edu)
    Rama Ranganathan (rama.ranganathan@utsouthwestern.edu)

Version: 7.0 (Modernized)
Date: 2024

Copyright (C) 2015-2024 Olivier Rivoire, Rama Ranganathan, Kimberly Reynolds

This program is free software distributed under the BSD 3-clause license.
Please see the file LICENSE for details.

Reference:
    Rivoire, O., Reynolds, K. A., and Ranganathan, R. (2016).
    Evolution-Based Functional Decomposition of Proteins.
    PLOS Computational Biology 12, e1004817.
"""

from __future__ import annotations

import os
import sys
import time
import copy
import random
import colorsys
import sqlite3
import subprocess
from pathlib import Path
from typing import Optional, Union, List, Tuple, Dict, Set, Any
from collections.abc import Sequence

import numpy as np
import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse import csr_matrix, diags
from scipy.stats import t, scoreatpercentile
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm

try:
    from Bio.PDB import PDBParser
    from Bio import pairwise2, SeqIO, Entrez
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False
    print("Warning: BioPython not available. Some functions will be limited.")

from pysca import settings


# ============================================================================
# CONSTANTS
# ============================================================================

# Standard amino acid code
STANDARD_AA_CODE = "ACDEFGHIKLMNPQRSTVWY"
STANDARD_AA_COUNT = 20

# Default background amino acid frequencies (from UniProt)
DEFAULT_AA_FREQ = np.array([
    0.073, 0.025, 0.050, 0.061, 0.042, 0.072, 0.023, 0.053, 0.064, 0.089,
    0.023, 0.043, 0.052, 0.040, 0.052, 0.073, 0.056, 0.063, 0.013, 0.033,
])

# Numerical tolerances
DEFAULT_TOLERANCE = 1e-12
ICA_TOLERANCE = 1e-15


# ============================================================================
# DATA CLASSES
# ============================================================================

class Unit:
    """
    Represents a unit (sector, sequence family, etc.) in the analysis.
    
    Attributes:
        name: String describing the unit (e.g., 'firmicutes')
        items: Set of member item indices
        col: Color code for plotting (0-1 for HSV)
        vect: Additional vector describing member items (e.g., sequence weights)
    """
    
    def __init__(self, name: str = "", items: Optional[Set[int]] = None, 
                 col: float = 0.0, vect: Any = None):
        self.name = name
        self.items = items if items is not None else set()
        self.col = col
        self.vect = vect
    
    def __repr__(self) -> str:
        return f"Unit(name='{self.name}', items={len(self.items)}, col={self.col})"


class Annot:
    """
    Sequence annotation container.
    
    Attributes:
        descr: Description (often the sequence header)
        species: Species string
        taxo: Taxonomy string
        seq: Sequence string
    """
    
    def __init__(self, descr: str, species: str, taxo: str, seq: str = ""):
        self.descr = descr
        self.species = species
        self.taxo = taxo
        self.seq = seq
    
    def __repr__(self) -> str:
        return f"Annot(descr='{self.descr[:30]}...', species='{self.species}')"


class Pair:
    """
    Represents a pair of positions with their coupling information.
    
    Attributes:
        pos: Pair of amino acid positions [i, j]
        DI: Direct information between the two positions
        dist: Physical distance between positions (if available)
    """
    
    def __init__(self, pos: List[int], DI: float, dist: Optional[float] = None):
        self.pos = pos
        self.DI = DI
        self.dist = dist


class Secton:
    """
    Represents a secton (group of positions).
    
    Attributes:
        pos: List of position indices
        num: Number of positions in the secton
    """
    
    def __init__(self, positions: List[int]):
        self.pos = positions
        self.num = len(positions)
    
    def dist(self, distmat: np.ndarray) -> np.ndarray:
        """Return distance matrix between positions in this secton."""
        return distmat[np.ix_(self.pos, self.pos)]
    
    def connected(self, distmat: np.ndarray, threshold: float) -> bool:
        """
        Check structural connectivity using graph theory.
        
        If M_ij is the adjacency matrix, M^n_ij > 0 for n = num_nodes
        indicates i and j are in the same connected component.
        """
        adj = self.dist(distmat) < threshold
        power = np.linalg.matrix_power(adj, self.num)
        return (power > 0).sum() / (self.num ** 2) == 1.0


# ============================================================================
# ALIGNMENT I/O
# ============================================================================

def read_alignment(filename: str, format: Optional[str] = None) -> Tuple[List[str], List[str]]:
    """
    Read multiple sequence alignment from file.
    
    Supports FASTA, Stockholm, and Clustal formats with automatic detection.
    
    Args:
        filename: Path to alignment file
        format: Format type ('fasta', 'stockholm', 'clustal') or None for auto-detect
        
    Returns:
        Tuple of (headers, sequences) where:
            headers: List of sequence headers/names
            sequences: List of sequence strings (uppercase)
            
    Example:
        >>> headers, sequences = read_alignment('alignment.fasta')
        >>> headers, sequences = read_alignment('pfam.stockholm', format='stockholm')
    """
    if format is None:
        format = _detect_alignment_format(filename)
    
    format_readers = {
        'fasta': _read_fasta,
        'stockholm': _read_stockholm,
        'clustal': _read_clustal,
    }
    
    if format not in format_readers:
        raise ValueError(
            f"Unknown format: {format}. "
            f"Supported formats: {list(format_readers.keys())}"
        )
    
    return format_readers[format](filename)


# Backward compatibility alias
readAlg = read_alignment


def _detect_alignment_format(filename: str) -> str:
    """
    Auto-detect alignment format by examining file header.
    
    Returns:
        Format string: 'fasta', 'stockholm', or 'clustal'
    """
    with open(filename, 'r') as f:
        first_lines = [f.readline().rstrip('\n\r') for _ in range(20)]
    
    # Check for Stockholm format (most specific)
    for line in first_lines:
        if line.startswith("# STOCKHOLM"):
            return 'stockholm'
        if line == "//" and any(l.startswith("#") for l in first_lines):
            return 'stockholm'
    
    # Check for Clustal format
    for line in first_lines:
        if line.upper().startswith("CLUSTAL"):
            return 'clustal'
    
    # Additional Clustal detection: conservation lines
    for i, line in enumerate(first_lines):
        if line and not line.startswith(("#", ">")):
            stripped = line.strip()
            if (line.startswith(" ") and len(stripped) > 5 and 
                set(stripped) <= set(" *:.")):
                if i > 0:
                    prev = first_lines[i-1]
                    if prev and not prev.startswith(("#", ">")):
                        parts = prev.split(None, 1)
                        if len(parts) == 2:
                            return 'clustal'
    
    # Default to FASTA
    return 'fasta'


def _read_fasta(filename: str) -> Tuple[List[str], List[str]]:
    """Read FASTA format alignment."""
    headers = []
    sequences = []
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            if line.startswith(">"):
                if current_seq:
                    sequences.append("".join(current_seq).upper())
                    current_seq = []
                headers.append(line[1:])
            elif line:
                current_seq.append(line)
        
        if current_seq:
            sequences.append("".join(current_seq).upper())
    
    return headers, sequences


def _read_stockholm(filename: str) -> Tuple[List[str], List[str]]:
    """
    Read Stockholm format alignment (Pfam format).
    
    Format:
    - Header: "# STOCKHOLM 1.0"
    - Sequence lines: "name    sequence_data" (whitespace separated)
    - Annotation lines start with "#"
    - Ends with "//"
    """
    headers = []
    sequences = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            
            if line.startswith("# STOCKHOLM"):
                continue
            if line == "//":
                break
            if line.startswith("#") or not line.strip():
                continue
            
            # Parse sequence line
            parts = line.split(None, 1)
            if len(parts) >= 2:
                seq_name, seq_data = parts[0], parts[1]
                seq_data = seq_data.replace(" ", "").replace("\t", "")
                
                if seq_data:
                    if seq_name not in sequences:
                        headers.append(seq_name)
                        sequences[seq_name] = []
                    sequences[seq_name].append(seq_data)
    
    # Combine fragments
    seq_list = ["".join(sequences[h]).upper() for h in headers]
    return headers, seq_list


def _read_clustal(filename: str) -> Tuple[List[str], List[str]]:
    """
    Read Clustal format alignment.
    
    Format:
    - Header: "CLUSTAL" or "CLUSTALW"
    - Sequence blocks separated by blank lines
    - Conservation lines with "*", ":", "." symbols
    """
    headers = []
    sequences = {}
    seen_headers = set()
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n\r')
            
            if line.upper().startswith("CLUSTAL") or not line.strip():
                continue
            
            # Skip conservation lines
            stripped = line.strip()
            if (line.startswith(" ") and len(stripped) > 0 and 
                set(stripped) <= set(" *:.")):
                continue
            
            # Parse sequence line
            parts = line.split(None, 1)
            if len(parts) >= 2:
                seq_name, seq_data = parts[0], parts[1]
                seq_data = seq_data.replace(" ", "").replace("\t", "")
                
                if seq_data:
                    if seq_name not in seen_headers:
                        headers.append(seq_name)
                        seen_headers.add(seq_name)
                        sequences[seq_name] = []
                    sequences[seq_name].append(seq_data)
    
    # Combine fragments
    seq_list = ["".join(sequences[h]).upper() for h in headers]
    return headers, seq_list


def parse_alignment_header(header: str, delimiter: str = "|") -> List[str]:
    """
    Parse alignment header with delimiter handling.
    
    Handles cases where delimiter appears inside fields (e.g., in {braces}).
    
    Args:
        header: Header string to parse
        delimiter: Delimiter character (default: "|")
        
    Returns:
        List of header fields
    """
    header_fields = header.split(delimiter)
    idx1 = [i for i, field in enumerate(header_fields) if "{" in field]
    idx2 = [i for i, field in enumerate(header_fields) if "}" in field]
    
    idx_loss = 0
    for i, j in zip(idx1, idx2):
        merged = delimiter.join(header_fields[(i - idx_loss):(j + 1 - idx_loss)])
        header_fields[(i - idx_loss):(j + 1 - idx_loss)] = [merged]
        idx_loss += (j - i)
    
    return header_fields


# Backward compatibility alias
parseAlgHeader = parse_alignment_header


# ============================================================================
# ALIGNMENT PROCESSING
# ============================================================================

def clean_alignment(alignment: List[str], code: str = STANDARD_AA_CODE, 
                   gap: str = "-") -> List[str]:
    """
    Replace invalid amino acid characters with gaps.
    
    Args:
        alignment: List of sequence strings
        code: Valid amino acid characters
        gap: Gap character for replacement
        
    Returns:
        Cleaned alignment with only valid AAs and gaps
    """
    valid_aa = set(code)
    return ["".join(aa if aa in valid_aa else gap for aa in seq) 
            for seq in alignment]


# Backward compatibility alias
clean_al = clean_alignment


def letters_to_numbers(msa_letters: List[str], 
                      code: str = STANDARD_AA_CODE) -> np.ndarray:
    """
    Convert letter-based alignment to numeric representation.
    
    Amino acids 1-20, gaps/other = 0.
    
    Args:
        msa_letters: List of sequence strings
        code: Amino acid code string
        
    Returns:
        Numeric alignment array (Nseq x Npos)
    """
    lett2index = {aa: i + 1 for i, aa in enumerate(code)}
    Nseq, Npos = len(msa_letters), len(msa_letters[0])
    msa_num = np.zeros((Nseq, Npos), dtype=np.int32)
    
    for s, seq in enumerate(msa_letters):
        for i, lett in enumerate(seq):
            if lett in lett2index:
                msa_num[s, i] = lett2index[lett]
    
    return msa_num


# Backward compatibility alias
lett2num = letters_to_numbers


def alignment_to_binary(alg: np.ndarray, N_aa: int = STANDARD_AA_COUNT) -> np.ndarray:
    """
    Convert numeric alignment to binary representation.
    
    Args:
        alg: Numeric alignment (M sequences x L positions)
        N_aa: Number of amino acids (default: 20)
        
    Returns:
        Binary array (M x (N_aa * L))
    """
    N_seq, N_pos = alg.shape
    Abin_tensor = np.zeros((N_aa, N_pos, N_seq), dtype=np.float32)
    
    for ia in range(N_aa):
        Abin_tensor[ia, :, :] = (alg == ia + 1).T
    
    Abin = Abin_tensor.reshape(N_aa * N_pos, N_seq, order="F").T
    return Abin


# Backward compatibility alias
alg2bin = alignment_to_binary


def alignment_to_binary_sparse(alg: np.ndarray, 
                                N_aa: int = STANDARD_AA_COUNT) -> csr_matrix:
    """
    Convert numeric alignment to sparse binary representation.
    
    Optimized version that builds sparse matrix directly.
    
    Args:
        alg: Numeric alignment (M sequences x L positions)
        N_aa: Number of amino acids (default: 20)
        
    Returns:
        Sparse binary matrix (M x (N_aa * L))
    """
    N_seq, N_pos = alg.shape
    rows, cols, data = [], [], []
    
    for ia in range(N_aa):
        matches = (alg == ia + 1)
        row_indices, pos_indices = np.where(matches)
        col_indices = ia * N_pos + pos_indices
        rows.extend(row_indices)
        cols.extend(col_indices)
        data.extend([1.0] * len(row_indices))
    
    return csr_matrix((data, (rows, cols)), shape=(N_seq, N_aa * N_pos))


# Backward compatibility alias
alg2binss = alignment_to_binary_sparse


def compute_sequence_weights(alg: List[str], max_seqid: float = 0.8, 
                            gaps: int = 1) -> np.ndarray:
    """
    Compute sequence weights for alignment.
    
    Weight = 1 / (number of sequences with similarity > max_seqid).
    
    Args:
        alg: List of sequence strings
        max_seqid: Maximum sequence identity threshold
        gaps: If 1, treat gaps as 21st amino acid; if 0, ignore gaps
        
    Returns:
        Sequence weights array (1 x Nseq)
    """
    codeaa = STANDARD_AA_CODE
    if gaps == 1:
        codeaa += "-"
    
    msa_num = letters_to_numbers(alg, code=codeaa)
    Nseq, Npos = msa_num.shape
    X2d = alignment_to_binary(msa_num, N_aa=len(codeaa))
    simMat = X2d.dot(X2d.T) / Npos
    seqw = np.array(1.0 / (simMat > max_seqid).sum(axis=0))
    seqw.shape = (1, Nseq)
    return seqw


# Backward compatibility alias
seqWeights = compute_sequence_weights


def filter_sequences(alg0: List[str], sref: Union[int, float] = 0.5,
                    max_fracgaps: float = 0.2, min_seqid: float = 0.2,
                    max_seqid: float = 0.8) -> Tuple[List[str], np.ndarray, List[int]]:
    """
    Filter alignment sequences by gap fraction and sequence identity.
    
    Args:
        alg0: Input alignment (list of sequences)
        sref: Reference sequence index, or 0.5 for auto-select
        max_fracgaps: Maximum fraction of gaps allowed
        min_seqid: Minimum sequence identity to reference
        max_seqid: Maximum sequence identity for weighting
        
    Returns:
        Tuple of (filtered_alignment, sequence_weights, kept_indices)
    """
    if sref == 0.5:
        sref = choose_reference_sequence(alg0)
    
    Nseq, Npos = len(alg0), len(alg0[0])
    
    # Filter by gap fraction (vectorized)
    gap_counts = np.array([seq.count("-") for seq in alg0])
    seqkeep0 = np.where(gap_counts / Npos < max_fracgaps)[0].tolist()
    
    print(f"Keeping {len(seqkeep0)} of {Nseq} sequences (after filtering for gaps)")
    
    # Filter by sequence identity to reference
    ref_seq = alg0[sref]
    identities = np.array([
        sum(alg0[s][i] == ref_seq[i] for i in range(Npos)) / Npos 
        for s in seqkeep0
    ])
    seqkeep = [seqkeep0[i] for i in np.where(identities > min_seqid)[0]]
    
    print(f"Keeping {len(seqkeep)} of {len(seqkeep0)} sequences "
          f"(after filtering for sequence similarity)")
    
    alg = [alg0[s] for s in seqkeep]
    seqw = compute_sequence_weights(alg, max_seqid)
    
    return alg, seqw, seqkeep


# Backward compatibility alias
filterSeq = filter_sequences


def filter_positions(alg: List[str], seqw: Union[int, List, np.ndarray] = 1,
                    max_fracgaps: float = 0.2) -> Tuple[List[str], List[int]]:
    """
    Filter alignment positions by gap fraction (weighted).
    
    Args:
        alg: Alignment (list of sequences)
        seqw: Sequence weights (default: uniform)
        max_fracgaps: Maximum fraction of gaps allowed at a position
        
    Returns:
        Tuple of (filtered_alignment, kept_position_indices)
    """
    Nseq, Npos = len(alg), len(alg[0])
    
    if isinstance(seqw, int) and seqw == 1:
        seqw = np.ones((1, Nseq))
    elif isinstance(seqw, list):
        seqw = np.array(seqw).reshape(1, -1)
    
    # Compute weighted gap fraction per position
    gaps_mat = np.array([[int(alg[s][i] == "-") for i in range(Npos)] 
                         for s in range(Nseq)])
    seqwn = seqw / seqw.sum()
    gaps_per_pos = seqwn.dot(gaps_mat)[0]
    
    # Select positions
    selpos = np.where(gaps_per_pos < max_fracgaps)[0].tolist()
    
    # Truncate alignment
    alg_tr = ["".join([alg[s][i] for i in selpos]) for s in range(Nseq)]
    
    return alg_tr, selpos


# Backward compatibility alias
filterPos = filter_positions


def choose_reference_sequence(alg: List[str]) -> int:
    """
    Choose reference sequence with mean pairwise identity closest to alignment mean.
    
    Args:
        alg: Alignment (list of sequences)
        
    Returns:
        Index of reference sequence
    """
    if len(alg) > 1000:
        seqw = compute_sequence_weights(alg)
        keep_seq = random_selection(seqw, 1000)
    else:
        keep_seq = list(range(len(alg)))
    
    alg_new = [alg[k] for k in keep_seq]
    num_alg_new = letters_to_numbers(alg_new)
    sim_mat = sequence_similarity(num_alg_new)
    
    # Extract upper triangle
    triu_i, triu_j = np.triu_indices_from(sim_mat, k=1)
    list_s = sim_mat[triu_i, triu_j]
    mean_sid = sim_mat.mean(axis=1)
    mean_diff = np.abs(mean_sid - np.mean(list_s))
    
    strseqnum = [i for i, k in enumerate(mean_diff) if k == min(mean_diff)]
    return keep_seq[strseqnum[0]]


# Backward compatibility alias
chooseRefSeq = choose_reference_sequence


# ============================================================================
# STATISTICAL FUNCTIONS
# ============================================================================

def compute_frequencies(alg: np.ndarray, seqw: Union[int, np.ndarray] = 1,
                        Naa: int = STANDARD_AA_COUNT, lbda: float = 0.0,
                        freq0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute amino acid frequencies for alignment.
    
    Args:
        alg: Numeric alignment (M x L)
        seqw: Sequence weights (1 x M) or 1 for uniform
        Naa: Number of amino acids
        lbda: Pseudocount parameter (0 = no pseudocounts)
        freq0: Background frequencies (default: uniform)
        
    Returns:
        Tuple of (freq1, freq2, freq0_reg) where:
            freq1: Single-site frequencies (Naa*L)
            freq2: Pairwise frequencies (Naa*L x Naa*L)
            freq0_reg: Regularized background frequencies (Naa)
    """
    if freq0 is None:
        freq0 = DEFAULT_AA_FREQ
    
    Nseq, Npos = alg.shape
    
    if isinstance(seqw, int) and seqw == 1:
        seqw = np.ones((1, Nseq))
    
    seqwn = seqw / seqw.sum()
    al2d = alignment_to_binary_sparse(alg, Naa)
    
    # Compute frequencies (keep sparse operations)
    freq1 = np.array(seqwn.dot(al2d)).flatten()
    freq2_sparse = al2d.T.dot(diags(seqwn[0], 0)).dot(al2d)
    freq2 = np.array(freq2_sparse.todense())
    
    # Background model
    block = np.outer(freq0, freq0)
    freq2_bkg = np.tile(block, (Npos, Npos))
    for i in range(Npos):
        freq2_bkg[Naa * i:Naa * (i + 1), Naa * i:Naa * (i + 1)] = np.diag(freq0)
    
    # Regularization
    freq1_reg = (1 - lbda) * freq1 + lbda * np.tile(freq0, Npos)
    freq2_reg = (1 - lbda) * freq2 + lbda * freq2_bkg
    freq0_reg = freq1_reg.reshape(Npos, Naa).mean(axis=0)
    
    return freq1_reg, freq2_reg, freq0_reg


# Backward compatibility alias
freq = compute_frequencies


def compute_position_weights(alg: np.ndarray, seqw: Union[int, np.ndarray] = 1,
                            lbda: float = 0.0, N_aa: int = STANDARD_AA_COUNT,
                            freq0: Optional[np.ndarray] = None,
                            tolerance: float = DEFAULT_TOLERANCE) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute position weights and conservation measures.
    
    Computes:
    - Wia: Position weights (derivative of relative entropy)
    - Dia: Relative entropy per position and amino acid
    - Di: Relative entropy per position
    
    Args:
        alg: Numeric alignment (M x L)
        seqw: Sequence weights
        lbda: Pseudocount parameter
        N_aa: Number of amino acids
        freq0: Background frequencies
        tolerance: Numerical tolerance
        
    Returns:
        Tuple of (Wia, Dia, Di)
    """
    if freq0 is None:
        freq0 = DEFAULT_AA_FREQ
    
    N_seq, N_pos = alg.shape
    
    if isinstance(seqw, int) and seqw == 1:
        seqw = np.ones((1, N_seq))
    
    freq1, freq2, _ = compute_frequencies(alg, seqw=seqw, Naa=N_aa, 
                                          lbda=lbda, freq0=freq0)
    
    # Overall gap fraction
    theta = 1 - freq1.sum() / N_pos
    if theta < tolerance:
        theta = 0.0
    
    # Background with gaps
    freqg0 = (1 - theta) * freq0
    freq0v = np.tile(freq0, N_pos)
    
    # Valid indices (0 < freq < 1)
    iok = np.where((freq1 > 0) & (freq1 < 1))[0]
    
    # Position weights (derivative of relative entropy)
    Wia = np.zeros(N_pos * N_aa)
    Wia[iok] = np.abs(np.log(
        (freq1[iok] * (1 - freq0v[iok])) / ((1 - freq1[iok]) * freq0v[iok])
    ))
    
    # Relative entropies
    Dia = np.zeros(N_pos * N_aa)
    Dia[iok] = (freq1[iok] * np.log(freq1[iok] / freq0v[iok]) + 
                (1 - freq1[iok]) * np.log((1 - freq1[iok]) / (1 - freq0v[iok])))
    
    # Per-position relative entropy
    Di = np.zeros(N_pos)
    for i in range(N_pos):
        freq1i = freq1[N_aa * i:N_aa * (i + 1)]
        aok = np.where(freq1i > 0)[0]
        flogf = freq1i[aok] * np.log(freq1i[aok] / freqg0[aok])
        Di[i] = flogf.sum()
        
        freqgi = 1 - freq1i.sum()
        if freqgi > tolerance:
            Di[i] += freqgi * np.log(freqgi / theta)
    
    return Wia, Dia, Di


# Backward compatibility alias
posWeights = compute_position_weights


def sequence_similarity(alg: np.ndarray) -> np.ndarray:
    """
    Compute sequence similarity matrix.
    
    Args:
        alg: Numeric alignment (M x L)
        
    Returns:
        Similarity matrix (M x M)
    """
    X2d = alignment_to_binary(alg)
    sim_mat = X2d.dot(X2d.T) / alg.shape[1]
    return sim_mat


# Backward compatibility alias
seqSim = sequence_similarity


# ============================================================================
# EIGENVALUE DECOMPOSITION
# ============================================================================

def compute_eigenvectors(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute eigenvectors and eigenvalues of symmetric matrix.
    
    Eigenvalues and vectors are sorted by decreasing eigenvalue.
    Eigenvector signs are fixed so mean component is non-negative.
    
    Args:
        M: Symmetric matrix
        
    Returns:
        Tuple of (eigenvectors, eigenvalues)
    """
    eigenvals, eigenvecs = np.linalg.eigh(M)
    idx = (-eigenvals).argsort()
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]
    
    # Fix sign: mean should be non-negative
    for k in range(eigenvecs.shape[1]):
        mean_val = np.mean(eigenvecs[:, k])
        if np.sign(mean_val) != 0:
            eigenvecs[:, k] = np.sign(mean_val) * eigenvecs[:, k]
    
    return eigenvecs, eigenvals


# Backward compatibility alias
eigenVect = compute_eigenvectors


def sparse_svd(X: csr_matrix, k: int = 6) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Singular value decomposition for sparse matrices (top k components).
    
    Args:
        X: Sparse matrix
        k: Number of components to compute
        
    Returns:
        Tuple of (u, s, v) where X ≈ u.dot(s).dot(v.T)
    """
    u, s, vt = scipy.sparse.linalg.svds(X, k=k)
    idx = (-s).argsort()
    s = s[idx]
    u = u[:, idx]
    
    # Fix sign
    for j in range(u.shape[1]):
        sign = np.sign(np.mean(u[:, j]))
        u[:, j] = sign * u[:, j]
    
    v = X.T.dot(u).dot(np.diag(1.0 / s))
    return u, s, v


# Backward compatibility alias
svdss = sparse_svd


# ============================================================================
# INDEPENDENT COMPONENT ANALYSIS
# ============================================================================

def basic_ica(x: np.ndarray, r0: float, Niter: int, 
              tolerance: float = ICA_TOLERANCE) -> Tuple[np.ndarray, List[float]]:
    """
    Basic ICA algorithm (Infomax, Bell & Sejnowski).
    
    Input data should be sphered (x.T.dot(x) = I).
    
    Args:
        x: Input matrix (L features x M samples)
        r0: Learning rate
        Niter: Number of iterations
        tolerance: Convergence tolerance
        
    Returns:
        Tuple of (unmixing_matrix, change_history)
    """
    L, M = x.shape
    w = np.eye(L)
    change = []
    r = r0 / M
    
    with np.errstate(over="raise"):
        try:
            for iteration in range(Niter):
                w_old = np.copy(w)
                u = w.dot(x)
                w += r * (M * np.eye(L) + 
                         (1.0 - 2.0 / (1.0 + np.exp(-u))).dot(u.T)).dot(w)
                delta = (w - w_old).ravel()
                val = delta.dot(delta.T)
                change.append(val)
                
                if np.isclose(val, 0, atol=tolerance):
                    break
                
                if iteration == Niter - 1:
                    print(f"basicICA failed to converge: {val}")
        except FloatingPointError as e:
            sys.exit(f"Error: basicICA {e}")
    
    return w, change


# Backward compatibility alias
basicICA = basic_ica


def rotate_ica(V: np.ndarray, kmax: int = 6, learnrate: float = 0.1,
               iterations: int = 100000) -> Tuple[np.ndarray, np.ndarray]:
    """
    ICA rotation with normalization.
    
    Args:
        V: Input matrix (N x kmax)
        kmax: Number of components
        learnrate: Learning rate
        iterations: Maximum iterations
        
    Returns:
        Tuple of (V_ica, W) where V_ica are independent components
    """
    V1 = V[:, :kmax].T
    W, changes = basic_ica(V1, learnrate, iterations)
    Vica = (W.dot(V1)).T
    
    # Normalize
    for n in range(kmax):
        imax = abs(Vica[:, n]).argmax()
        Vica[:, n] = (np.sign(Vica[imax, n]) * Vica[:, n] / 
                     np.linalg.norm(Vica[:, n]))
    
    return Vica, W


# Backward compatibility alias
rotICA = rotate_ica


# ============================================================================
# SCA CORE FUNCTIONS
# ============================================================================

def compute_sca_matrix(alg: np.ndarray, seqw: Union[int, np.ndarray] = 1,
                      norm: str = "frob", lbda: float = 0.0,
                      freq0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the SCA correlation matrix.
    
    This is the core SCA computation that produces the positional correlation
    matrix from the alignment.
    
    Args:
        alg: Numeric alignment (M x L)
        seqw: Sequence weights
        norm: Matrix norm type ('frob' or 'spec')
        lbda: Pseudocount parameter
        freq0: Background frequencies
        
    Returns:
        Tuple of (Csca, tX, Proj) where:
            Csca: SCA positional correlation matrix (L x L)
            tX: Projected alignment (M x L)
            Proj: Projector (N_aa*L)
    """
    if freq0 is None:
        freq0 = np.ones(STANDARD_AA_COUNT) / (STANDARD_AA_COUNT + 1)
    
    N_seq, N_pos = alg.shape
    N_aa = STANDARD_AA_COUNT
    
    if isinstance(seqw, int) and seqw == 1:
        seqw = np.ones((1, N_seq))
    
    freq1, freq2, freq0 = compute_frequencies(alg, Naa=N_aa, seqw=seqw, 
                                             lbda=lbda, freq0=freq0)
    Wpos = compute_position_weights(alg, seqw, lbda)[0]
    tildeC = np.outer(Wpos, Wpos) * (freq2 - np.outer(freq1, freq1))
    
    # Positional correlations (this is the bottleneck - O(N_pos^2) SVD operations)
    Cspec = np.zeros((N_pos, N_pos))
    Cfrob = np.zeros((N_pos, N_pos))
    P = np.zeros((N_pos, N_pos, N_aa))
    
    for i in range(N_pos):
        for j in range(i, N_pos):
            block = tildeC[N_aa * i:N_aa * (i + 1), N_aa * j:N_aa * (j + 1)]
            u, s, vt = np.linalg.svd(block)
            Cspec[i, j] = s[0]
            Cfrob[i, j] = np.sqrt(np.sum(s ** 2))
            P[i, j, :] = np.sign(np.mean(u[:, 0])) * u[:, 0]
            P[j, i, :] = np.sign(np.mean(u[:, 0])) * vt[0, :].T
    
    # Make symmetric
    Cspec += np.triu(Cspec, 1).T
    Cfrob += np.triu(Cfrob, 1).T
    
    # Projector
    al2d = alignment_to_binary_sparse(alg)
    tX = np.zeros((N_seq, N_pos))
    Proj = Wpos * freq1
    ProjMat = np.zeros((N_pos, N_aa))
    
    for i in range(N_pos):
        Projati = Proj[N_aa * i:N_aa * (i + 1)]
        norm_val = np.sqrt(np.sum(Projati ** 2))
        if norm_val > 0:
            Projati /= norm_val
        ProjMat[i, :] = Projati
        
        # Use sparse matrix slicing
        tX[:, i] = al2d[:, N_aa * i:N_aa * (i + 1)].dot(Projati.T).toarray().flatten()
    
    if norm == "frob":
        Cspec = Cfrob
    
    return Cspec, tX, Proj


# Backward compatibility alias
scaMat = compute_sca_matrix


def sequence_projection(msa_num: np.ndarray, seqw: np.ndarray,
                       kseq: int = 15, kica: int = 6) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """
    Compute sequence projections with different weighting schemes.
    
    Args:
        msa_num: Numeric alignment
        seqw: Sequence weights
        kseq: Number of eigenvectors
        kica: Number of ICA components
        
    Returns:
        Tuple of (Useq, Uica) where each is a list of 3 projections:
            [0]: No weights
            [1]: Sequence weights
            [2]: Sequence + position weights
    """
    posw, Dia, Di = compute_position_weights(msa_num, seqw)
    Useq = []
    
    # 1 - Raw (no weights)
    X2d = alignment_to_binary_sparse(msa_num)
    Useq.append(sparse_svd(X2d, k=kseq)[0])
    
    # 2 - With sequence weights
    X2dw = csr_matrix(np.diag(np.sqrt(seqw[0]))).dot(X2d)
    u, s, v = sparse_svd(X2dw, k=kseq)
    Useq.append(X2d.dot(v).dot(np.diag(1.0 / s)))
    
    # 3 - With sequence and position weights
    X2dp = X2d.dot(csr_matrix(np.diag(posw)))
    X2dpw = csr_matrix(np.diag(np.sqrt(seqw[0]))).dot(X2dp)
    u, s, v = sparse_svd(X2dpw, k=kseq)
    Useq.append(X2dp.dot(v).dot(np.diag(1.0 / s)))
    
    # Fix signs
    for U in Useq:
        for j in range(U.shape[1]):
            U[:, j] = np.sign(np.mean(U[:, j])) * U[:, j]
    
    # ICA rotation
    Uica = [rotate_ica(U, kmax=kica)[0] for U in Useq]
    
    return Useq, Uica


# Backward compatibility alias
seqProj = sequence_projection


# ============================================================================
# RANDOMIZATION
# ============================================================================

def randomize_alignment(frq: np.ndarray, Mseq: int) -> np.ndarray:
    """
    Generate random alignment preserving position-specific frequencies.
    
    Args:
        frq: Frequency matrix (Npos x (Naa+1)) including gaps
        Mseq: Number of sequences to generate
        
    Returns:
        Random alignment (Mseq x Npos)
    """
    Npos = frq.shape[0]
    msa_rand = np.zeros((Mseq, Npos), dtype=np.int32)
    
    for i in range(Npos):
        Maa = np.random.multinomial(Mseq, frq[i, :])
        col = np.empty(Mseq, dtype=np.int32)
        idx = 0
        for aa, count in enumerate(Maa):
            col[idx:idx+count] = aa
            idx += count
        np.random.shuffle(col)
        msa_rand[:, i] = col
    
    return msa_rand


# Backward compatibility alias
randAlg = randomize_alignment


def randomize_sca(msa_num: np.ndarray, Ntrials: int, seqw: Union[int, np.ndarray] = 1,
                  norm: str = "frob", lbda: float = 0.0, Naa: int = STANDARD_AA_COUNT,
                  kmax: int = 6, tolerance: float = DEFAULT_TOLERANCE) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Randomize alignment and compute SCA spectrum for statistical testing.
    
    Args:
        msa_num: Numeric alignment
        Ntrials: Number of randomization trials
        seqw: Sequence weights
        norm: Matrix norm type
        lbda: Pseudocount parameter
        Naa: Number of amino acids
        kmax: Number of eigenvectors to keep
        tolerance: Numerical tolerance
        
    Returns:
        Tuple of (Vrand, Lrand, Crand) where:
            Vrand: Eigenvectors for each trial (Ntrials x Npos x kmax)
            Lrand: Eigenvalues for each trial (Ntrials x Npos)
            Crand: Average correlation matrix (Npos x Npos)
    """
    Nseq, Npos = msa_num.shape
    
    if isinstance(seqw, int) and seqw == 1:
        seqw = np.ones((1, Nseq))
    
    Mseq = int(np.round(seqw.sum()))
    Crnd = np.zeros((Npos, Npos))
    
    # Compute frequencies including gaps
    f1, f2, f0 = compute_frequencies(msa_num, Naa=20, seqw=seqw, 
                                    lbda=lbda, freq0=np.ones(20) / 21)
    fr1 = np.reshape(f1, (Npos, Naa))
    fr0 = (1.0 - fr1.sum(axis=1)).reshape(Npos, 1)
    
    # Handle roundoff errors
    fr0[fr0 < tolerance] = 0
    fr0[(fr0 > 1) * ((fr0 - tolerance) < 1)] = 1
    fr1[fr1 < tolerance] = 0
    fr1[(fr1 > 1) * ((fr1 - tolerance) < 1)] = 1
    
    fr01 = np.concatenate((fr0, fr1), axis=1)
    
    # Multiple randomizations
    Vrand = np.zeros((Ntrials, Npos, kmax))
    Lrand = np.zeros((Ntrials, Npos))
    
    for t in range(Ntrials):
        msa_rand = randomize_alignment(fr01, Mseq)
        Csca = compute_sca_matrix(msa_rand, norm=norm, lbda=lbda)[0]
        Crnd += Csca
        V, L = compute_eigenvectors(Csca)
        Vrand[t, :, :] = V[:, :kmax]
        Lrand[t, :] = L
    
    Crnd = Crnd / Ntrials
    return Vrand, Lrand, Crnd


# Backward compatibility alias
randomize = randomize_sca


# ============================================================================
# SECTOR ANALYSIS
# ============================================================================

def choose_k_positions(Lsca: np.ndarray, Lrand: np.ndarray) -> int:
    """
    Determine number of significant eigenmodes.
    
    Args:
        Lsca: Eigenvalues of SCA matrix
        Lrand: Eigenvalues from randomized matrices (Ntrials x Npos)
        
    Returns:
        Number of significant modes
    """
    threshold = Lrand[:, 1].mean() + 3 * Lrand[:, 1].std()
    return int(np.sum(Lsca > threshold))


# Backward compatibility alias
chooseKpos = choose_k_positions


def compute_ic_list(Vpica: np.ndarray, kpos: int, Csca: np.ndarray,
                    p_cut: float = 0.95) -> Tuple[List[Unit], List[int], List[int], 
                                                  List[float], List[np.ndarray], List]:
    """
    Identify positions contributing to each independent component.
    
    Args:
        Vpica: Independent components (Npos x kpos)
        kpos: Number of components
        Csca: SCA correlation matrix
        p_cut: Statistical cutoff (CDF of t-distribution)
        
    Returns:
        Tuple of (ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits)
    """
    Npos = len(Vpica)
    cutoff = []
    scaled_pdf = []
    all_fits = []
    
    # Fit t-distribution and compute cutoffs
    for k in range(kpos):
        pd = t.fit(Vpica[:, k])
        all_fits.append(pd)
        iqr = scoreatpercentile(Vpica[:, k], 75) - scoreatpercentile(Vpica[:, k], 25)
        binwidth = 2 * iqr * (len(Vpica[:, k]) ** (-0.33))
        nbins = int(round((np.max(Vpica[:, k]) - np.min(Vpica[:, k])) / binwidth))
        h_params = np.histogram(Vpica[:, k], nbins)
        x_dist = np.linspace(np.min(h_params[1]), np.max(h_params[1]), num=100)
        area_hist = Npos * (h_params[1][2] - h_params[1][1])
        scaled_pdf.append(area_hist * t.pdf(x_dist, pd[0], pd[1], pd[2]))
        cd = t.cdf(x_dist, pd[0], pd[1], pd[2])
        tmp = scaled_pdf[k].argmax()
        
        if abs(np.max(Vpica[:, k])) > abs(np.min(Vpica[:, k])):
            tail = cd[tmp:]
        else:
            cd = 1 - cd
            tail = cd[:tmp]
        
        diff = np.abs(tail - p_cut)
        x_pos = diff.argmin()
        cutoff.append(x_dist[x_pos + tmp])
    
    # Select significant positions
    ic_init = []
    for k in range(kpos):
        ic_init.append(np.where(Vpica[:, k] > cutoff[k])[0].tolist())
    
    # Construct sorted, non-redundant list
    sortedpos = []
    icsize = []
    ics = []
    Csca_nodiag = Csca.copy()
    np.fill_diagonal(Csca_nodiag, 0)
    
    for k in range(kpos):
        icpos_tmp = list(ic_init[k])
        for kprime in [kp for kp in range(kpos) if kp != k]:
            tmp = [v for v in icpos_tmp if v in ic_init[kprime]]
            for i in tmp:
                remsec = (np.linalg.norm(Csca_nodiag[i, ic_init[k]]) < 
                         np.linalg.norm(Csca_nodiag[i, ic_init[kprime]]))
                if remsec:
                    icpos_tmp.remove(i)
        
        sortedpos += sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        icsize.append(len(icpos_tmp))
        
        s = Unit()
        s.items = sorted(icpos_tmp, key=lambda i: -Vpica[i, k])
        s.col = k / kpos
        s.vect = -Vpica[s.items, k]
        ics.append(s)
    
    return ics, icsize, sortedpos, cutoff, scaled_pdf, all_fits


# Backward compatibility alias
icList = compute_ic_list


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def random_selection(seqw: np.ndarray, Mtot: int, keepSeq: List[int] = None) -> List[int]:
    """
    Random selection of sequences with weights, without replacement.
    
    Args:
        seqw: Sequence weights (1 x Nseq)
        Mtot: Total number to select
        keepSeq: Optional list of sequences to always keep
        
    Returns:
        List of selected sequence indices
    """
    if keepSeq is None:
        keepSeq = []
    
    random.seed(0)  # For reproducibility
    return weighted_random_list(seqw[0], Mtot, keepSeq)


# Backward compatibility alias
randSel = random_selection


def weighted_random_list(weights: np.ndarray, Nmax: int, keepList: List[int]) -> List[int]:
    """Generate weighted random list without replacement."""
    Ntot = min(int(np.sum(weights > 0)), Nmax)
    wlist = list(weights)
    selection = list(keepList)
    
    for k in keepList:
        wlist[k] = 0
        Ntot -= 1
    
    for k in range(Ntot):
        i = weighted_random_selection(wlist)
        selection.append(i)
        wlist[i] = 0
    
    return selection


def weighted_random_selection(weights: List[float]) -> int:
    """Select random index with probability given by weights."""
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i
    return len(weights) - 1


# ============================================================================
# LEGACY FUNCTION ALIASES (for backward compatibility)
# ============================================================================

# These maintain backward compatibility with existing code
# All new code should use the modern function names above

# Alignment I/O
readAlg = read_alignment
parseAlgHeader = parse_alignment_header
clean_al = clean_alignment

# Alignment processing
lett2num = letters_to_numbers
alg2bin = alignment_to_binary
alg2binss = alignment_to_binary_sparse
seqWeights = compute_sequence_weights
filterSeq = filter_sequences
filterPos = filter_positions
chooseRefSeq = choose_reference_sequence

# Statistical functions
freq = compute_frequencies
posWeights = compute_position_weights
seqSim = sequence_similarity

# Linear algebra
eigenVect = compute_eigenvectors
svdss = sparse_svd
basicICA = basic_ica
rotICA = rotate_ica

# SCA functions
scaMat = compute_sca_matrix
seqProj = sequence_projection

# Randomization
randAlg = randomize_alignment
randomize = randomize_sca

# Sector analysis
chooseKpos = choose_k_positions
icList = compute_ic_list

# Utilities
randSel = random_selection


