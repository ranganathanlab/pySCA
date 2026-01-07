#!/usr/bin/env python3
"""
scaProcessMSA_py3_big.py

Large-MSA-safe SCA pre-processing pipeline with optional MMseqs2 preclustering.

Improvements vs earlier version:
  - Structured logging to both console and optional log file (--log)
  - Consistent final summary (M, M', L) regardless of weighting path
  - Explicit Matlab output path and confirmation message
  - Safer handling of --output (always treated as a basename)
  - Optional verbosity controls (--quiet / --verbose)

Input formats:
  FASTA / Stockholm / Clustal (via Biopython AlignIO inside mmseqs_precluster_msa.py and scaTools.readAlg)

Requires:
  numpy, scipy, biopython
  mmseqs (if --precluster)
  scaTools (Python 3 version)

"""
from __future__ import annotations

import argparse
import gzip
import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.io import savemat

# IMPORTANT: package-friendly import
from pysca import scaTools as sca

# BioPython for pairwise protein alignment
try:
    from Bio.Align import PairwiseAligner
    HAS_BIOALIGN = True
except ImportError:
    HAS_BIOALIGN = False

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWY-")


def setup_logger(log_path: Optional[str], verbose: bool, quiet: bool) -> logging.Logger:
    logger = logging.getLogger("sca-process-msa")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # console handler
    ch = logging.StreamHandler()
    if quiet:
        ch.setLevel(logging.WARNING)
    elif verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)

    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    # file handler
    if log_path:
        fp = Path(log_path)
        fp.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(fp, mode="w", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.info(f"Logging to file: {fp}")

    return logger


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def load_weights_tsv(path: str) -> Dict[str, float]:
    """
    Load sequence weights from TSV file.
    
    Format: sequence_id <TAB> weight
    Optimized for large files with efficient parsing.
    """
    w: Dict[str, float] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Use split with maxsplit=1 for efficiency (only split on first tab)
            parts = line.split("\t", 1)
            if len(parts) >= 2:
                sid, val = parts[0], parts[1]
                try:
                    w[sid] = float(val)
                except ValueError:
                    # Skip invalid weight values
                    continue
    return w


def filter_nonstandard(headers: List[str], seqs: List[str]) -> Tuple[List[str], List[str]]:
    """
    Filter sequences to keep only those with standard amino acids.
    
    Optimized for large MSAs using vectorized operations where possible.
    """
    if not seqs:
        return [], []
    
    # Convert to uppercase once (more efficient than per-sequence)
    seqs_upper = [s.upper() for s in seqs]
    
    # Vectorized filtering: check if all characters in each sequence are allowed
    # For very large alignments, this is faster than set operations per sequence
    valid_mask = []
    for sU in seqs_upper:
        # Use set intersection for fast check (O(1) average for set operations)
        if set(sU) <= AA_ALLOWED:
            valid_mask.append(True)
        else:
            valid_mask.append(False)
    
    # Use list comprehension with zip for efficient filtering
    hd_out = [h for h, valid in zip(headers, valid_mask) if valid]
    alg_out = [sU for sU, valid in zip(seqs_upper, valid_mask) if valid]
    
    return hd_out, alg_out


def write_fasta(path: Path, headers: List[str], seqs: List[str]) -> None:
    """
    Write alignment to FASTA file.
    
    Optimized for large MSAs by batching writes and using efficient string formatting.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    
    # For very large alignments, batch the writes for better I/O performance
    # Build lines in memory first, then write in chunks
    if len(seqs) > 10000:
        # Batch writing for large files
        with path.open("w", encoding="utf-8", buffering=8192*4) as f:  # Larger buffer
            batch_lines = []
            for h, s in zip(headers, seqs):
                batch_lines.append(f">{h}\n{s}\n")
                # Write in batches to reduce I/O overhead
                if len(batch_lines) >= 1000:
                    f.writelines(batch_lines)
                    batch_lines = []
            # Write remaining lines
            if batch_lines:
                f.writelines(batch_lines)
    else:
        # For smaller files, simple approach is fine
        with path.open("w", encoding="utf-8") as f:
            for h, s in zip(headers, seqs):
                f.write(f">{h}\n{s}\n")


def parse_refpos_file(path: str) -> List[str]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return [line.strip() for line in f if line.strip()]


def run_precluster(alignment: str, out_prefix: Path, cluster_id: float, coverage: float, cov_mode: int,
                   keep_tmp: bool, logger: logging.Logger, keep_sequence_id: Optional[str] = None,
                   keep_sequence_ids: Optional[List[str]] = None) -> Tuple[str, str]:
    """
    Runs mmseqs_precluster_msa.py (expected alongside this script) and returns:
      reps_fasta_path, weights_tsv_path
    
    Args:
        keep_sequence_id: If provided, ensures this sequence ID is always included as a representative
        keep_sequence_ids: If provided, list of sequence IDs to always include as representatives
    """
    helper = Path(__file__).with_name("mmseqs_precluster_msa.py")
    if not helper.exists():
        raise FileNotFoundError(f"Expected helper script next to this file: {helper}")

    cmd = [
        "python3", str(helper),
        alignment,
        "--out-prefix", str(out_prefix),
        "--min-seq-id", str(cluster_id),
        "-c", str(coverage),
        "--cov-mode", str(cov_mode),
    ]
    if keep_tmp:
        cmd.append("--keep-tmp")
    # Combine keep_sequence_id and keep_sequence_ids into a single list
    all_keep_ids = []
    if keep_sequence_id:
        all_keep_ids.append(keep_sequence_id)
    if keep_sequence_ids:
        all_keep_ids.extend(keep_sequence_ids)
    
    # Pass all sequence IDs to retain (mmseqs_precluster_msa.py accepts multiple --keep-sequence-id flags)
    for seq_id in all_keep_ids:
        cmd.extend(["--keep-sequence-id", seq_id])
    
    if all_keep_ids:
        logger.info(f"Ensuring {len(all_keep_ids)} sequence(s) are retained in preclustered alignment: {all_keep_ids}")

    logger.info("=" * 80)
    logger.info("Using MMseqs2 for sequence preclustering")
    logger.info(f"  Identity threshold: {cluster_id}")
    logger.info(f"  Coverage threshold: {coverage}")
    logger.info(f"  Coverage mode: {cov_mode}")
    logger.info("=" * 80)
    # Suppress MMseqs2 verbose output - only show summary information
    result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    # Extract and log only useful summary information from stdout
    if result.stdout:
        stdout_lines = result.stdout.strip().split('\n')
        for line in stdout_lines:
            line_lower = line.lower()
            # Log only key summary lines (clusters, representatives, completion messages)
            if any(keyword in line_lower for keyword in ['clusters:', 'representative', 'wrote:']):
                # Clean up the line and log it
                clean_line = line.strip()
                if clean_line:
                    logger.info(f"  {clean_line}")
    # Only log stderr if there are actual errors (not warnings)
    if result.stderr:
        stderr_lines = result.stderr.strip().split('\n')
        error_lines = [line.strip() for line in stderr_lines 
                      if line.strip() and ('error' in line.lower() or 'failed' in line.lower())]
        if error_lines:
            for line in error_lines:
                logger.error(f"MMseqs2 error: {line}")

    reps_fa = str(out_prefix) + "_reps.fasta"
    wts = str(out_prefix) + "_weights.tsv"
    return reps_fa, wts


def main(argv: Optional[List[str]] = None) -> int:
    """
    Main workflow for SCA MSA preprocessing.
    
    Workflow steps:
    1. Load and validate input alignment
    2. Optionally run MMseqs2 preclustering (for large MSAs >50k sequences)
    3. Filter non-standard amino acids
    4. Initial trim of highly gapped positions
    5. Identify reference sequence (PDB, file, or auto-selected)
    6. Create ATS (Alignment-To-Structure) mapping (if PDB/reference provided)
    7. Filter sequences by gap fraction and sequence identity
    8. Filter positions by gap fraction
    9. Compute sequence weights (MMseqs2 cluster sizes or scaTools.seqWeights)
    10. Write processed alignment and database
    
    Key features:
    - Automatically enables preclustering for alignments >50k sequences
    - Ensures reference sequence is always retained (at its original position in input MSA)
    - Normalizes MMseqs2 cluster-size weights to sum = number of clusters
    - Supports multiple alignment formats (FASTA, Stockholm, Clustal)
    
    Returns:
        int: Exit code (0 for success)
    """
    p = argparse.ArgumentParser()
    p.add_argument("alignment", help="Input alignment (FASTA/Stockholm/Clustal; optionally .gz).")

    # Reference options
    p.add_argument("-s", "--pdb", dest="pdbid", help="PDB identifier or path (passed to scaTools.pdbSeq)")
    p.add_argument("-c", "--chainID", dest="chainID", default="A", help="Chain ID in the PDB for the reference sequence")
    p.add_argument("-f", "--species", dest="species", help="Species of reference (for MSAsearch heuristic)")
    p.add_argument("-r", "--refseq", dest="refseq", help="Reference sequence FASTA file")
    p.add_argument("-o", "--refpos", dest="refpos", help="Reference positions file")
    p.add_argument("-i", "--refindex", dest="i_ref", type=int, help="Reference sequence index (0-based)")

    p.add_argument("-p", "--parameters", dest="parameters", default=[0.3, 0.2, 0.15, 0.85],
                   type=float, nargs=4,
                   help="Filtering params: [max_gap_pos, max_gap_seq, min_SID, max_SID]")
    p.add_argument("-t", "--truncate", action="store_true", default=False, help="Truncate to PDB positions (if using PDB)")
    p.add_argument("-m", "--matlab", action="store_true", default=False, help="Also write Matlab .mat workspace")
    p.add_argument("--output", dest="outputfile", default=None, help="Output base name (basename only)")

    # Big-MSA controls
    p.add_argument("--precluster", action="store_true", default=None,
                   help="Run MMseqs2 preclustering and use cluster reps + cluster-size weights. "
                        "Automatically enabled for alignments with >50,000 sequences.")
    p.add_argument("--no-precluster", dest="precluster", action="store_false",
                   help="Disable automatic preclustering even for large alignments.")
    p.add_argument("--cluster-id", type=float, default=0.85, help="MMseqs2 identity threshold (default 0.85).")
    p.add_argument("--cluster-coverage", type=float, default=0.8, help="MMseqs2 coverage threshold (default 0.8).")
    p.add_argument("--cluster-cov-mode", type=int, default=0, choices=[0, 1, 2], help="MMseqs2 cov-mode (default 0).")
    p.add_argument("--keep-mmseqs-tmp", action="store_true", default=False, help="Keep mmseqs tmp dir for debugging.")
    p.add_argument("--keep-sequences", type=int, nargs="+", default=None,
                   help="Indices (1-based) of sequences from the full MSA to always retain after MMseqs2 clustering. "
                        "First sequence is 1, second is 2, etc. These sequences will be added as representatives "
                        "if not already selected by MMseqs2.")
    p.add_argument("--keep-sequences-file", type=str, default=None,
                   help="Text file containing sequence indices (1-based) to retain, one per line or space-separated. "
                        "Indices in the file are 1-based (first sequence is 1, second is 2, etc.). "
                        "Can be combined with --keep-sequences.")

    p.add_argument("--initial-trim-gap", type=float, default=None, help="Initial trim threshold for position gaps. If not specified, defaults to 0.8. Filters positions with gap fraction >= this value.")
    p.add_argument("--save-msa-num", action="store_true", default=False, help="Save msa_num to compressed npz.")

    # Logging / verbosity
    p.add_argument("--log", default=None, help="Write a full debug log to this file (e.g., Outputs/run.log).")
    p.add_argument("--verbose", action="store_true", default=False, help="Verbose console logging.")
    p.add_argument("--quiet", action="store_true", default=False, help="Only warnings/errors to console.")

    args = p.parse_args(argv)
    
    # If initial_trim_gap not specified, use default of 0.8
    initial_trim_was_calculated = False
    if args.initial_trim_gap is None:
        args.initial_trim_gap = 0.8
        initial_trim_was_calculated = True
    
    logger = setup_logger(args.log, args.verbose, args.quiet)
    
    # Configure scaTools logger to use the same handlers (so ggsearch36 output goes to log file)
    sca_logger = logging.getLogger("pysca.scaTools")
    sca_logger.setLevel(logging.DEBUG)
    sca_logger.handlers.clear()  # Remove any existing handlers
    sca_logger.propagate = False  # Don't propagate to root logger (we'll handle it)
    # Add the same handlers as the main logger
    for handler in logger.handlers:
        sca_logger.addHandler(handler)

    ensure_dir(Path("Inputs"))
    ensure_dir(Path("Outputs"))

    in_path = Path(args.alignment)
    if not in_path.exists():
        raise FileNotFoundError(str(in_path))

    # Sanitize output name (basename only)
    base = Path(args.outputfile).name if args.outputfile else in_path.stem
    work_prefix = Path("Outputs") / base

    logger.info("=== sca-process-msa run ===")
    logger.info(f"Input alignment: {in_path}")
    logger.info(f"Output base: {base}")
    logger.info(f"Parameters: pos_gap={args.parameters[0]}, seq_gap={args.parameters[1]}, minSID={args.parameters[2]}, maxSID={args.parameters[3]}")
    if initial_trim_was_calculated:
        logger.info(f"Initial trim gap threshold (position gaps): {args.initial_trim_gap} (default)")
    else:
        logger.info(f"Initial trim gap threshold (position gaps): {args.initial_trim_gap} (user-specified)")

    # Quick check of alignment size to auto-enable preclustering if needed
    # This prevents users from accidentally running O(N²) seqWeights on huge alignments
    if args.precluster is None:
        logger.info("Checking alignment size to determine if preclustering is needed...")
        try:
            # Quick read just to count sequences (we'll reload after preclustering if needed)
            headers_check, _ = sca.readAlg(str(in_path))
            n_seqs = len(headers_check)
            if n_seqs > 50000:
                args.precluster = True
                logger.info(f"Large alignment detected ({n_seqs} sequences). "
                           f"Automatically enabling MMseqs2 preclustering for efficiency. "
                           f"Use --no-precluster to disable.")
            else:
                args.precluster = False
                logger.info(f"Alignment size: {n_seqs} sequences. Preclustering not needed.")
        except Exception as e:
            logger.warning(f"Could not check alignment size: {e}. Proceeding without auto-enable.")
            args.precluster = False
    elif args.precluster:
        logger.info("MMseqs2 preclustering explicitly enabled by user.")
    else:
        logger.info("MMseqs2 preclustering explicitly disabled by user.")

    # Load original alignment first to find reference sequence and process keep-sequences indices
    logger.info(f"Loading alignment from {in_path}...")
    headers_original, sequences_original = sca.readAlg(str(in_path))
    if not sequences_original:
        raise ValueError("No sequences read from alignment.")
    logger.info(f"Loaded alignment: N={len(headers_original)} sequences, L={len(sequences_original[0])} positions")
    
    # Process --keep-sequences and --keep-sequences-file: convert 1-based indices to sequence IDs
    keep_sequence_ids_list = None
    keep_indices_1based = []  # Store as 1-based from user input
    
    # Collect indices from command-line argument (1-based)
    if args.keep_sequences is not None:
        keep_indices_1based.extend(args.keep_sequences)
    
    # Collect indices from file (1-based)
    if args.keep_sequences_file is not None:
        file_path = Path(args.keep_sequences_file)
        if not file_path.exists():
            raise FileNotFoundError(f"Keep-sequences file not found: {file_path}")
        logger.info(f"Reading sequence indices to retain from file: {file_path} (1-based indices)")
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):  # Skip empty lines and comments
                    continue
                # Try to parse as space-separated or comma-separated numbers
                for part in line.replace(',', ' ').split():
                    try:
                        idx_1based = int(part)
                        keep_indices_1based.append(idx_1based)
                    except ValueError:
                        logger.warning(f"Skipping non-numeric value in keep-sequences file: {part}")
        logger.info(f"Read {len(keep_indices_1based)} sequence indices from file")
    
    # Convert all 1-based indices to sequence IDs (convert to 0-based for array access)
    if keep_indices_1based:
        keep_sequence_ids_list = []
        unique_indices_1based = sorted(set(keep_indices_1based))
        for idx_1based in unique_indices_1based:
            idx_0based = idx_1based - 1  # Convert to 0-based for array access
            if idx_0based < 0 or idx_0based >= len(headers_original):
                raise ValueError(f"Invalid sequence index {idx_1based} (1-based) for --keep-sequences "
                               f"(valid range: 1-{len(headers_original)})")
            seq_id = headers_original[idx_0based].split()[0]  # Get first part of header (before whitespace)
            keep_sequence_ids_list.append(seq_id)
        # Remove duplicates while preserving order
        seen = set()
        keep_sequence_ids_list = [x for x in keep_sequence_ids_list if not (x in seen or seen.add(x))]
        logger.info(f"Will retain {len(keep_sequence_ids_list)} user-specified sequences after preclustering "
                   f"(1-based indices: {unique_indices_1based})")
    
    # Clean sequences: strip whitespace and ensure consistent gap representation
    sequences_cleaned = [s.replace(" ", "").replace("\t", "").replace(".", "-").replace("\n", "").replace("\r", "") for s in sequences_original]
    
    # Normalize + remove non-standard sequences
    headers_original, sequences_original = filter_nonstandard(headers_original, sequences_cleaned)
    logger.info(f"After non-standard AA filter: N={len(headers_original)} sequences")
    
    if len(headers_original) == 0:
        # Diagnostic: check what characters are actually in the original sequences
        if sequences_cleaned and len(sequences_cleaned) > 0:
            sample_seq = sequences_cleaned[0]
            unique_chars = set("".join(sequences_cleaned[:min(10, len(sequences_cleaned))]))  # Check first 10 sequences
            unexpected = unique_chars - AA_ALLOWED
            raise ValueError(
                f"All sequences were filtered out! "
                f"Found unexpected characters: {sorted(unexpected)}. "
                f"Allowed characters: {sorted(AA_ALLOWED)}. "
                f"Sample sequence (first 100 chars): {sample_seq[:100]}"
            )
        else:
            raise ValueError("All sequences were filtered out! No sequences after reading alignment.")

    # Initial trim of highly gapped positions (for reference search)
    sequences_original, poskeep_init = sca.filterPos(sequences_original, [1], args.initial_trim_gap)
    logger.info(f"After initial trim: L={len(sequences_original[0])} positions (kept {len(poskeep_init)})")

    # Find reference sequence on TRIMMED alignment (before MMseqs2 preclustering)
    # This gives us i_ref and ref_header_id for retention during preclustering
    ref_header_id = None  # Will store the header ID to ensure it's retained during preclustering
    seq_pdb_stored = None
    ats_pdb_stored = None
    dist_pdb_stored = None
    i_ref_original = None  # Reference index in trimmed original alignment
    
    if args.i_ref is not None:
        # User provided reference index
        if args.i_ref >= len(headers_original):
            raise ValueError(f"Reference index {args.i_ref} is out of range (alignment has {len(headers_original)} sequences)")
        i_ref_original = args.i_ref
        ref_header_id = headers_original[i_ref_original].split()[0]
        logger.info(f"Using provided reference index: i_ref={i_ref_original} ({headers_original[i_ref_original]})")
    elif args.pdbid is not None:
        # Find reference sequence using PDB (on trimmed alignment)
        logger.info(f"Using PDB reference: {args.pdbid} chain {args.chainID}")
        seq_pdb, ats_pdb, dist_pdb = sca.pdbSeq(args.pdbid, args.chainID)
        seq_pdb_stored = seq_pdb
        ats_pdb_stored = ats_pdb
        dist_pdb_stored = dist_pdb
        
        # Search on TRIMMED alignment
        if args.species:
            try:
                logger.info(f"Finding reference sequence using species-based best match: {args.species}")
                i_ref_original = sca.MSAsearch(headers_original, sequences_original, seq_pdb, args.species)
            except Exception as e:
                logger.warning(f"Species-based MSAsearch failed ({e}); falling back to global MSAsearch.")
                i_ref_original = sca.MSAsearch(headers_original, sequences_original, seq_pdb)
        else:
            logger.info("Finding reference sequence using global MSAsearch")
            i_ref_original = sca.MSAsearch(headers_original, sequences_original, seq_pdb)
        
        logger.info(f"Reference index: i_ref={i_ref_original} ({headers_original[i_ref_original]})")
        ref_header_id = headers_original[i_ref_original].split()[0]
        
        # Log reference sequence header and pairwise alignment
        ref_header = headers_original[i_ref_original]
        ref_seq_from_aln = sequences_original[i_ref_original]
        ref_seq_ungapped = ref_seq_from_aln.replace("-", "").replace(".", "")
        seq_pdb_ungapped = seq_pdb.replace("-", "").replace(".", "")
        
        logger.info("=" * 80)
        logger.info("Reference Sequence Information:")
        logger.info(f"Selected reference sequence header: {ref_header}")
        logger.info("=" * 80)
        
        # Perform pairwise alignment using BioPython's PairwiseAligner
        if HAS_BIOALIGN:
            try:
                aligner = PairwiseAligner()
                aligner.mode = 'global'
                aligner.substitution_matrix = None
                aligner.match_score = 2.0
                aligner.mismatch_score = -1.0
                aligner.open_gap_score = -0.5
                aligner.extend_gap_score = -0.1
                
                alignments = aligner.align(seq_pdb_ungapped, ref_seq_ungapped)
                if alignments:
                    alignment = alignments[0]
                    aln_str = str(alignment)
                    lines = aln_str.split('\n')
                    formatted_lines = [line for line in lines if line.strip()]
                    
                    logger.info("=" * 80)
                    logger.info("Pairwise Alignment (PDB sequence vs Reference sequence):")
                    logger.info(f"PDB sequence ({args.pdbid} chain {args.chainID}): {len(seq_pdb_ungapped)} residues")
                    logger.info(f"Reference sequence ({ref_header}): {len(ref_seq_ungapped)} residues")
                    logger.info(f"Alignment score: {alignment.score:.2f}")
                    logger.info("-" * 80)
                    for line in formatted_lines:
                        logger.info(line)
                    logger.info("=" * 80)
            except Exception as e:
                logger.warning(f"Could not perform pairwise alignment with BioPython: {e}")
    elif args.refseq is not None:
        # Find reference sequence from file
        h_tmp, s_tmp = sca.readAlg(args.refseq)
        ref_seq = s_tmp[0].upper()
        i_ref_original = sca.MSAsearch(headers_original, sequences_original, ref_seq)
        ref_header_id = headers_original[i_ref_original].split()[0]
        logger.info(f"Reference index: i_ref={i_ref_original} ({headers_original[i_ref_original]})")
    else:
        # No reference specified - choose one
        i_ref_original = sca.chooseRefSeq(sequences_original)
        ref_header_id = headers_original[i_ref_original].split()[0]
        logger.info(f"No reference supplied; chose i_ref={i_ref_original} ({headers_original[i_ref_original]})")

    # Run preclustering on ORIGINAL untrimmed MSA (if requested)
    # Then re-execute all steps (trimming, filtering) on the preclustered MSA
    weights_by_id: Optional[Dict[str, float]] = None
    
    if args.precluster:
        logger.info("=" * 80)
        logger.info("MMseqs2 PRECLUSTERING ENABLED")
        logger.info("=" * 80)
        if ref_header_id:
            logger.info(f"Ensuring reference sequence '{ref_header_id}' is retained during preclustering")
            logger.info(f"Ensuring reference sequence '{ref_header_id}' is retained in preclustered alignment")
        
        # Precluster the ORIGINAL untrimmed alignment (from file)
        reps_fa, wts_tsv = run_precluster(
            str(in_path), work_prefix,  # Use original input file (untrimmed)
            cluster_id=args.cluster_id,
            coverage=args.cluster_coverage,
            cov_mode=args.cluster_cov_mode,
            keep_tmp=args.keep_mmseqs_tmp,
            logger=logger,
            keep_sequence_id=ref_header_id,  # Ensure reference is retained if we found it
            keep_sequence_ids=keep_sequence_ids_list  # User-specified sequences to retain
        )
        logger.info(f"Preclustered reps alignment: {reps_fa}")
        logger.info(f"Precluster weights: {wts_tsv}")
        weights_by_id = load_weights_tsv(wts_tsv)
        
        # Load the preclustered alignment
        logger.info(f"Loading preclustered alignment from {reps_fa}...")
        headers_full, sequences_full = sca.readAlg(reps_fa)
        if not sequences_full:
            raise ValueError("No sequences read from preclustered alignment.")
        
        # Re-execute all steps on preclustered alignment:
        # 1. Clean sequences
        sequences_cleaned = [s.replace(" ", "").replace("\t", "").replace(".", "-").replace("\n", "").replace("\r", "") for s in sequences_full]
        
        # 2. Filter non-standard AAs
        headers_full, sequences_full = filter_nonstandard(headers_full, sequences_cleaned)
        logger.info(f"After non-standard AA filter: N={len(headers_full)} sequences")
        
        if len(headers_full) == 0:
            raise ValueError("All sequences were filtered out from preclustered alignment!")
        
        # 3. Initial trim of highly gapped positions
        sequences_full, poskeep_init = sca.filterPos(sequences_full, [1], args.initial_trim_gap)
        logger.info(f"After initial trim: L={len(sequences_full[0])} positions (kept {len(poskeep_init)})")
        
        # 4. Find reference sequence in preclustered alignment by header ID (don't re-search)
        if ref_header_id:
            i_ref_new = None
            for idx, h in enumerate(headers_full):
                if h.split()[0] == ref_header_id:
                    i_ref_new = idx
                    break
            if i_ref_new is not None:
                logger.info(f"Reference sequence found in preclustered alignment at index: i_ref={i_ref_new} ({headers_full[i_ref_new]})")
                # Update i_ref_original to point to the preclustered alignment
                i_ref_original = i_ref_new
            else:
                raise ValueError(f"Reference sequence '{ref_header_id}' not found in preclustered alignment! "
                               f"This should not happen if preclustering retained it correctly.")
    else:
        # No preclustering - use original trimmed alignment
        headers_full = headers_original
        sequences_full = sequences_original

    logger.info(f"Final alignment: N={len(headers_full)} sequences, L={len(sequences_full[0])} positions")
    
    # Use the reference sequence we found earlier (i_ref_original)
    # It was found in the trimmed original alignment, and if preclustering was done,
    # we've already updated it to point to the preclustered alignment
    i_ref = i_ref_original if i_ref_original is not None else args.i_ref
    distmat = None
    ats = None
    sequences = None

    if i_ref is not None:
        if i_ref >= len(headers_full):
            raise ValueError(f"Reference index {i_ref} is out of range (alignment has {len(headers_full)} sequences)")
        
        if args.pdbid is not None:
            # Use PDB info for makeATS
            if seq_pdb_stored is None:
                seq_pdb_stored, ats_pdb_stored, dist_pdb_stored = sca.pdbSeq(args.pdbid, args.chainID)
            sequences, ats = sca.makeATS(sequences_full, ats_pdb_stored, seq_pdb_stored, i_ref, args.truncate)
        elif args.refseq is not None:
            # Use reference sequence file
            h_tmp, s_tmp = sca.readAlg(args.refseq)
            ref_seq = s_tmp[0].upper()
            ats_tmp = parse_refpos_file(args.refpos) if args.refpos else list(range(1, len(ref_seq) + 1))
            sequences, ats = sca.makeATS(sequences_full, ats_tmp, ref_seq, i_ref, args.truncate)
        else:
            # No PDB or refseq - use sequence from alignment
            ref_seq = sequences_full[i_ref]
            ats_tmp = parse_refpos_file(args.refpos) if args.refpos else list(range(len(sequences_full[0])))
            sequences, ats = sca.makeATS(sequences_full, ats_tmp, ref_seq, i_ref, args.truncate)
        
        logger.info(f"Using reference index: i_ref={i_ref} ({headers_full[i_ref]})")
        
        # Distance matrix remap for PDB case
        if args.pdbid is not None and dist_pdb_stored is not None:
            idx = {pos: j for j, pos in enumerate(ats_pdb_stored)}
            Lats = len(ats)
            dist_new = np.full((Lats, Lats), 1000.0, dtype=np.float32)  # Use float32 to save memory
            
            # Vectorized validity check using NumPy
            valid_mask = np.array([(p != "-" and p in idx) for p in ats], dtype=bool)
            valid_idx = np.where(valid_mask)[0]
            
            if len(valid_idx) > 0:
                # Vectorized index mapping
                pdb_idx = np.array([idx[ats[i]] for i in valid_idx], dtype=np.int32)
                # Use advanced indexing for efficient remapping
                dist_new[np.ix_(valid_idx, valid_idx)] = dist_pdb_stored[np.ix_(pdb_idx, pdb_idx)].astype(np.float32)
            
            distmat = dist_new
    else:
        # No reference was found - this shouldn't happen, but handle gracefully
        raise ValueError("No reference sequence found! This should not happen.")

    assert sequences is not None and ats is not None

    # Filter sequences
    # Important: i_ref must be valid for the sequences returned by makeATS
    # makeATS doesn't change the sequence order, so i_ref should still be valid
    logger.info(f"Filtering sequences (pre-filter): N={len(sequences)}, L={len(sequences[0]) if sequences else 0}")
    if len(sequences) == 0:
        raise ValueError("No sequences returned from makeATS! Check alignment and reference sequence.")
    if i_ref >= len(sequences):
        raise ValueError(f"Reference index {i_ref} is out of range for {len(sequences)} sequences after makeATS!")
    
    # Check reference sequence gap fraction before filtering
    ref_gap_frac = sequences[i_ref].count('-') / len(sequences[i_ref])
    logger.info(f"Reference sequence gap fraction: {ref_gap_frac:.3f} (max allowed: {args.parameters[1]})")
    if ref_gap_frac >= args.parameters[1]:
        logger.warning(f"Reference sequence has gap fraction {ref_gap_frac:.3f} >= max_fracgaps {args.parameters[1]}. "
                      f"This may cause all sequences to be filtered out!")
    
    alg0, seqw0, seqkeep = sca.filterSeq(
        sequences, i_ref,
        max_fracgaps=args.parameters[1],
        min_seqid=args.parameters[2],
        max_seqid=args.parameters[3],
    )
    
    # Ensure reference sequence is always kept (it may have been filtered out if its self-identity > max_seqid)
    # seqkeep contains the indices of kept sequences from the original alignment
    i_ref_filtered = None
    i_ref_was_moved = False  # Track if reference sequence was moved to index 0
    for idx, orig_idx in enumerate(seqkeep):
        if orig_idx == i_ref:
            i_ref_filtered = idx
            break
    
    if i_ref_filtered is None:
        # Reference sequence was filtered out - add it back
        logger.warning(f"Reference sequence (index {i_ref}) was filtered out by filterSeq. "
                      f"Adding it back to ensure it's always retained.")
        # Add reference sequence to the beginning of the filtered alignment
        alg0.insert(0, sequences[i_ref])
        headers_full_ref = headers_full[i_ref]
        headers = [headers_full_ref] + [headers_full[s] for s in seqkeep]
        # Recompute weights including the reference sequence
        seqw0 = sca.seqWeights(alg0, args.parameters[3])
        # Update seqkeep to include reference at position 0
        seqkeep.insert(0, i_ref)
        i_ref_filtered = 0
        i_ref_was_moved = True  # Mark that reference was moved to index 0
        logger.info(f"Reference sequence added back at index {i_ref_filtered} (original position: {i_ref})")
        logger.info(f"NOTE: Reference sequence was moved to index 0 (was filtered out, original position: {i_ref})")
    
    logger.info(f"After filterSeq: N={len(alg0)} (kept {len(seqkeep)})")
    if i_ref_was_moved:
        logger.info(f"Reference sequence is at index {i_ref_filtered} in filtered alignment (MOVED from original position {i_ref} to index 0)")
    else:
        logger.info(f"Reference sequence is at index {i_ref_filtered} in filtered alignment (original position: {i_ref}, not moved)")
    
    if len(alg0) == 0:
        raise ValueError(f"All sequences were filtered out! Check filtering parameters: "
                        f"max_fracgaps={args.parameters[1]}, min_seqid={args.parameters[2]}, max_seqid={args.parameters[3]}. "
                        f"Reference sequence gap fraction: {ref_gap_frac:.3f}")

    # Store i_ref_original: position in the original input MSA (before filtering)
    # i_ref_original is already set earlier (position in sequences_full/headers_full)
    # i_ref_filtered is the position in the final filtered alignment (alg0/alg1)
    # Note: position filtering (filterPos) doesn't change sequence order, so i_ref_filtered is valid for alg1

    # Filter positions
    alg1, iposkeep = sca.filterPos(alg0, seqw0, args.parameters[0])
    ats = [ats[i] for i in iposkeep]
    if distmat is not None:
        distmat = distmat[np.ix_(iposkeep, iposkeep)]
    logger.info(f"After filterPos: L={len(alg1[0]) if alg1 else 0} (kept {len(iposkeep)})")

    # Final weights
    if weights_by_id is not None:
        # Biopython record.id is token up to whitespace; use same convention here
        # Optimized: pre-split headers once if needed, but split() is already efficient
        seqw_raw = np.array([weights_by_id.get(h.split()[0], 1.0) for h in headers], dtype=np.float32)
        
        # Normalize cluster-size weights so sum(seqw) = M (number of clusters) = M_eff
        # Raw cluster sizes sum to total original sequences, but for SCA we need sum(seqw) = M_eff
        M_clusters = float(len(headers))
        sum_raw = float(seqw_raw.sum())
        if sum_raw > M_clusters * 1.5:  # Heuristic: if sum is much larger than M, these are raw cluster sizes
            # Normalize: multiply by M / sum_raw so sum(seqw) = M
            scale = M_clusters / sum_raw
            seqw = (seqw_raw * scale).astype(np.float32)
            effseqs = M_clusters  # M_eff = number of clusters after normalization
            logger.info(f"Using MMseqs2 cluster-size weights. "
                       f"Normalized from sum(raw)={sum_raw:.0f} to sum(seqw)={effseqs:.0f} (M_eff = number of clusters).")
        else:
            # Weights already appear normalized (or are not cluster sizes)
            seqw = seqw_raw
            effseqs = float(seqw.sum())
            logger.info(f"Using MMseqs2 cluster-size weights (sum(seqw)={effseqs:.2f}).")
    else:
        logger.info("Computing sequence weights via scaTools.seqWeights (may be expensive).")
        if len(alg1) > 50000:
            logger.info(f"Large alignment ({len(alg1)} sequences). This may take several minutes...")
        seqw = np.asarray(sca.seqWeights(alg1), dtype=np.float32)  # Use float32 to save memory
        effseqs = float(seqw.sum())

    # Final summary (always)
    Nseq_final = len(alg1)
    Npos_final = len(alg1[0]) if alg1 else 0
    logger.info("")  # Empty line for readability
    logger.info("=== Final alignment parameters ===")
    logger.info(f"Number of sequences: M = {Nseq_final}")
    logger.info(f"Number of effective sequences: M' = {effseqs:.2f}")
    logger.info(f"Number of alignment positions: L = {Npos_final}")
    if i_ref_original is not None:
        if i_ref_was_moved:
            logger.info(f"Reference sequence: original position {i_ref_original} in input MSA, MOVED to index 0 in filtered alignment")
        else:
            logger.info(f"Reference sequence: position {i_ref_original} in input MSA, at index {i_ref_filtered} in filtered alignment (not moved)")
    if distmat is not None:
        # Optimized: use list comprehension with enumerate (already efficient)
        structPos = [i for (i, k) in enumerate(ats) if k != "-"]
        logger.info(f"Number of positions in the ats: {len(ats)}")
        logger.info(f"Number of structure positions mapped: {len(structPos)}")
        logger.info(f"Size of the distance matrix: {distmat.shape[0]} x {distmat.shape[1]}")

    # Write processed alignment
    out_fa = Path("Inputs") / f"{base}_processed.fasta"
    write_fasta(out_fa, headers, alg1)
    logger.info(f"Wrote processed alignment FASTA: {out_fa}")

    # Build DB (gzip pickle)
    # Note: seqw0 might be a 2D array from filterSeq, ensure it's properly summed
    seqw0_sum = float(np.sum(seqw0)) if isinstance(seqw0, np.ndarray) else float(seqw0)
    
    # Store reference sequence header ID for finding it in sca-core
    ref_header_id = None
    if i_ref_original is not None and i_ref_filtered is not None and i_ref_filtered < len(headers):
        ref_header_id = headers[i_ref_filtered].split()[0]  # First token of header
    
    D = {
        "alg": alg1,
        "hd": headers,
        "seqw": seqw,
        "Nseq": Nseq_final,
        "Npos": Npos_final,
        "ats": ats,
        "effseqs": effseqs,
        "limitseqs": False,
        "effseqsPrelimit": int(seqw0_sum),
        "i_ref": int(i_ref_filtered) if i_ref_filtered is not None else None,  # Position in final processed alignment (alg)
        "i_ref_original": int(i_ref_original) if i_ref_original is not None else None,  # Original position in input MSA (alg_original)
        "i_ref_was_moved": bool(i_ref_was_moved),  # Flag indicating if reference was moved to index 0
        "ref_header_id": ref_header_id,  # Reference sequence header ID (for finding in sca-core)
        "trim_parameters": list(map(float, args.parameters)),
        "truncate_flag": bool(args.truncate),
        "preclustered": bool(weights_by_id is not None),
        "cluster_id": float(args.cluster_id) if weights_by_id is not None else None,
        "keep_sequences": keep_sequence_ids_list,  # Store sequence IDs to retain (for sca-core subsampling)
        # Store the original input MSA (that i_ref_original refers to)
        # This is the preclustered MSA if preclustering was done, otherwise the original input MSA
        "alg_original": sequences_full,  # Original input MSA sequences
        "hd_original": headers_full,  # Original input MSA headers
    }
    if distmat is not None:
        D["pdbid"] = args.pdbid
        D["pdb_chainID"] = args.chainID
        D["distmat"] = distmat
    if args.refseq is not None:
        D["refseq"] = args.refseq
    if args.refpos is not None:
        D["refpos"] = args.refpos

    if args.save_msa_num:
        logger.info("Converting alignment to numeric representation...")
        msa_num = sca.lett2num(alg1)
        npz_path = Path("Outputs") / f"{base}_msa_num.npz"
        logger.info(f"Saving numeric MSA to {npz_path}...")
        np.savez_compressed(npz_path, msa_num=msa_num)
        logger.info(f"Saved numeric MSA: {npz_path} (shape: {msa_num.shape})")

    db = {"sequence": D}

    if args.matlab:
        mat_path = Path("Outputs") / f"{base}.mat"
        # Convert ats to MATLAB cell array format (numpy object array of strings)
        # and remove None values (MATLAB can't handle None)
        mat_db = dict(db)  # Shallow copy
        if "sequence" in mat_db:
            mat_seq = dict(mat_db["sequence"])
            # Convert ats to MATLAB cell array format
            if "ats" in mat_seq:
                ats = mat_seq["ats"]
                ats_str = [str(a) for a in ats]
                mat_seq["ats"] = np.array(ats_str, dtype=object)
            # Ensure i_ref_was_moved is stored as a logical (boolean) for MATLAB
            if "i_ref_was_moved" in mat_seq:
                mat_seq["i_ref_was_moved"] = bool(mat_seq["i_ref_was_moved"])
            # Convert i_ref and i_ref_original to appropriate types for MATLAB if needed
            # (they should already be int/None, but ensure they're proper types)
            if "i_ref" in mat_seq and mat_seq["i_ref"] is not None:
                mat_seq["i_ref"] = int(mat_seq["i_ref"])
            if "i_ref_original" in mat_seq and mat_seq["i_ref_original"] is not None:
                mat_seq["i_ref_original"] = int(mat_seq["i_ref_original"])
            # Remove None values (MATLAB can't handle None)
            mat_seq_clean = {k: v for k, v in mat_seq.items() if v is not None}
            mat_db["sequence"] = mat_seq_clean
        # Wrap in a structure named after the output base name
        # Sanitize base name for MATLAB (must start with letter, only alphanumeric/underscore)
        import re
        mat_struct_name = re.sub(r'^[^a-zA-Z]', 'x', base)  # Ensure starts with letter
        mat_struct_name = re.sub(r'[^a-zA-Z0-9_]', '_', mat_struct_name)  # Replace invalid chars with underscore
        # This creates: <mat_struct_name>.sequence
        mat_struct = {mat_struct_name: mat_db}
        savemat(mat_path, mat_struct, appendmat=False, oned_as="column")
        logger.info(f"Wrote Matlab workspace: {mat_path} (structure name: {mat_struct_name})")

    out_db = Path("Outputs") / f"{base}.db.gz"
    import pickle
    with gzip.open(out_db, "wb") as f:
        pickle.dump(db, f, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info(f"Wrote database: {out_db}")
    logger.info("=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
