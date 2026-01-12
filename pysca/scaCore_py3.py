#!/usr/bin/env python3
"""pysca.scaCore_py3

Memory-safer Python 3 rewrite of the legacy scaCore.py script with integrated
sector identification functionality.

Key features:
- Python 3 compatible, package-friendly imports.
- Structured logging (console + optional log file).
- Optional skipping of sequence-correlation calculations (simMat / Useq / Uica),
  since simMat is O(M^2) memory and will overflow for large MSAs.
- Intelligent automatic subsampling for seqSim (1.5 × effective sequences, capped).
- Uses provided sequence weights from the processed DB; if preclustered weights
  look like raw cluster sizes, they are normalized so sum(weights)=#clusters,
  matching classic "effective sequence" semantics.

Outputs (db['sca']):
  Di, Dia, Csca, tX, Proj, Vrand, Lrand, (Crand optional)
  and optionally simMat/Useq/Uica if enabled.

Integrated Sector Identification:
  When --do-sector-id is enabled, performs sector identification after SCA core
  calculations and stores results in db['sector']. This integrates the functionality
  of scaSectorID_py3.py into a single workflow step.
  
  Outputs (db['sector'] when --do-sector-id):
    kpos, ic_cutoff, eigvals, V, Vica, W, ic_list, pd, scaled_pd, ic_pos, ic_ats

The script reads .db or .db.gz created by sca-process-msa (scaProcessMSA_py3_big.py).

Integration with scaSectorID:
  This script now contains the full functionality of scaSectorID_py3.py when
  --do-sector-id is enabled. The standalone scaSectorID_py3.py script is still
  available for legacy workflows or when you need to run sector ID separately
  on existing SCA results.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.io import savemat
from scipy.stats import kurtosis

from pysca import scaTools as sca


def setup_logger(log_path: Optional[str], verbose: bool, quiet: bool) -> logging.Logger:
    logger = logging.getLogger("sca-core")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

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

    if log_path:
        fp = Path(log_path)
        fp.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(fp, mode="w", encoding="utf-8")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
        logger.info(f"Logging to file: {fp}")

    return logger


def load_db(path: Path) -> Dict[str, Any]:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            return pickle.load(f)
    with open(path, "rb") as f:
        return pickle.load(f)


def save_db(path: Path, db: Dict[str, Any]) -> None:
    if path.suffix == ".gz":
        with gzip.open(path, "wb") as f:
            pickle.dump(db, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(path, "wb") as f:
            pickle.dump(db, f, protocol=pickle.HIGHEST_PROTOCOL)


def _mat_sanitize(x: Any) -> Any:
    """Make nested structures savemat-friendly.
    
    - None -> empty float array
    - dict -> dict with sanitized values
    - list/tuple -> numeric ndarray if possible else object ndarray
    - numpy object/structured arrays -> elementwise sanitize
    - numpy scalars -> python scalars
    - unknown objects -> string
    """
    if x is None:
        return np.array([], dtype=float)
    
    if isinstance(x, (str, bytes, int, float, bool)):
        return x
    
    if isinstance(x, np.generic):
        return x.item()
    
    if isinstance(x, dict):
        return {k: _mat_sanitize(v) for k, v in x.items()}
    
    if isinstance(x, np.ndarray):
        if x.dtype.fields is not None or x.dtype == object:
            out = np.empty(x.shape, dtype=object)
            it = np.nditer(x, flags=["multi_index", "refs_ok", "zerosize_ok"], op_flags=["readonly"])
            for a in it:
                out[it.multi_index] = _mat_sanitize(a.item())
            return out
        return x
    
    if isinstance(x, (list, tuple)):
        # Special handling for list of lists (like ic_ats: list[list[str]])
        # MATLAB expects nested cell arrays, so we need to convert each inner list to object array
        if len(x) > 0 and isinstance(x[0], (list, tuple)):
            # This is a nested list (e.g., ic_ats)
            # Convert to array of object arrays for MATLAB cell array of cell arrays
            return np.array([np.array([str(item) for item in sublist], dtype=object) for sublist in x], dtype=object)
        try:
            arr = np.asarray(x)
            if arr.dtype == object:
                return _mat_sanitize(arr)
            return arr
        except Exception:
            return np.array([_mat_sanitize(v) for v in x], dtype=object)
    
    return str(x)


def prepare_matlab_db(db: Dict[str, Any]) -> Dict[str, Any]:
    """Prepare database for MATLAB export: convert ats to cell array of strings and sanitize structures."""
    mat_db = {}
    # Copy top-level structure
    for key, value in db.items():
        if key == "sequence" and isinstance(value, dict):
            # Deep copy sequence dict and convert ats
            mat_seq = dict(value)  # Shallow copy
            if "ats" in mat_seq:
                ats = mat_seq["ats"]
                # Convert to list of strings if not already
                ats_str = [str(a) for a in ats]
                # Create object array for MATLAB cell array
                mat_seq["ats"] = np.array(ats_str, dtype=object)
            mat_db[key] = mat_seq
        else:
            mat_db[key] = value
    
    # Sanitize the entire structure for MATLAB compatibility
    return _mat_sanitize(mat_db)


def derive_outbase(db_path: Path, out_arg: Optional[str]) -> str:
    if out_arg:
        return Path(out_arg).name
    # strip .db or .db.gz
    name = db_path.name
    if name.endswith(".db.gz"):
        return name[:-6]
    if name.endswith(".db"):
        return name[:-3]
    return db_path.stem


def _sector_as_index_lists(sectors_raw: Any, Lpos: int) -> list:
    """Normalize output from sca.t() into list-of-lists of 0-based indices.
    
    This is a simplified version of the helper function from scaSectorID_py3.py.
    Supports common formats returned by sca.t().
    """
    # tuple: first element is often the sectors object
    if isinstance(sectors_raw, tuple) and len(sectors_raw) > 0:
        sectors_raw = sectors_raw[0]
    
    # dict: keys are components, values are indices/masks/labels
    if isinstance(sectors_raw, dict):
        out = []
        for k in sorted(sectors_raw.keys()):
            one = _sector_as_index_lists(sectors_raw[k], Lpos)
            if len(one) == 1:
                out.append(one[0])
            else:
                out.extend(one)
        return out
    
    # list/tuple: indices
    if isinstance(sectors_raw, (list, tuple)):
        if len(sectors_raw) == 0:
            return []
        if all(isinstance(x, (list, tuple, np.ndarray)) for x in sectors_raw):
            out = []
            for x in sectors_raw:
                xarr = np.asarray(x)
                if xarr.ndim == 2:
                    out.extend(_sector_as_index_lists(xarr, Lpos))
                else:
                    xs = [int(i) for i in list(x)]
                    out.append(sorted({i for i in xs if 0 <= i < Lpos}))
            return out
        if isinstance(sectors_raw[0], (int, np.integer)):
            xs = [int(i) for i in sectors_raw]
            return [sorted({i for i in xs if 0 <= i < Lpos})]
    
    arr = np.asarray(sectors_raw)
    
    # 2D mask
    if arr.ndim == 2:
        if arr.shape[1] != Lpos and arr.shape[0] == Lpos:
            arr = arr.T
        if arr.shape[1] != Lpos:
            raise ValueError(f"Bad sector mask shape {arr.shape} for L={Lpos}")
        if arr.dtype != bool:
            arr = (arr != 0)
        return [np.flatnonzero(arr[k]).astype(int).tolist() for k in range(arr.shape[0])]
    
    # 1D labels (length Lpos)
    if arr.ndim == 1 and arr.shape[0] == Lpos:
        labels = arr.astype(int)
        uniq = sorted(set(labels.tolist()))
        # drop common "no-sector" labels
        uniq2 = [u for u in uniq if u not in (-1, 0)]
        out = []
        for u in uniq2:
            out.append(np.flatnonzero(labels == u).astype(int).tolist())
        # if 0 is actually used as a sector label
        if len(out) == 0 and 0 in uniq:
            out = [np.flatnonzero(labels == u).astype(int).tolist() for u in uniq]
        return out
    
    raise TypeError(f"Unrecognized sector format: type={type(sectors_raw)}, ndim={getattr(arr, 'ndim', None)}, shape={getattr(arr, 'shape', None)}")


def ensure_msa_num(seq: Dict[str, Any], logger: logging.Logger) -> np.ndarray:
    # Prefer explicit msa_num if present; else compute from alg
    msa_num = seq.get("msa_num", None)
    if msa_num is not None:
        return msa_num
    alg = seq.get("alg")
    if alg is None:
        raise KeyError("DB missing 'sequence.alg' and 'sequence.msa_num'.")
    logger.info("msa_num not found in DB; computing from alg (this may take a bit). ")
    return sca.lett2num(alg)


def normalize_cluster_weights(seqw: np.ndarray, logger: logging.Logger) -> Tuple[np.ndarray, float, float]:
    """If weights look like raw cluster sizes, normalize so sum=number_of_clusters.

    Returns: (seqw_norm, represented, scale)
    represented = sum(raw weights)
    scale = multiplier applied to raw weights
    """
    M = float(len(seqw))
    represented = float(np.sum(seqw))
    if M <= 0:
        return seqw, represented, 1.0

    # heuristic: if represented is much larger than M, assume raw cluster-size weights
    if represented > 2.0 * M:
        scale = M / represented
        logger.info(
            f"Detected cluster-size weights: sum(raw)={represented:.0f} >> M={M:.0f}. Normalizing so sum(weights)=M."
        )
        return (seqw * scale).astype(float), represented, scale
    return seqw.astype(float), represented, 1.0


def main(argv: Optional[list[str]] = None) -> int:
    """
    Main workflow for SCA core calculations and optional independent component identification.
    
    Workflow steps:
    1. Load database from sca-process-msa
    2. Normalize sequence weights (if from MMseqs2 preclustering)
    3. Compute positional weights (Di, Dia) using freq() and posWeights()
    4. Compute SCA matrix (Csca) and projectors (tX, Proj)
    5. Compute eigendecomposition of Csca (V, L) - stored in db['sca']
    6. Run randomization trials (Vrand, Lrand) for statistical comparison
    7. Optionally compute sequence correlations (simMat, Useq, Uica) with subsampling
    8. Optionally perform independent component identification:
       - Select kpos eigenmodes (auto or user-specified)
       - Perform ICA rotation
       - Identify significant positions using t-distribution fits
       - Extract independent components
    
    Key features:
    - Automatically enables sequence correlations if preclustered
    - Intelligent subsampling for seqSim (1.5 × M_eff, capped)
    - Integrates sector identification into single workflow
    - Supports float32 storage for memory efficiency
    
    Returns:
        int: Exit code (0 for success)
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("database", help="Database file from sca-process-msa (.db or .db.gz)")
    ap.add_argument("-n", dest="norm", default="frob",
                    help="Norm type for reducing SCA matrix: 'spec' or 'frob' (default: frob)")
    ap.add_argument("-t", "--Ntrials", dest="Ntrials", default=10, type=int,
                    help="Number of randomization trials (default 10)")
    ap.add_argument("-l", dest="lbda", default=0.03, type=float,
                    help="Regularization parameter (default 0.01)")
    ap.add_argument("--do-seqcorr", action="store_true", default=False,
                    help="Compute sequence correlation matrices (simMat/Useq/Uica). WARNING: simMat is O(M^2) memory.")
    ap.add_argument("--max-seqcorr-M", type=int, default=5000,
                    help="Safety limit: if M exceeds this and --do-seqcorr not set, skip seq-corr (default 5000)")
    ap.add_argument("--max-seqcorr-seqs", type=int, default=None,
                    help="Maximum sequences for seqSim subsampling (for large MSAs). If None, uses all sequences.")
    ap.add_argument("--seqcorr-ref", type=int, default=None,
                    help="Reference sequence index for seqSim subsampling (always retained).")
    ap.add_argument("--seqcorr-mmseqs2", action="store_true", default=False,
                    help="Use MMseqs2 clustering for seqSim subsampling (better diversity).")
    ap.add_argument("--seqcorr-keep-indices", type=int, nargs="+", default=None,
                   help="Additional sequence indices (1-based) to always retain in seqSim subsampling. "
                        "Can be combined with sequences stored in database from sca-process-msa.")
    ap.add_argument("--seqcorr-max-cap", type=int, default=50000,
                    help="Maximum sequences cap for seqSim subsampling (default 50000).")
    ap.add_argument("--no-auto-seqcorr-subsample", action="store_true", default=False,
                    help="Disable automatic seqSim subsampling (1.5 × M_eff). Uses all sequences unless --max-seqcorr-seqs is specified.")
    ap.add_argument("--store-crand", action="store_true", default=False,
                    help="Store randomized correlation matrices Crand in DB (can inflate DB size). Default: False.")
    ap.add_argument("--float32", action="store_true", default=False,
                    help="Store large matrices as float32 (saves memory). Default: False.")
    ap.add_argument("--do-sector-id", action="store_true", default=False,
                    help="Perform sector identification after SCA core calculations.")
    ap.add_argument("--kpos", type=int, default=0,
                    help="Number of significant eigenmodes for sector ID (0=auto, based on Lrand). Default: 0")
    ap.add_argument("--kica", type=int, default=None,
                    help="Number of independent components to compute via ICA (must be >= 2 and <= kpos). Default: kpos")
    ap.add_argument("--ic-cutoff", type=float, default=0.95,
                    help="Cutoff for selecting significant positions per IC (default 0.95).")
    ap.add_argument("--kmax-cap", type=int, default=10,
                    help="Safety cap on kpos if automatic selection returns larger (default 10).")
    ap.add_argument("-m", "--matlab", action="store_true", dest="matfile", default=False,
                    help="Write Matlab workspace .mat in Outputs/")
    ap.add_argument("--output", dest="outputfile", default=None,
                    help="Output base name (basename only). Default: derived from input DB name.")
    ap.add_argument("--no-db", action="store_true", default=False,
                    help="Skip writing database file (.db.gz). Only write MATLAB file if --matlab is specified.")
    ap.add_argument("--log", default=None, help="Write full debug log to this file.")
    ap.add_argument("--verbose", action="store_true", default=False)
    ap.add_argument("--quiet", action="store_true", default=False)

    args = ap.parse_args(argv)
    logger = setup_logger(args.log, args.verbose, args.quiet)

    db_path = Path(args.database)
    db = load_db(db_path)

    if "sequence" not in db:
        raise KeyError("Database missing 'sequence' key.")
    seq = db["sequence"]

    outbase = derive_outbase(db_path, args.outputfile)
    outdir = Path("Outputs")
    outdir.mkdir(parents=True, exist_ok=True)

    # Pull alignment and weights
    msa_num = ensure_msa_num(seq, logger)
    alg = seq.get("alg")
    hd = seq.get("hd", None)

    seqw_in = np.asarray(seq.get("seqw"), dtype=float)
    # seqw can be 1D (M,) or 2D (1, M) or (M, 1)
    # Flatten to 1D for validation
    seqw_flat = seqw_in.flatten()
    M_expected = msa_num.shape[0]
    if len(seqw_flat) != M_expected:
        raise ValueError(f"seqw shape does not match number of sequences: "
                        f"seqw has {len(seqw_flat)} elements (shape {seqw_in.shape}), "
                        f"but alignment has {M_expected} sequences")
    # Ensure seqw_in is 1D for consistency
    if seqw_in.ndim != 1:
        seqw_in = seqw_flat

    # Get reference sequence index from database (original position in input MSA)
    i_ref_original = seq.get("i_ref", None)
    ref_header_id = seq.get("ref_header_id", None)
    
    # Find reference sequence in current alignment by matching header ID
    # i_ref_original is the position in the original input MSA
    # We need to find where that sequence is in the current alignment
    i_ref_current = None
    if ref_header_id is not None and hd is not None:
        # Find reference sequence by matching header ID (first token)
        for idx, header in enumerate(hd):
            if header.split()[0] == ref_header_id:
                i_ref_current = idx
                break
        
        if i_ref_current is not None:
            logger.info(f"Reference sequence found at index {i_ref_current} in current alignment "
                       f"(original position in input MSA: {i_ref_original}, header ID: {ref_header_id})")
        else:
            logger.warning(f"Reference sequence (header ID: {ref_header_id}, original position: {i_ref_original}) "
                          f"not found in current alignment. Subsampling may not preserve reference sequence.")
    elif i_ref_original is not None:
        # Fallback: try using original position if it's still valid (no filtering happened)
        if 0 <= i_ref_original < M:
            i_ref_current = i_ref_original
            logger.info(f"Using original reference position {i_ref_original} (no header ID stored, assuming no filtering)")
        else:
            logger.warning(f"Reference sequence (original position {i_ref_original}) not found in current alignment. "
                          f"Subsampling may not preserve reference sequence.")

    # Normalize cluster weights if needed (Option 1 semantics)
    seqw, represented, scale = normalize_cluster_weights(seqw_in, logger)

    M, L = msa_num.shape
    logger.info("=== SCA core parameters ===")
    logger.info(f"M (sequences) = {M}")
    logger.info(f"L (positions) = {L}")
    logger.info(f"Sum weights (M_eff) = {seqw.sum():.2f}")
    if represented > 0 and scale != 1.0:
        logger.info(f"Represented sequences (sum raw cluster sizes) = {represented:.0f}")
    logger.info(f"Regularization = {args.lbda}")
    logger.info(f"Randomization trials = {args.Ntrials}")
    logger.info(f"Matrix dtype preference = {'float32' if args.float32 else 'float64'}")

    # Sequence correlations (optional)
    # Auto-enable if preclustering was used (alignment is already reduced to manageable size)
    # Since preclustering reduces alignment size and we use 1.5×M_eff subsampling,
    # sequence correlations are safe and useful to compute
    was_preclustered = seq.get("preclustered", False)
    if was_preclustered and not args.do_seqcorr:
        args.do_seqcorr = True
        logger.info("Preclustered alignment detected. Automatically enabling sequence correlations "
                   f"(alignment size {M} is manageable after preclustering, will subsample to ~1.5×M_eff).")
    
    D_sca: Dict[str, Any] = {}
    if args.do_seqcorr or (M <= args.max_seqcorr_M):
        logger.info("Computing sequence similarity matrix (simMat) and sequence projections...")
        
        # Determine reference sequence index for subsampling
        # Use user-specified index if provided, otherwise use the one from database
        if args.seqcorr_ref is not None:
            i_ref_seqcorr = args.seqcorr_ref - 1  # Convert from 1-based to 0-based
            if i_ref_seqcorr < 0 or i_ref_seqcorr >= M:
                raise ValueError(f"User-specified reference index {args.seqcorr_ref} (0-based: {i_ref_seqcorr}) "
                               f"is out of range (0-{M-1})")
            logger.info(f"Using user-specified reference sequence index {i_ref_seqcorr} (1-based: {args.seqcorr_ref})")
        elif i_ref_current is not None:
            i_ref_seqcorr = i_ref_current
            logger.info(f"Using reference sequence from database at index {i_ref_seqcorr} "
                       f"(original position in input MSA: {i_ref_original})")
        else:
            i_ref_seqcorr = None
            logger.warning("No reference sequence index available. Subsampling will not preserve a specific reference sequence.")
        
        if i_ref_seqcorr is not None:
            logger.info(f"Will retain reference sequence at index {i_ref_seqcorr} during subsampling")
        
        # Use automatic subsampling (1.5 × effective sequences) if enabled and seqw available
        # However, for seqProj computation, we need to be more conservative due to memory
        # seqProj is more memory-intensive than simMat (does SVD on large sparse matrices)
        auto_sub = not args.no_auto_seqcorr_subsample
        
        # Calculate suggested max_seqs, but cap more aggressively for memory safety
        # seqProj does SVD on large sparse matrices which is memory-intensive
        # Use a more conservative cap: 10,000 sequences for seqProj safety
        max_seqs_for_sim = args.max_seqcorr_seqs
        if max_seqs_for_sim is None and auto_sub:
            m_eff = seqw.sum()
            suggested_max = int(np.round(1.5 * m_eff))
            # Cap more aggressively: seqProj is memory-intensive
            seqproj_cap = 10000  # Conservative cap for seqProj memory safety
            if suggested_max > seqproj_cap:
                max_seqs_for_sim = seqproj_cap
                logger.info(f"Automatic subsampling suggests {suggested_max} sequences (1.5 × M_eff={m_eff:.0f}), "
                           f"but capping at {seqproj_cap} for memory safety (seqProj is memory-intensive).")
            elif suggested_max > M:
                max_seqs_for_sim = M  # Can't use more than we have
                logger.info(f"Automatic subsampling suggests {suggested_max} sequences, "
                           f"but only {M} available. Using all {M} sequences.")
                if M > seqproj_cap:
                    logger.warning(f"Warning: Using all {M} sequences for seqProj may cause memory issues. "
                                 f"Consider using --max-seqcorr-seqs {seqproj_cap} to limit memory usage.")
            else:
                max_seqs_for_sim = suggested_max
        
        # Collect keep_indices from multiple sources:
        # 1. Command-line argument (1-based, convert to 0-based)
        # 2. Database (if stored from sca-process-msa) - stored as sequence IDs, need to find indices
        keep_indices_combined = []
        
        # From command-line (1-based, convert to 0-based)
        if args.seqcorr_keep_indices is not None:
            keep_indices_combined.extend([idx - 1 for idx in args.seqcorr_keep_indices])
            logger.info(f"User specified {len(args.seqcorr_keep_indices)} sequences to retain (1-based indices: {args.seqcorr_keep_indices})")
        
        # From database (if available) - stored as sequence IDs from sca-process-msa
        keep_sequences_db = seq.get("keep_sequences", None)
        if keep_sequences_db is not None:
            # keep_sequences_db should be a list of sequence IDs (header first parts)
            # We need to find their indices in the current alignment
            if isinstance(keep_sequences_db, (list, tuple)):
                found_count = 0
                for seq_id in keep_sequences_db:
                    # Find this sequence ID in the current alignment
                    for idx, h in enumerate(hd):
                        if h.split()[0] == seq_id:
                            keep_indices_combined.append(idx)
                            found_count += 1
                            break
                if found_count > 0:
                    logger.info(f"Found {found_count} sequences from database (from sca-process-msa) to retain during subsampling")
        
        # Validate and filter keep_indices
        # If alignment was preclustered, indices from original alignment may be invalid
        keep_indices_validated = None
        if keep_indices_combined:
            keep_indices_validated = [idx for idx in keep_indices_combined if 0 <= idx < M]
            if len(keep_indices_validated) < len(keep_indices_combined):
                invalid = [idx for idx in keep_indices_combined if idx < 0 or idx >= M]
                logger.warning(f"Filtered out {len(invalid)} invalid keep_indices: {invalid} "
                             f"(valid range: 0-{M-1}). Using {len(keep_indices_validated)} valid indices.")
            if len(keep_indices_validated) == 0:
                keep_indices_validated = None
                logger.warning("All keep_indices were invalid. Proceeding without keep_indices.")
            else:
                # Remove duplicates
                keep_indices_validated = sorted(set(keep_indices_validated))
                logger.info(f"Will retain {len(keep_indices_validated)} additional sequences during subsampling "
                           f"(0-based indices: {keep_indices_validated}, 1-based: {[i+1 for i in keep_indices_validated]})")
        
        if auto_sub:
            logger.info(f"Computing sequence similarity matrix with automatic subsampling (M_eff={seqw.sum():.2f})...")
        else:
            logger.info(f"Computing sequence similarity matrix without automatic subsampling (M_eff={seqw.sum():.2f})...")
        simMat, selected_indices = sca.seqSim(
            msa_num,
            max_seqs=max_seqs_for_sim,  # Use calculated/capped value
            i_ref=i_ref_seqcorr,
            seqw=seqw,  # For automatic calculation and weighted selection
            use_mmseqs2=args.seqcorr_mmseqs2,
            max_seqs_cap=args.seqcorr_max_cap,  # Keep original cap for seqSim internal logic
            keep_indices=keep_indices_validated,  # Use validated indices
            auto_subsample=False if max_seqs_for_sim is not None else auto_sub,  # Disable auto if we set max_seqs
            return_indices=True
        )
        logger.info(f"Using {len(selected_indices)} sequences for simMat computation "
                   f"(from {M} total, M_eff={seqw.sum():.2f})")
        
        # Get subsampled alignment and weights for seqProj
        # seqProj needs msa_num and seqw, not simMat
        msa_num_subsampled = msa_num[selected_indices, :]
        # seqw is (1, M) or (M,), need to subsample it
        if seqw.ndim == 1:
            seqw_subsampled = seqw[selected_indices]
        else:
            seqw_subsampled = seqw[0, selected_indices] if seqw.shape[0] == 1 else seqw[selected_indices, 0]
            # Ensure it's 1D for seqProj (which expects it to be reshapable to (1, M'))
            seqw_subsampled = seqw_subsampled.flatten()
        
        # seqProj expects seqw as (1, M) shape
        if seqw_subsampled.ndim == 1:
            seqw_subsampled = seqw_subsampled.reshape(1, -1)
        
        Useq, Uica = sca.seqProj(msa_num_subsampled, seqw_subsampled)
        if args.float32:
            simMat = np.asarray(simMat, dtype=np.float32)
            Useq = np.asarray(Useq, dtype=np.float32)
            Uica = np.asarray(Uica, dtype=np.float32)
        D_sca["simMat"] = simMat
        D_sca["Useq"] = Useq
        D_sca["Uica"] = Uica
        logger.info(f"simMat shape: {getattr(simMat,'shape',None)}")
    else:
        logger.warning(
            f"Skipping sequence-correlation calculations because M={M} exceeds safety limit {args.max_seqcorr_M}. "
            f"Use --do-seqcorr to force (may require lots of RAM), or use --max-seqcorr-seqs for subsampling."
        )

    # scaTools (legacy) expects seqw as a 2D row vector (1, M)
    # seqw is already a 1D array from normalize_cluster_weights, reshape if needed
    if seqw.ndim == 1:
        seqw = seqw.reshape(1, -1)
    elif seqw.ndim == 2:
        # If it's 2D, ensure it's (1, M) not (M, 1)
        if seqw.shape[0] != 1 and seqw.shape[1] == 1:
            seqw = seqw.T
    # Ensure float dtype (should already be float from normalize_cluster_weights)
    seqw = np.asarray(seqw, dtype=float)

    # Positional weights and SCA matrix
    logger.info("Computing positional weights (Di, Dia)...")
    if M > 100000:
        logger.info("Large alignment detected; positional weights computation may take a few minutes...")
    pw = sca.posWeights(msa_num, seqw=seqw, lbda=args.lbda)

    # scaTools.posWeights() historically varies in return signature across versions:
    #  - some return (Di, Dia)
    #  - yours returns (Wia, Dia, Di)
    if isinstance(pw, tuple) and len(pw) == 3:
        Wia, Dia, Di = pw
    elif isinstance(pw, tuple) and len(pw) == 2:
        Di, Dia = pw
        Wia = None
    else:
        raise TypeError(f"Unexpected return from sca.posWeights: {type(pw)}")

    logger.info("Computing SCA matrix (Csca), projected alignment (tX), and projector (Proj)...")
    if L > 1000:
        logger.info(f"Large alignment ({L} positions); SCA matrix computation may take several minutes...")
    # scaTools.scaMat signature in this codebase:
    # scaMat(alg, seqw=1, norm="frob", lbda=0, freq0=...)
    Csca, tX, Proj = sca.scaMat(msa_num, seqw=seqw, norm=args.norm, lbda=args.lbda)
    
    # Compute eigendecomposition of Csca (always, not just for sector ID)
    logger.info("Computing eigendecomposition of Csca...")
    Vfull, eigvals = sca.eigenVect(Csca)

    if args.float32:
        Di = np.asarray(Di, dtype=np.float32)
        Dia = np.asarray(Dia, dtype=np.float32)
        Csca = np.asarray(Csca, dtype=np.float32)
        tX = np.asarray(tX, dtype=np.float32)
        Proj = np.asarray(Proj, dtype=np.float32)
        Vfull = np.asarray(Vfull, dtype=np.float32)
        eigvals = np.asarray(eigvals, dtype=np.float32)

    # Randomization trials
    logger.info(f"Computing randomized trials (eigenvalues/vectors vs randomized alignments, Ntrials={args.Ntrials})...")
    if args.Ntrials > 10:
        logger.info(f"Note: This may take several minutes for {args.Ntrials} trials on large alignments.")
    
    import time
    start_time = time.time()
    Vrand, Lrand, Crand = sca.randomize(msa_num, args.Ntrials, seqw, norm=args.norm, lbda=args.lbda)
    elapsed = time.time() - start_time
    logger.info(f"Randomization complete: {args.Ntrials} trials in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    if args.float32:
        Vrand = np.asarray(Vrand, dtype=np.float32)
        Lrand = np.asarray(Lrand, dtype=np.float32)
        if args.store_crand:
            Crand = np.asarray(Crand, dtype=np.float32)

    # Populate sca dict
    D_sca.update({
        "lbda": float(args.lbda),
        "Ntrials": int(args.Ntrials),
        "Di": Di,
        "Dia": Dia,
        "Csca": Csca,
        "tX": tX,
        "Proj": Proj,
        "Vrand": Vrand,
        "Lrand": Lrand,
        "V": Vfull,  # Full eigenvectors of Csca (renamed from Vfull)
        "L": eigvals,  # Full eigenvalues of Csca (renamed from eigvals)
    })
    if args.store_crand:
        D_sca["Crand"] = Crand

    # Carry weight metadata forward for clarity
    D_sca["seqw_sum_raw"] = float(represented)
    D_sca["seqw_scale"] = float(scale)

    db["sca"] = D_sca

    # Sector identification (optional)
    if args.do_sector_id:
        logger.info("=== Independent Component identification ===")
        logger.info(f"Csca shape: {Csca.shape}")
        logger.info(f"cutoff: {args.ic_cutoff:.3f}")
        
        # Eigendecomposition (already computed above and stored in db['sca'])
        # Retrieve V and L from sca field (renamed from Vfull and eigvals)
        Vfull = db["sca"].get("V")
        eigvals = db["sca"].get("L")
        if Vfull is None or eigvals is None:
            logger.warning("V or L not found in db['sca'], recomputing...")
            Vfull, eigvals = sca.eigenVect(Csca)
        else:
            logger.info("Using precomputed eigendecomposition from db['sca']")
        
        # Choose kpos
        kpos_user = int(args.kpos)
        kpos_auto = None
        kpos_was_auto = False
        
        # Always compute kpos_auto for reference, even if user specified kpos
        # Use chooseKpos which sets threshold to mean of second eigenvalue of Lrand
        if Lrand is not None and Lrand.ndim == 2 and Lrand.shape[1] == eigvals.shape[0]:
            kpos_auto = sca.chooseKpos(eigvals, Lrand)
            if kpos_auto <= 0:
                kpos_auto = 1
        else:
            # Fallback: pick a small default
            kpos_auto = int(min(6, eigvals.shape[0]))
        
        # Now determine which kpos to use
        if kpos_user <= 0:
            # Use auto-estimated value
            kpos = kpos_auto
            kpos_was_auto = True
        else:
            # Use user-specified value
            kpos = kpos_user
            kpos_was_auto = False
        
        # Apply safety cap
        kpos_before_cap = kpos
        kpos = min(kpos, int(args.kmax_cap), Vfull.shape[1])
        
        # Clean log output
        logger.info(f"kpos_auto (auto-estimated): {kpos_auto}")
        if kpos_user > 0:
            logger.info(f"kpos_user (user-specified): {kpos_user}")
        if kpos < kpos_before_cap:
            logger.info(f"kpos capped from {kpos_before_cap} to {kpos} (--kmax-cap={args.kmax_cap})")
        logger.info(f"Using kpos={kpos} eigenmodes {'(auto-selected)' if kpos_was_auto else '(user-specified)'}")
        
        V = Vfull[:, :kpos]
        top_eigvals = np.asarray(eigvals[:kpos])
        logger.info(f"Computed top eigenmodes: kpos={kpos}")
        
        # Determine number of independent components to compute
        kica = args.kica
        if kica is None:
            kica = kpos  # Default: use all eigenmodes
        else:
            kica = int(kica)
            if kica < 2:
                raise ValueError(f"--kica must be >= 2 (got {kica})")
            if kica > kpos:
                raise ValueError(f"--kica ({kica}) cannot exceed --kpos ({kpos})")
        
        logger.info(f"Computing kica={kica} independent components (from {kpos} eigenmodes)")
        
        # ICA rotation
        logger.info("Performing ICA rotation...")
        Vica, W = sca.rotICA(V, kmax=kica)
        logger.info("ICA rotation complete")
        
        # Compute excess kurtosis for each independent component
        logger.info("Computing excess kurtosis for each independent component...")
        ic_kurtosis = []
        for k in range(kica):
            kurt = kurtosis(Vica[:, k], fisher=True)  # Fisher=True gives excess kurtosis (kurtosis - 3)
            ic_kurtosis.append(float(kurt))
            logger.info(f"IC {k+1}: excess kurtosis = {kurt:.4f}")
            print(f"IC {k+1}: excess kurtosis = {kurt:.4f}")
        
        # Identify IC position sets (use kica, not kpos)
        out = sca.icList(Vica, kica, Csca, p_cut=args.ic_cutoff)
        ic_list, icsize, sortedpos, cutoff, scaled_pdf, all_fits = out
        
        # Get ATS labels for mapping position indices to ATS numbering
        Lpos = int(Csca.shape[0])
        ats = seq.get("ats", list(range(Lpos)))
        ats_list = [str(a) for a in ats]
        
        # Log t-distribution information and significant positions
        logger.info("=== T-distribution fits and cutoffs ===")
        for k in range(kica):
            fit_params = all_fits[k]  # (df, loc, scale) tuple from t_dist.fit
            cutoff_val = cutoff[k]
            logger.info(f"IC {k+1}: t-dist fit (df={fit_params[0]:.2f}, loc={fit_params[1]:.4f}, scale={fit_params[2]:.4f}), "
                       f"cutoff={cutoff_val:.4f}, positions above cutoff: {icsize[k]}")
            
            # Extract significant positions for this IC (in ATS numbering)
            if k < len(ic_list):
                ic_unit = ic_list[k]
                if hasattr(ic_unit, 'items'):
                    # Get position indices from Unit object (items is a set or list)
                    pos_indices = list(ic_unit.items)
                elif isinstance(ic_unit, (list, tuple)):
                    # Fallback: if ic_list is already a list of position indices
                    pos_indices = list(ic_unit)
                else:
                    pos_indices = []
                
                # Sort positions by their value along this IC (descending order)
                if len(pos_indices) > 0:
                    # Get values for these positions along IC k
                    pos_values = Vica[pos_indices, k]
                    # Sort by value (descending: highest values first)
                    sorted_pairs = sorted(zip(pos_indices, pos_values), key=lambda x: x[1], reverse=True)
                    pos_indices_sorted = [idx for idx, val in sorted_pairs]
                    
                    # Map to ATS labels
                    pos_ats = [ats_list[i] for i in pos_indices_sorted]
                    # Format as sprintf('%g+', positions) - positions separated by '+'
                    pos_str = '+'.join(str(p) for p in pos_ats) + '+'
                    logger.info(f"IC {k+1} significant positions (ATS): {pos_str}")
        
        # Independent component extraction
        logger.info("Extracting independent components...")
        sectors_raw = sca.t(Vica, ic_list)
        logger.debug(f"sca.t returned type={type(sectors_raw)}")
        
        # Convert to index lists
        Lpos = int(Csca.shape[0])
        ic_pos = _sector_as_index_lists(sectors_raw, Lpos)
        
        # Sort positions in each independent component by their value along the corresponding IC
        # (independent components correspond to ICs, so ic_pos[k] corresponds to IC k)
        ic_pos_sorted = []
        for k, sector in enumerate(ic_pos):
            if len(sector) > 0 and k < kica:
                # Get values for positions in this independent component along IC k
                pos_values = Vica[sector, k]
                # Sort by value (descending: highest values first)
                sorted_pairs = sorted(zip(sector, pos_values), key=lambda x: x[1], reverse=True)
                sector_sorted = [idx for idx, val in sorted_pairs]
                ic_pos_sorted.append(sector_sorted)
            else:
                # If no IC info available, keep original order
                ic_pos_sorted.append(list(sector))
        
        ic_pos = ic_pos_sorted
        
        # Map to ATS labels
        ats = seq.get("ats", list(range(Lpos)))
        ats_list = [str(a) for a in ats]
        ic_ats = [[ats_list[i] for i in s] for s in ic_pos]
        
        # Store sector results
        db.setdefault("sector", {})
        db["sector"]["kpos"] = kpos
        # Store kpos_auto as a scalar integer (not array) for MATLAB compatibility
        # Always store kpos_auto if it was computed (even if user specified kpos)
        logger.debug(f"Storing kpos_auto: was_auto={kpos_was_auto}, kpos_auto={kpos_auto}, type={type(kpos_auto)}")
        if kpos_auto is not None:
            # Ensure it's a Python int, not numpy scalar (which can become empty array in MATLAB)
            kpos_auto_stored = int(kpos_auto) if not isinstance(kpos_auto, int) else kpos_auto
            db["sector"]["kpos_auto"] = kpos_auto_stored
            logger.debug(f"  Stored kpos_auto as: {kpos_auto_stored} (type: {type(kpos_auto_stored)})")
        else:
            db["sector"]["kpos_auto"] = None
            logger.debug(f"  Stored kpos_auto as: None (could not compute)")
        db["sector"]["kpos_was_auto"] = kpos_was_auto
        db["sector"]["kica"] = kica
        db["sector"]["ic_cutoff"] = float(args.ic_cutoff)
        # V and L (full eigenvectors and eigenvalues) are stored in db['sca']
        # Only store top kpos results in sector
        db["sector"]["top_eigvals"] = top_eigvals.astype(np.float32) if args.float32 else top_eigvals  # Top kpos eigenvalues
        db["sector"]["V"] = V.astype(np.float32) if args.float32 else V  # Top kpos eigenvectors
        db["sector"]["Vica"] = Vica.astype(np.float32) if args.float32 else Vica
        if hasattr(W, "dtype"):
            db["sector"]["W"] = W.astype(np.float32) if args.float32 else W
        else:
            db["sector"]["W"] = W
        
        # Compute IC correlation analysis: sum of triu(Csca[ic_pos, ic_pos], 1) / n_pos for each IC
        logger.info("Computing IC correlation analysis (sum of triu(Csca[ic_pos, ic_pos], 1) / n_pos)...")
        ic_corr_mean = []
        for k in range(kica):
            # Get positions for this IC
            if k < len(ic_list):
                ic_unit = ic_list[k]
                if hasattr(ic_unit, 'items'):
                    ic_positions = list(ic_unit.items)
                elif isinstance(ic_unit, (list, tuple)):
                    ic_positions = list(ic_unit)
                else:
                    ic_positions = []
            else:
                ic_positions = []
            
            if len(ic_positions) > 0:
                # Extract submatrix for this IC
                ic_submatrix = Csca[np.ix_(ic_positions, ic_positions)]
                # Get upper triangle (excluding diagonal, k=1 means start from first off-diagonal)
                triu_values = np.triu(ic_submatrix, k=1)
                # Sum all values in upper triangle
                triu_sum = np.sum(triu_values)
                # Divide by number of positions in IC
                ic_corr_mean_val = triu_sum / len(ic_positions)
            else:
                ic_corr_mean_val = 0.0
            
            ic_corr_mean.append(float(ic_corr_mean_val))
            logger.info(f"IC {k+1}: mean correlation = {ic_corr_mean_val:.6f} (n_pos = {len(ic_positions)})")
            print(f"IC {k+1}: mean correlation = {ic_corr_mean_val:.6f} (n_pos = {len(ic_positions)})")
        
        db["sector"]["ic_corr_mean"] = ic_corr_mean
        db["sector"]["ic_list"] = ic_list
        db["sector"]["icsize"] = icsize
        db["sector"]["sortedpos"] = sortedpos
        # T-distribution information - structured format: list of dicts, one per IC
        t_dist_struct = []
        for k in range(kica):
            fit_params = all_fits[k]  # (df, loc, scale) tuple
            t_dist_struct.append({
                "df": float(fit_params[0]),
                "loc": float(fit_params[1]),
                "scale": float(fit_params[2]),
                "cutoff": float(cutoff[k]),
                "scaled_pdf": np.asarray(scaled_pdf[k], dtype=np.float32 if args.float32 else np.float64)
            })
        db["sector"]["t_dist"] = t_dist_struct
        db["sector"]["ic_pos"] = ic_pos
        db["sector"]["ic_ats"] = ic_ats
        db["sector"]["ic_kurtosis"] = ic_kurtosis
        logger.info("Independent Component identification complete")

    # Write outputs
    if not args.no_db:
        out_db = outdir / f"{outbase}.db.gz"
        logger.info(f"Writing updated database: {out_db}")
        logger.info("This may take a moment for large databases...")
        save_db(out_db, db)
        # Get file size for user information
        db_size_mb = out_db.stat().st_size / (1024 * 1024)
        logger.info(f"Database written: {db_size_mb:.1f} MB")
    else:
        logger.info("Skipping database write (--no-db specified)")

    if args.matfile:
        mat_path = outdir / f"{outbase}.mat"
        # Prepare MATLAB-friendly format: convert ats to cell array of strings
        mat_db = prepare_matlab_db(db)
        # Wrap in a structure named after the output base name
        # Sanitize base name for MATLAB (must start with letter, only alphanumeric/underscore)
        import re
        mat_struct_name = re.sub(r'^[^a-zA-Z]', 'x', outbase)  # Ensure starts with letter
        mat_struct_name = re.sub(r'[^a-zA-Z0-9_]', '_', mat_struct_name)  # Replace invalid chars with underscore
        # This creates: <mat_struct_name>.sequence, <mat_struct_name>.sca, <mat_struct_name>.sector
        mat_struct = {mat_struct_name: mat_db}
        savemat(mat_path, mat_struct, appendmat=False, oned_as="column")
        logger.info(f"Wrote Matlab workspace: {mat_path} (structure name: {mat_struct_name})")

    logger.info("=== Done scaCore ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
