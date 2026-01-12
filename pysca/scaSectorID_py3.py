#!/usr/bin/env python3
"""pysca.scaSectorID_py3

Python 3 + logging rewrite of legacy scaSectorID.

This script performs sector identification on SCA results. Note that this
functionality is now also integrated into scaCore_py3.py (use --do-sector-id).

Expected DB contents
--------------------
- db['sequence'] from sca-process-msa
- db['sca'] from sca-core (must include Csca; may include random-trial eigen info)

Outputs
-------
Adds/overwrites db['sector'].
Canonical sector representation is saved as indices into `ats`:
  - db['sector']['ic_pos'] : list[list[int]] (0-based indices into `ats`)
  - db['sector']['ic_ats'] : list[list[str]] (ATS labels; works when ats are strings)

Matlab export
-------------
If -m is specified, writes a .mat file containing only the keys:
  sequence, sca, sector
This mirrors the original workflow and avoids savemat failures on Python-only objects.

Integration with scaCore:
-------------------------
The sector identification functionality in this script is fully integrated into
scaCore_py3.py. You can run both SCA core calculations and sector identification
in a single command:

  scaCore_py3.py input.db --do-sector-id

This script remains available for:
- Legacy workflows
- Running sector ID separately on existing SCA results
- When you need the MATLAB-specific sanitization (_mat_sanitize)
- When you prefer the in-place database update (this script) vs. separate output file (scaCore)
"""

from __future__ import annotations

import argparse
import gzip
import logging
import pickle
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from scipy.io import savemat
from scipy.stats import kurtosis, spearmanr

from pysca import scaTools as sca


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


def _sector_as_index_lists(sectors_raw: Any, Lpos: int):
    """Normalize output from sca.t(...) into list-of-lists of 0-based indices into ats.

    Supports:
      - list/tuple of index lists
      - tuple where first element is sectors (common in some versions)
      - dict mapping component -> indices or mask
      - 2D mask array (k x L) or (L x k), bool or numeric
      - 1D label vector of length L with sector IDs (0..k-1 or 1..k)
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
            raise ValueError("Bad sector mask shape %r for L=%r" % (arr.shape, Lpos))
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

    raise TypeError("Unrecognized sector format: type=%r, ndim=%r, shape=%r" % (type(sectors_raw), getattr(arr, "ndim", None), getattr(arr, "shape", None)))


def setup_logger(log_path: Optional[str], verbose: bool, quiet: bool) -> logging.Logger:
    logger = logging.getLogger("sca-sector")
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
        logger.info("Logging to file: %s" % fp)

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
    name = db_path.name
    if name.endswith(".db.gz"):
        return name[:-6]
    if name.endswith(".db"):
        return name[:-3]
    return db_path.stem


def main(argv: Optional[list[str]] = None) -> int:
    """
    Standalone independent component identification from SCA results.
    
    This script performs independent component identification on existing SCA core results.
    Note: This functionality is also integrated into scaCore_py3.py (use --do-sector-id).
    
    Workflow steps:
    1. Load database with SCA core results (must contain Csca)
    2. Retrieve or compute eigendecomposition (V, L from db['sca'])
    3. Select kpos eigenmodes (auto-based on Lrand or user-specified)
    4. Extract top kpos eigenvectors
    5. Perform ICA rotation (kica components, default kica=kpos)
    6. Identify significant positions using t-distribution fits
    7. Extract independent components and map to ATS numbering
    8. Compute Spearman correlation matrix between ICs
    9. Store results in db['sector']
    
    Use this script when:
    - Running independent component ID separately on existing SCA results
    - Preferring in-place database updates (this script) vs. separate output (scaCore)
    - Using legacy workflows
    
    Returns:
        int: Exit code (0 for success)
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("database", help="Database file produced by sca-core (.db/.db.gz)")
    ap.add_argument("-k", "--kpos", dest="kpos", type=int, default=0,
                    help="Number of significant eigenmodes (0 => choose automatically if possible)")
    ap.add_argument("--kica", type=int, default=None,
                    help="Number of independent components to compute via ICA (must be >= 2 and <= kpos). Default: kpos")
    ap.add_argument("-p", "--ic-cutoff", dest="cutoff", type=float, default=0.95,
                    help="Cutoff for selecting significant positions per IC (default 0.95)")
    ap.add_argument("--kmax-cap", type=int, default=10,
                    help="Safety cap on kpos if automatic selection returns larger (default 10)")
    ap.add_argument("--float32", action="store_true", default=False,
                    help="Store heavy arrays as float32")
    ap.add_argument("-m", "--matlab", action="store_true", dest="matfile", default=False)
    ap.add_argument("--output", dest="outputfile", default=None,
                    help="Output base name (basename only). Default: derived from input DB name.")
    ap.add_argument("--no-db", action="store_true", default=False,
                    help="Skip writing database file (.db.gz). Only write MATLAB file if --matlab is specified.")
    ap.add_argument("--log", default=None)
    ap.add_argument("--verbose", action="store_true", default=False)
    ap.add_argument("--quiet", action="store_true", default=False)

    args = ap.parse_args(argv)
    logger = setup_logger(args.log, args.verbose, args.quiet)

    db_path = Path(args.database)
    db = load_db(db_path)

    if "sequence" not in db or "sca" not in db:
        raise KeyError("Database must contain 'sequence' and 'sca'. Run sca-core first.")

    Csca = db["sca"].get("Csca", None)
    if Csca is None:
        raise KeyError("db['sca']['Csca'] not found.")
    Csca = np.asarray(Csca)

    logger.info("=== Independent Component identification ===")
    logger.info("Csca shape: %s" % (Csca.shape,))

    # Optional random trial eigenvalues information
    Lrand = db["sca"].get("Lrand", None)
    if Lrand is not None:
        Lrand = np.asarray(Lrand)
        logger.info("Random trials: %s" % (Lrand.shape,))

    logger.info("cutoff: %.3f" % args.cutoff)

    # Eigendecomposition (retrieve from db['sca'] if available, otherwise compute)
    Vfull = db["sca"].get("V")  # Full eigenvectors (renamed from Vfull)
    eigvals = db["sca"].get("L")  # Full eigenvalues (renamed from eigvals)
    if Vfull is None or eigvals is None:
        logger.info("V or L not found in db['sca'], computing eigendecomposition...")
        Vfull, eigvals = sca.eigenVect(Csca)
        # Store in sca field for future use
        db.setdefault("sca", {})
        db["sca"]["V"] = Vfull.astype(np.float32) if args.float32 else Vfull
        db["sca"]["L"] = eigvals.astype(np.float32) if args.float32 else eigvals
    else:
        logger.info("Using precomputed eigendecomposition from db['sca']")

    # Choose kpos
    kpos_user = int(args.kpos)
    kpos_auto = None
    kpos_was_auto = False
    
    # Always compute kpos_auto for reference, even if user specified kpos
    if Lrand is not None and Lrand.ndim == 2 and Lrand.shape[1] == eigvals.shape[0]:
        thr = np.quantile(Lrand, args.cutoff, axis=0)
        # Compute kpos_auto: count eigenvalues that exceed threshold
        count = np.sum(eigvals > thr)
        # Convert to Python int (handle both numpy scalar and Python int)
        if hasattr(count, 'item'):
            kpos_auto = int(count.item())
        else:
            kpos_auto = int(count)
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
    logger.info("kpos_auto (auto-estimated): %d" % kpos_auto)
    if kpos_user > 0:
        logger.info("kpos_user (user-specified): %d" % kpos_user)
    if kpos < kpos_before_cap:
        logger.info("kpos capped from %d to %d (--kmax-cap=%d)" % (kpos_before_cap, kpos, args.kmax_cap))
    logger.info("Using kpos=%d eigenmodes %s" % (kpos, "(auto-selected)" if kpos_was_auto else "(user-specified)"))

    V = Vfull[:, :kpos]
    top_eigvals = np.asarray(eigvals[:kpos])

    logger.info("Computed top eigenmodes: kpos=%d" % kpos)

    # Determine number of independent components to compute
    kica = args.kica
    if kica is None:
        kica = kpos  # Default: use all eigenmodes
    else:
        kica = int(kica)
        if kica < 2:
            raise ValueError("--kica must be >= 2 (got %d)" % kica)
        if kica > kpos:
            raise ValueError("--kica (%d) cannot exceed --kpos (%d)" % (kica, kpos))

    logger.info("Computing kica=%d independent components (from %d eigenmodes)" % (kica, kpos))

    # ICA rotation
    Vica, W = sca.rotICA(V, kmax=kica)
    logger.info("ICA rotation complete")

    # Compute excess kurtosis for each independent component
    logger.info("Computing excess kurtosis for each independent component...")
    ic_kurtosis = []
    for k in range(kica):
        kurt = kurtosis(Vica[:, k], fisher=True)  # Fisher=True gives excess kurtosis (kurtosis - 3)
        ic_kurtosis.append(float(kurt))
        logger.info("IC %d: excess kurtosis = %.4f" % (k+1, kurt))
        print("IC %d: excess kurtosis = %.4f" % (k+1, kurt))

    # Identify IC position sets (use kica, not kpos)
    out = sca.icList(Vica, kica, Csca, p_cut=args.cutoff)
    ic_list, icsize, sortedpos, cutoff, scaled_pdf, all_fits = out
    
    try:
        n_ic = len(ic_list)
    except Exception:
        n_ic = kica
    logger.info("Identified %d IC position sets" % n_ic)
    
    # Get ATS labels for mapping position indices to ATS numbering
    Lpos = int(Csca.shape[0])
    ats = db["sequence"].get("ats", list(range(Lpos)))
    ats_list = [str(a) for a in ats]
    
    # Log t-distribution information and significant positions
    logger.info("=== T-distribution fits and cutoffs ===")
    for k in range(kica):
        fit_params = all_fits[k]  # (df, loc, scale) tuple from t_dist.fit
        cutoff_val = cutoff[k]
        logger.info("IC %d: t-dist fit (df=%.2f, loc=%.4f, scale=%.4f), cutoff=%.4f, positions above cutoff: %d" %
                   (k+1, fit_params[0], fit_params[1], fit_params[2], cutoff_val, icsize[k]))
        
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
                logger.info("IC %d significant positions (ATS): %s" % (k+1, pos_str))

    # Independent component extraction
    logger.info("Extracting independent components...")
    sectors_raw = sca.t(Vica, ic_list)
    logger.debug("sca.t returned type=%r" % type(sectors_raw))

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
    ats = db["sequence"].get("ats", list(range(Lpos)))
    ats_list = [str(a) for a in ats]
    ic_ats = [[ats_list[i] for i in s] for s in ic_pos]

    logger.info("Identified %d independent components" % len(ic_pos))
    for i, sect in enumerate(ic_pos):
        logger.info("  Independent component %d: %d positions" % (i+1, len(sect)))

    # Store
    db.setdefault("sector", {})
    db["sector"]["kpos"] = kpos
    # Store kpos_auto as a scalar integer (not array) for MATLAB compatibility
    # Always store kpos_auto if it was computed (even if user specified kpos)
    logger.debug("Storing kpos_auto: was_auto=%s, kpos_auto=%s, type=%s" % (kpos_was_auto, kpos_auto, type(kpos_auto)))
    if kpos_auto is not None:
        # Ensure it's a Python int, not numpy scalar (which can become empty array in MATLAB)
        kpos_auto_stored = int(kpos_auto) if not isinstance(kpos_auto, int) else kpos_auto
        db["sector"]["kpos_auto"] = kpos_auto_stored
        logger.debug("  Stored kpos_auto as: %s (type: %s)" % (kpos_auto_stored, type(kpos_auto_stored)))
    else:
        db["sector"]["kpos_auto"] = None
        logger.debug("  Stored kpos_auto as: None (could not compute)")
    db["sector"]["kpos_was_auto"] = kpos_was_auto
    db["sector"]["kica"] = kica
    db["sector"]["ic_cutoff"] = float(args.cutoff)
    # V and L (full eigenvectors and eigenvalues) are stored in db['sca']
    # Only store top kpos results in sector
    db["sector"]["top_eigvals"] = top_eigvals.astype(np.float32) if args.float32 else top_eigvals  # Top kpos eigenvalues
    db["sector"]["V"] = V.astype(np.float32) if args.float32 else V  # Top kpos eigenvectors
    db["sector"]["Vica"] = Vica.astype(np.float32) if args.float32 else Vica
    db["sector"]["W"] = W.astype(np.float32) if (args.float32 and hasattr(W, "dtype")) else W
    
    # Compute Spearman rank correlation matrix between ICs
    logger.info("Computing Spearman rank correlations between ICs...")
    ic_corr = np.zeros((kica, kica))
    for i in range(kica):
        for j in range(kica):
            corr, _ = spearmanr(Vica[:, i], Vica[:, j])
            ic_corr[i, j] = corr
    db["sector"]["ic_corr_spearman"] = ic_corr.astype(np.float32 if args.float32 else np.float64)
    logger.info("IC correlation matrix shape: %s" % (ic_corr.shape,))
    
    # Format and log the correlation matrix
    logger.info("Spearman rank correlation matrix between ICs:")
    # Header row
    header = "IC" + "".join("%8d" % (i+1) for i in range(kica))
    logger.info(header)
    # Data rows
    for i in range(kica):
        row_str = "%2d" % (i+1) + "".join("%8.4f" % ic_corr[i, j] for j in range(kica))
        logger.info(row_str)
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

    # Persist updated DB (in-place)
    if not args.no_db:
        logger.info("Writing updated database: %s" % db_path)
        logger.info("This may take a moment for large databases...")
        save_db(db_path, db)
        # Get file size for user information
        db_size_mb = db_path.stat().st_size / (1024 * 1024)
        logger.info("Database written: %.1f MB" % db_size_mb)
    else:
        logger.info("Skipping database write (--no-db specified)")

    # Matlab export
    if args.matfile:
        outbase = derive_outbase(db_path, args.outputfile)
        mat_path = str(db_path.with_name(outbase + ".mat"))
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
        logger.info("Wrote Matlab workspace: %s (structure name: %s)" % (mat_path, mat_struct_name))

    logger.info("=== Done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

