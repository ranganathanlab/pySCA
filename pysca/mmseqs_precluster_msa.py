#!/usr/bin/env python3
"""
mmseqs_precluster_msa.py

Precluster a multiple sequence alignment (FASTA/Stockholm/Clustal) with MMseqs2 at a given identity threshold,
then emit a reduced alignment containing only cluster representatives PLUS a weights file giving cluster sizes.

Why:
  - For very large MSAs (e.g., 330k x 680), Python-object overhead and any O(N^2) weighting becomes infeasible.
  - MMseqs2 clustering reduces redundancy cheaply; cluster sizes provide principled sequence weights.

Inputs:
  alignment file (FASTA/Stockholm/Clustal). Can be gzipped if Biopython can read it via open() wrapper.

Outputs:
  <out_prefix>_reps.fasta        aligned representatives (still aligned; taken directly from input alignment)
  <out_prefix>_weights.tsv       tab-separated: rep_id <TAB> cluster_size

Requires:
  - mmseqs in PATH
  - biopython

Notes:
  - Clustering is performed on *ungapped* sequences (alignment gaps removed).
  - Representative IDs are taken from the alignment record.id (Biopython splits at whitespace).
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord


def detect_format(path: str) -> str:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                return "fasta"
            if line.upper().startswith("CLUSTAL"):
                return "clustal"
            if line.startswith("#") and "STOCKHOLM" in line.upper():
                return "stockholm"
            break
    # fall back (AlignIO can often infer poorly); default to fasta
    return "fasta"


def read_alignment(path: str):
    fmt = detect_format(path)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        aln = AlignIO.read(handle, fmt)
    return aln


def write_ungapped_fasta(aln, out_fa: Path) -> None:
    with out_fa.open("w", encoding="utf-8") as f:
        for rec in aln:
            sid = rec.id
            seq = str(rec.seq).replace(".", "-")  # stockholm sometimes has '.'
            ungapped = seq.replace("-", "")
            if not ungapped:
                continue
            f.write(f">{sid}\n{ungapped}\n")


def run_mmseqs_cluster(ungapped_fa: Path, tmpdir: Path, min_seq_id: float, coverage: float, cov_mode: int) -> Path:
    """
    Returns path to tsv with columns: rep  member
    """
    ensure = lambda p: p.mkdir(parents=True, exist_ok=True)
    ensure(tmpdir)

    db = tmpdir / "db"
    clu = tmpdir / "clu"
    repdb = tmpdir / "repdb"
    tsv = tmpdir / "clu.tsv"

    # mmseqs createdb input.fasta db
    # Suppress verbose output - MMseqs2 is very chatty
    subprocess.run(
        ["mmseqs", "createdb", str(ungapped_fa), str(db)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    # clustering
    # --min-seq-id sets identity threshold
    # -c sets coverage threshold, --cov-mode sets coverage mode (0=bidirectional, 1=query, 2=target)
    subprocess.run(
        [
            "mmseqs", "cluster",
            str(db), str(clu), str(tmpdir / "tmp"),
            "--min-seq-id", str(min_seq_id),
            "-c", str(coverage),
            "--cov-mode", str(cov_mode),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )

    # emit rep->member TSV
    subprocess.run(
        ["mmseqs", "createtsv", str(db), str(db), str(clu), str(tsv)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    return tsv


def compute_cluster_sizes(tsv_path: Path) -> Dict[str, int]:
    """
    Compute cluster sizes from MMseqs2 TSV output.
    
    The TSV file contains lines: rep_id <TAB> member_id
    This function counts the number of members per representative.
    
    Args:
        tsv_path: Path to MMseqs2 cluster TSV file (rep <TAB> member format)
    
    Returns:
        Dictionary mapping representative ID to cluster size (number of members)
    """
    # TSV lines: rep \t member
    sizes: Dict[str, int] = {}
    with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rep, member = line.split("\t")[:2]
            sizes[rep] = sizes.get(rep, 0) + 1
    return sizes


def write_reps_alignment(aln, rep_sizes: Dict[str, int], out_fa: Path, out_weights: Path, keep_sequence_ids: Optional[List[str]] = None) -> None:
    """
    Write representative alignment and weights file.
    
    This function writes the cluster representatives (with aligned sequences from input)
    to a FASTA file and creates a TSV file with representative IDs and their cluster sizes.
    
    Args:
        aln: BioPython Alignment object (source of aligned sequences)
        rep_sizes: Dictionary mapping representative ID to cluster size
        out_fa: Output FASTA file path for representatives
        out_weights: Output TSV file path for weights (rep_id <TAB> cluster_size)
        keep_sequence_ids: Optional list of sequence IDs to ensure are included as representatives
                          (useful for retaining reference sequences)
    
    Notes:
        - Aligned sequences are taken directly from the input alignment (gaps preserved)
        - If keep_sequence_ids are provided, they are added to representatives if not already present
        - Cluster sizes are used as sequence weights in downstream SCA analysis
    """
    rep_set = set(rep_sizes.keys())
    
    # If keep_sequence_ids is specified, ensure all are in the representatives
    if keep_sequence_ids:
        for keep_sequence_id in keep_sequence_ids:
            # Check if the sequence exists in the alignment
            found_in_aln = False
            for rec in aln:
                if rec.id.split()[0] == keep_sequence_id:  # Match first part of ID
                    found_in_aln = True
                    if keep_sequence_id not in rep_set:
                        # Add it as a representative with cluster size 1 (or merge with existing cluster)
                        rep_set.add(keep_sequence_id)
                        if keep_sequence_id not in rep_sizes:
                            rep_sizes[keep_sequence_id] = 1
                        print(f"Ensuring sequence '{keep_sequence_id}' is retained as representative")
                    break
            if not found_in_aln:
                print(f"Warning: Sequence ID '{keep_sequence_id}' not found in alignment. Cannot ensure retention.")
    
    # build dict id->aligned seq
    reps: List[Tuple[str, str]] = []
    for rec in aln:
        rec_id = rec.id.split()[0]  # Get first part of ID (before whitespace)
        if rec_id in rep_set:
            seq = str(rec.seq).replace(".", "-").upper()
            reps.append((rec_id, seq))

    # Some reps may be missing if IDs don't match; warn via count
    if len(reps) != len(rep_set):
        missing = rep_set - {r[0] for r in reps}
        print(f"Warning: {len(missing)} representatives not found in alignment by ID. Example: {next(iter(missing)) if missing else 'n/a'}")

    with out_fa.open("w", encoding="utf-8") as f:
        for rid, seq in reps:
            f.write(f">{rid}\n{seq}\n")

    with out_weights.open("w", encoding="utf-8") as f:
        for rid, _seq in reps:
            f.write(f"{rid}\t{rep_sizes.get(rid, 1)}\n")


def main():
    """
    Main entry point for MMseqs2 MSA preclustering.
    
    This script uses MMseqs2 to cluster sequences at a given identity threshold,
    then outputs a reduced alignment containing only cluster representatives
    and a weights file with cluster sizes.
    
    Workflow:
    1. Detect alignment format (FASTA/Stockholm/Clustal)
    2. Read alignment using BioPython
    3. Write ungapped sequences to temporary FASTA (MMseqs2 requires ungapped)
    4. Run MMseqs2 clustering (createdb, cluster, createtsv)
    5. Compute cluster sizes from TSV output
    6. Write aligned representatives FASTA and weights TSV
    
    Key features:
    - Clustering is performed on ungapped sequences
    - Aligned sequences in output are taken from input alignment (gaps preserved)
    - Supports retaining specific sequences via --keep-sequence-id
    - Cluster sizes become sequence weights for downstream analysis
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("alignment", help="Input MSA (FASTA/Stockholm/Clustal; optionally .gz).")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (path without extension).")
    ap.add_argument("--min-seq-id", type=float, default=0.85, help="MMseqs2 identity threshold (default 0.85).")
    ap.add_argument("-c", "--coverage", type=float, default=0.8, help="MMseqs2 coverage threshold (default 0.8).")
    ap.add_argument("--cov-mode", type=int, default=0, choices=[0,1,2], help="MMseqs2 cov-mode (default 0, bidirectional).")
    ap.add_argument("--keep-tmp", action="store_true", help="Keep MMseqs2 temp directory (for debugging).")
    ap.add_argument("--keep-sequence-id", type=str, action="append", default=None,
                   help="Sequence ID to ensure is retained as representative (first part of header, before whitespace). "
                        "Can be specified multiple times to retain multiple sequences.")
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    tmpdir = out_prefix.parent / (out_prefix.name + "_mmseqs_tmp")
    if tmpdir.exists():
        shutil.rmtree(tmpdir)

    aln = read_alignment(args.alignment)

    ungapped_fa = tmpdir / "ungapped.fasta"
    tmpdir.mkdir(parents=True, exist_ok=True)
    write_ungapped_fasta(aln, ungapped_fa)

    tsv = run_mmseqs_cluster(ungapped_fa, tmpdir, args.min_seq_id, args.coverage, args.cov_mode)
    rep_sizes = compute_cluster_sizes(tsv)

    out_reps = Path(str(out_prefix) + "_reps.fasta")
    out_wts = Path(str(out_prefix) + "_weights.tsv")
    # args.keep_sequence_id is now a list (or None) due to action="append"
    keep_ids = args.keep_sequence_id if args.keep_sequence_id else None
    write_reps_alignment(aln, rep_sizes, out_reps, out_wts, keep_sequence_ids=keep_ids)

    print(f"Wrote: {out_reps}")
    print(f"Wrote: {out_wts}")
    print(f"Clusters: {len(rep_sizes)}; total members counted: {sum(rep_sizes.values())}")

    if not args.keep_tmp:
        shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == "__main__":
    main()
