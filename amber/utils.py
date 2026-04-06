"""
amber.utils
============
Shared helper functions used across the AMBER pipeline.

Includes DNA-sequence utilities, a terminal progress bar, and
functions for reading/writing AMBER's tab-delimited input/output files.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from os.path import exists
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pysam

# ---------------------------------------------------------------------------
# Pre-compiled regex for MD-tag tokenisation (compile once, reuse everywhere)
# ---------------------------------------------------------------------------
MD_TOKEN_RE = re.compile(r"(\d+|[A-Za-z]|\^[A-Za-z]+)")

# ---------------------------------------------------------------------------
# DNA complement table – used by reverse_complement()
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    """Return the reverse-complement of a DNA string.

    Uses ``str.translate`` + ``str.maketrans`` for speed — roughly 3–5×
    faster than the original dict-lookup approach on typical aDNA read
    lengths (30–150 bp).

    Parameters
    ----------
    seq : str
        Input DNA sequence (upper- or lower-case).

    Returns
    -------
    str
        Reverse-complemented sequence in the same case as *seq*.
    """
    return seq.translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# Terminal progress bar
# ---------------------------------------------------------------------------
def print_progress_bar(
    iteration: int,
    total: int,
    prefix: str = "",
    suffix: str = "",
    decimals: int = 1,
    length: int = 50,
    fill: str = "█",
) -> None:
    """Display or update a single-line progress bar on *stdout*.

    Parameters
    ----------
    iteration : int
        Current step (0-based or 1-based — just be consistent).
    total : int
        Total number of steps.  Must be > 0.
    prefix, suffix : str
        Optional text flanking the bar.
    decimals : int
        Number of decimal places in the percentage.
    length : int
        Character width of the bar itself.
    fill : str
        Character used for the filled portion of the bar.
    """
    if total <= 0:
        return
    pct = f"{100 * iteration / total:.{decimals}f}"
    filled = int(length * iteration // total)
    bar = fill * filled + "-" * (length - filled)
    sys.stdout.write(f"\r{prefix} |{bar}| {pct}% {suffix}")
    sys.stdout.flush()
    if iteration >= total:
        sys.stdout.write("\n")


# ---------------------------------------------------------------------------
# BAM-file helpers
# ---------------------------------------------------------------------------
def ensure_bam_index(filepath: str) -> None:
    """Create a ``.bai`` index for a BAM file if one does not exist.

    Parameters
    ----------
    filepath : str
        Path to a coordinate-sorted BAM file.
    """
    if not exists(filepath + ".bai"):
        print("  Creating BAM index …", end=" ", flush=True)
        pysam.index(filepath)
        print("done.")


def load_contigs_to_exclude(exclude_path: Optional[str]) -> Set[str]:
    """Read a newline-delimited list of contig names to skip.

    The unmapped-read sentinel ``*`` is always excluded.

    Parameters
    ----------
    exclude_path : str or None
        Path to the exclusion file.  ``None`` means "exclude nothing".

    Returns
    -------
    set[str]
        Set of contig names to ignore during analysis.
    """
    excluded: Set[str] = {"*"}
    if exclude_path is not None:
        with open(exclude_path) as fh:
            for line in fh:
                name = line.strip()
                if name:
                    excluded.add(name)
    return excluded


def parse_bamfile_list(path: str) -> List[Tuple[str, str]]:
    """Parse the user-supplied tab-separated BAM-file list.

    Each non-empty line must have at least two tab-separated columns:
    ``sample_name<TAB>path_to_bam``.

    Parameters
    ----------
    path : str
        Path to the BAM-file list.

    Returns
    -------
    list of (sample_name, bam_path) tuples

    Raises
    ------
    SystemExit
        If the file is empty or contains more than 6 entries (plotting
        limit imposed by the colour palette).
    """
    entries: List[Tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"  WARNING: skipping malformed line: {line!r}")
                continue
            entries.append((parts[0].strip(), parts[1].strip()))

    if not entries:
        sys.exit("ERROR: BAM-file list is empty.")
    if len(entries) > 6:
        sys.exit("ERROR: A maximum of 6 BAM files is supported per run.")
    return entries


# ---------------------------------------------------------------------------
# MD-tag parsing (optimised)
# ---------------------------------------------------------------------------
def reconstruct_reference(
    read_seq: str, md_tag: str
) -> Tuple[str, str]:
    """Rebuild the reference sequence from a read and its MD tag.

    Only handles reads whose CIGAR is a pure match (no indels, clips, or
    padding).  The caller must pre-filter on CIGAR.

    Parameters
    ----------
    read_seq : str
        The observed query sequence from the BAM record.
    md_tag : str
        The ``MD:Z:`` tag value.

    Returns
    -------
    (ref_seq, read_seq_aligned) : tuple[str, str]
        The reconstructed reference and the aligned portion of the read
        (identical to *read_seq* when there are no deletions).
    """
    tokens = MD_TOKEN_RE.findall(md_tag)
    ref_parts: list[str] = []
    read_parts: list[str] = []
    pos = 0  # current position in read_seq

    for tok in tokens:
        if tok[0] == "^":
            # Deletion from reference — skip (read has no bases here).
            continue
        if tok.isdigit():
            n = int(tok)
            ref_parts.append(read_seq[pos : pos + n])
            read_parts.append(read_seq[pos : pos + n])
            pos += n
        else:
            # Single-base mismatch(es) — each letter is one mismatch.
            for base in tok:
                ref_parts.append(base)
                read_parts.append(read_seq[pos])
                pos += 1

    return "".join(ref_parts), "".join(read_parts)


# ---------------------------------------------------------------------------
# Dataclass for per-sample results (lightweight container)
# ---------------------------------------------------------------------------
@dataclass
class SampleResult:
    """Container holding all analysis results for a single BAM file.

    Fields are populated by the coverage and damage analysis modules and
    consumed by the plotting and output-writing modules.
    """

    name: str

    # -- Coverage ----------------------------------------------------------
    window_depths: np.ndarray = field(default_factory=lambda: np.array([]))
    avg_depth: float = 0.0

    # -- Mismatch rates (keyed by read length) ----------------------------
    mismatch_lengths: np.ndarray = field(default_factory=lambda: np.array([]))
    mismatch_means: np.ndarray = field(default_factory=lambda: np.array([]))
    mismatch_ci95: np.ndarray = field(default_factory=lambda: np.array([]))

    # -- Read-length distribution -----------------------------------------
    readlen_lengths: np.ndarray = field(default_factory=lambda: np.array([]))
    readlen_counts: np.ndarray = field(default_factory=lambda: np.array([]))

    # -- DNA damage -------------------------------------------------------
    damage_ct: np.ndarray = field(default_factory=lambda: np.zeros(31))
    damage_cpg: np.ndarray = field(default_factory=lambda: np.zeros(31))
    damage_other: np.ndarray = field(default_factory=lambda: np.zeros(31))
