"""
amber.indels
==============
Indel-rate analysis.

Iterates over every mapped read in the BAM file to compute:

1. **Per-read-length indel counts** — the mean number of indel events 
   (insertions + deletions) per read, grouped by read length, with 95 % 
   confidence intervals.
2. **Fragment-length distribution** — counts of reads at each length.

Performance notes
-----------------
* Uses pysam's numeric ``cigartuples`` attribute to quickly scan for 
  indel operations without string parsing.
* Accumulates per-read metrics in flat lists and uses vectorised NumPy 
  operations for the final aggregation step.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Set

import numpy as np
import pysam

from .utils import print_progress_bar

# ── CIGAR operation codes for indels ──────────────────────────────────────
# 1 = I (insertion), 2 = D (deletion)
_INDEL_CIGAR_OPS = frozenset({1, 2})

# Maximum read length tracked (longer reads are binned here).
_MAX_RL = 300


@dataclass
class IndelResult:
    """Stores the aggregated results of the indel analysis."""
    name: str
    readlen_lengths: np.ndarray
    readlen_counts: np.ndarray
    indel_lengths: np.ndarray  # <-- Add this line
    indel_means: np.ndarray
    indel_ci95: np.ndarray


def analyse_indels(
    filepath: str,
    excluded: Set[str],
    total_reads: int,
) -> IndelResult:
    """Run the indel-rate analysis on a single BAM.

    Parameters
    ----------
    filepath : str
        Path to a sorted, indexed BAM file.
    excluded : set[str]
        Contig names to skip.
    total_reads : int
        Total mapped reads (from ``analyse_coverage``), used only for the
        progress bar.

    Returns
    -------
    IndelResult
        Aggregated indel rates and read length distributions.
    """
    bam = pysam.AlignmentFile(filepath, "rb")

    # ── Pre-allocate data structures ──────────────────────────────────────
    # Indel count per read length: we collect lists of per-read values
    # and convert to mean/CI after the loop.
    indel_accum: List[List[int]] = [[] for _ in range(_MAX_RL + 1)]

    # Read-length histogram
    readlen_counts = np.zeros(_MAX_RL + 1, dtype=np.int64)

    # ── Main read loop ────────────────────────────────────────────────────
    print("  Calculating indel events per read …")
    read_counter = 0

    for read in bam.fetch():
        # -- Basic quality filters ----------------------------------------
        if read.is_unmapped or read.mapping_quality == 0:
            continue
        if read.reference_name in excluded:
            continue
            
        read_counter += 1
        if read_counter % 50_000 == 0:
            print_progress_bar(
                min(read_counter, total_reads),
                total_reads,
                prefix="  Progress:",
            )

        # Get read length (query_length is the sequence length)
        read_length = read.query_length
        if read_length == 0 or read.cigartuples is None:
            continue

        rl_bin = min(read_length, _MAX_RL)

        # -- Count indels from CIGAR tuples -------------------------------
        # cigartuples format: list of (operation, length)
        indel_count = sum(
            1 for op, _ in read.cigartuples if op in _INDEL_CIGAR_OPS
        )

        readlen_counts[rl_bin] += 1
        indel_accum[rl_bin].append(indel_count)

    print_progress_bar(total_reads, total_reads, prefix="  Progress:")
    bam.close()

    # ── Aggregate indel statistics ────────────────────────────────────────
    id_lengths: list[int] = []
    id_means: list[float] = []
    id_ci95: list[float] = []

    for rl in range(_MAX_RL + 1):
        vals = indel_accum[rl]
        if len(vals) > 10:  # Only calculate stats for lengths with >10 reads
            arr = np.array(vals, dtype=np.float64)
            id_lengths.append(rl)
            id_means.append(float(arr.mean()))
            
            # Calculate 95% confidence interval
            id_ci95.append(float(1.96 * arr.std() / np.sqrt(len(arr))))

    # ── Build read-length arrays (only lengths with reads) ───────────────
    rl_mask = readlen_counts > 0
    rl_lengths = np.where(rl_mask)[0]
    rl_vals = readlen_counts[rl_mask]

# ── Package results ──────────────────────────────────────────────────
    result = IndelResult(
        name="",
        readlen_lengths=rl_lengths,
        readlen_counts=rl_vals,
        indel_lengths=np.array(id_lengths, dtype=np.int64),  
        indel_means=np.array(id_means, dtype=np.float64),
        indel_ci95=np.array(id_ci95, dtype=np.float64),
    )

    return result