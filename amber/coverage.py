"""
amber.coverage
===============
Depth-of-coverage analysis for BAM files.

Computes per-window coverage across all included contigs using
``pysam.idxstats`` for contig metadata and ``pysam.fetch`` for
high-speed read aggregation. The genome is divided into exactly 10,000 windows.

Performance notes
-----------------
* Uses `bam.fetch()` instead of `bam.pileup()`. Iterating over reads rather 
  than individual covered bases yields a ~10x to 100x speedup.
* Replaced inner-loop division with pre-calculated multiplication.
* Pulled contig dictionary lookups out of the read-iteration loop.
"""

from __future__ import annotations

from typing import Dict, Set, Tuple

import numpy as np
import pysam

from .utils import print_progress_bar


def _parse_idxstats(
    filepath: str, excluded: Set[str]
) -> Tuple[Dict[str, int], Dict[str, float], int, int]:
    # (This function remains exactly the same as your updated version)
    contig_base_offsets: Dict[str, int] = {}
    contig_lengths: Dict[str, float] = {}
    total_reads = 0
    total_sites = 0

    raw = pysam.idxstats(filepath)
    for line in raw.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        
        contig_name = parts[0]
        if contig_name in excluded:
            continue
            
        length = int(parts[1])
        mapped = int(parts[2])

        if length == 0:
            continue

        contig_base_offsets[contig_name] = total_sites
        contig_lengths[contig_name] = float(length)

        total_reads += mapped
        total_sites += length

    return contig_base_offsets, contig_lengths, total_reads, total_sites


def analyse_coverage(
    filepath: str,
    excluded: Set[str],
) -> Tuple[np.ndarray, float, int]:
    
    bam = pysam.AlignmentFile(filepath, "rb")

    # -- Gather contig metadata from the index ----------------------------
    contig_base_offsets, contig_lengths, total_reads, total_sites = _parse_idxstats(
        filepath, excluded
    )
    if total_reads == 0 or total_sites == 0:
        raise RuntimeError("No reads or valid sites found on the analysed contigs.")

    TOTAL_WINDOWS = 10000
    window_size = total_sites / float(TOTAL_WINDOWS)
    
    # Pre-calculate a multiplier (multiplication is dramatically faster than division in Python)
    win_multiplier = TOTAL_WINDOWS / total_sites 
    
    depth_arr = np.zeros(TOTAL_WINDOWS, dtype=np.float64)

    # -- Fetch pass (Iterating reads instead of bases) --------------------
    print(f"  Calculating depth across 10,000 windows (approx {window_size:.0f} bp per window) …")
    
    reads_processed = 0

    # Iterate contigs in the outer loop so we don't do dict lookups per read
    for contig, offset in contig_base_offsets.items():
        
        # Iterate over READS, not bases. 
        for read in bam.fetch(contig=contig):
            if read.is_unmapped or read.mapping_quality < 1:
                continue
            
            # The length of the alignment on the reference genome
            ref_len = read.reference_length
            if ref_len is None:
                continue

            global_start = offset + read.reference_start
            global_end = offset + read.reference_end - 1
            
            # Map start and end to windows
            win_start_idx = int(global_start * win_multiplier)
            win_end_idx = int(global_end * win_multiplier)
            
            # Guard against bounds
            if win_start_idx >= TOTAL_WINDOWS: win_start_idx = TOTAL_WINDOWS - 1
            if win_end_idx >= TOTAL_WINDOWS: win_end_idx = TOTAL_WINDOWS - 1

            # 99.9% of the time, a read falls entirely within one large window
            if win_start_idx == win_end_idx:
                depth_arr[win_start_idx] += ref_len
            else:
                # If the read crosses a window boundary, split its bases proportionally
                window_boundary = (win_start_idx + 1) * window_size
                bases_in_first = int(window_boundary - global_start)
                bases_in_second = ref_len - bases_in_first
                
                depth_arr[win_start_idx] += bases_in_first
                depth_arr[win_end_idx] += max(0, bases_in_second) # max(0) safeguards against edge case float math

            # Progress reporting (every 1M reads instead of 1M bases)
            reads_processed += 1
            if reads_processed % 1_000_000 == 0:
                # Approximate progress based on reads to save computation
                print_progress_bar(reads_processed, total_reads, prefix="  Progress:")

    print_progress_bar(total_reads, total_reads, prefix="  Progress:")
    bam.close()

    # -- Convert summed base counts → average depth per window ------------
    depth_arr /= window_size

    # -- Compute summary statistics on covered windows --------------------
    covered_mask = depth_arr > 0
    if not np.any(covered_mask):
        raise RuntimeError("No reads found on the analysed contigs.")
    avg_depth = float(np.mean(depth_arr[covered_mask]))

    print(f"\n  Average depth (covered windows): {avg_depth:.4f}×")
    return depth_arr, avg_depth, total_reads