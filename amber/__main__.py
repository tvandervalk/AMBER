"""
amber.__main__
===============
Command-line entry point for the AMBER pipeline.

Run with::

    python -m amber --bamfiles bamfile_list.txt --output my_run

This module wires together the analysis stages:

1. Parse CLI arguments and load the BAM-file list.
2. For each BAM file, run coverage analysis (``amber.coverage``).
3. For each BAM file, run damage/mismatch analysis (``amber.damage``).
4. For each BAM file, run indel analysis (``amber.indels``).
5. Generate the multi-panel PDF plot (``amber.plotting``) and the
   tab-delimited text report (``amber.output``).
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import List

from .utils import (
    SampleResult,
    ensure_bam_index,
    load_contigs_to_exclude,
    parse_bamfile_list,
)
from .coverage import analyse_coverage
from .damage import analyse_damage
from .indels import analyse_indels  # <-- Added import for indels
from .plotting import generate_plots
from .output import write_report

import pysam


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line argument parser.

    Replaces the deprecated ``optparse`` module used in AMBER v1.
    """
    parser = argparse.ArgumentParser(
        prog="amber",
        description=(
            "AMBER — Assess Mapping Biases and Evaluate Read Reliability.\n"
            "Computes mismatch rates, DNA-damage patterns, fragment-length "
            "distributions, coverage histograms, and indel rates from BAM files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bamfiles",
        required=True,
        help=(
            "Tab-separated file listing BAM files to analyse. "
            "Column 1 = sample name, column 2 = path to BAM."
        ),
    )
    parser.add_argument(
        "--output",
        default="amber",
        help="Output file-name prefix (default: 'amber').",
    )
    parser.add_argument(
        "--exclude",
        default=None,
        help=(
            "File listing contig names to exclude (one per line). "
            "Default: none."
        ),
    )
    parser.add_argument(
        "--errorbars",
        action="store_true",
        default=False,
        help="Show 95%% confidence intervals on the mismatch and indel plots.",
    )
    parser.add_argument(
        "--counts",
        action="store_true",
        default=False,
        help="Plot raw counts instead of fractions.",
    )
    return parser


def process_sample(
    sample_name: str,
    bam_path: str,
    excluded: set,
) -> SampleResult:
    """Run the full analysis pipeline for a single BAM file.

    Parameters
    ----------
    sample_name : str
        Human-readable label for this sample.
    bam_path : str
        Path to a coordinate-sorted BAM file.
    excluded : set[str]
        Contig names to skip.

    Returns
    -------
    SampleResult
        Fully populated result object.
    """
    print(f"\n{'─' * 70}")
    print(f"Processing sample: {sample_name}")
    t0 = time.time()

    # Validate BAM is sorted
    bam = pysam.AlignmentFile(bam_path, "rb")
    hd = bam.header.get("HD", {})
    if hd.get("SO") == "unsorted":
        sys.exit("ERROR: BAM file is unsorted. Only sorted BAMs are supported.")
    bam.close()

    # Ensure .bai index exists
    ensure_bam_index(bam_path)

    # Stage 1 — Coverage
    window_depths, avg_depth, total_reads = analyse_coverage(bam_path, excluded)

    # Stage 2 — Damage / mismatch
    result = analyse_damage(bam_path, excluded, total_reads)

    # Stage 3 — Indels (NEW)
    indel_result = analyse_indels(bam_path, excluded, total_reads)

    # Merge data into a single result object
    result.name = sample_name
    result.window_depths = window_depths
    result.avg_depth = avg_depth
    
    # Merge the new indel properties
    result.indel_lengths = indel_result.indel_lengths  # <-- CHANGE THIS LINE
    result.indel_means = indel_result.indel_means
    result.indel_ci95 = indel_result.indel_ci95

    elapsed = time.time() - t0
    print(f"  Sample completed in {elapsed:.1f} s")
    return result


def main(argv: list[str] | None = None) -> None:
    """Main entry point — parse arguments and run the pipeline."""
    args = build_parser().parse_args(argv)

    # ── Load inputs ──────────────────────────────────────────────────────
    bam_entries = parse_bamfile_list(args.bamfiles)
    excluded = load_contigs_to_exclude(args.exclude)

    # ── Process each BAM file ────────────────────────────────────────────
    results: List[SampleResult] = []
    for sample_name, bam_path in bam_entries:
        result = process_sample(sample_name, bam_path, excluded)
        results.append(result)

    # ── Generate outputs ─────────────────────────────────────────────────
    print(f"\n{'─' * 70}")
    print("Generating outputs …")
    generate_plots(
        results,
        args.output,
        show_errorbars=args.errorbars,
        show_counts=args.counts,
    )
    write_report(results, args.output, show_counts=args.counts)
    print("\nFinished.")


if __name__ == "__main__":
    main()