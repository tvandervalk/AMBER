"""
amber.output
==============
Write per-sample statistics to the AMBER text report.

Generates the same tab-delimited output format as the original AMBER
script so that downstream parsers continue to work unchanged.
"""

from __future__ import annotations

from typing import IO, List

import numpy as np

from .utils import SampleResult


def write_report(
    results: List[SampleResult],
    output_prefix: str,
    show_counts: bool = False,
) -> None:
    """Write the full statistics report for all samples.

    Parameters
    ----------
    results : list[SampleResult]
        One entry per BAM file.
    output_prefix : str
        The text file is saved as ``{output_prefix}.txt``.
    show_counts : bool
        If True, emit raw read counts in the read-length table; if
        False, emit percentages.
    """
    txt_path = f"{output_prefix}.txt"
    with open(txt_path, "w") as fh:
        for res in results:
            _write_sample(fh, res, show_counts)
    print(f"  Report saved → {txt_path}")


def _write_sample(
    fh: IO[str],
    res: SampleResult,
    show_counts: bool,
) -> None:
    """Emit one sample's statistics block."""
    sep = "-" * 100

    fh.write(f"{sep}\n")
    fh.write(f"sample: {res.name}\n")
    fh.write(f"average depth per window: {res.avg_depth:.4f}\n\n")

    # ── Coverage histogram ────────────────────────────────────────────────
    capped = np.clip(res.window_depths, 0, res.avg_depth * 5)
    hist, bin_edges = np.histogram(capped, bins=100)

    fh.write("WINDOW_DEPTH\tNUMBER_OF_WINDOWS\n")
    for count, edge in zip(hist, bin_edges[1:-1]):
        fh.write(f"{edge:.4f}\t{count}\n")
    fh.write(f">{bin_edges[-1]:.4f}\t{hist[-1]}\n")
    fh.write(f"{sep}\n")

    # ── Mismatch rate per read length ─────────────────────────────────────
    fh.write("READ_LENGTH\tMISMATCH_RATE\n")
    for rl, mm in zip(res.mismatch_lengths, res.mismatch_means):
        fh.write(f"{rl}\t{mm:.5f}\n")
    fh.write(f"{sep}\n")

    # ── Read-length distribution ──────────────────────────────────────────
    fh.write("READ_LENGTH\tREAD_COUNTS\n")
    total = res.readlen_counts.sum()
    for rl, cnt in zip(res.readlen_lengths, res.readlen_counts):
        if show_counts:
            fh.write(f"{rl}\t{cnt}\n")
        else:
            pct = (cnt / total) * 100.0 if total > 0 else 0.0
            fh.write(f"{rl}\t{pct:.5f}\n")
    fh.write(f"{sep}\n")

    # ── DNA damage ────────────────────────────────────────────────────────
    fh.write("DISTANCE_FROM_READ_END\tC-to-T\tCpG-to-TG\tother\n")
    for i in range(31):
        ct = res.damage_ct[i] if i < len(res.damage_ct) else 0.0
        cpg = res.damage_cpg[i] if i < len(res.damage_cpg) else 0.0
        oth = res.damage_other[i] if i < len(res.damage_other) else 0.0
        fh.write(f"{i}\t{ct:.4f}\t{cpg:.4f}\t{oth:.4f}\n")
    fh.write(f"\n{sep}\n")
