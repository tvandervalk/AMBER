"""
amber.plotting
===============
Matplotlib-based visualisation of AMBER analysis results.

Produces a single PDF with five subplots in a 3x2 grid:

1. **Mismatch rate vs read length** — scatter or error-bar plot.
2. **Positional DNA-damage frequency** — C→T, CpG→TpG, and "other".
3. **Fragment-length distribution** — filled area chart.
4. **Coverage histogram** — depth distribution of 1 kb windows.
5. **Indel rate vs read length** — mean indels per read.
"""

from __future__ import annotations

from typing import List

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from .utils import SampleResult

# Colour palette (max 6 samples)
COLOURS = ["dimgray", "darkred", "lightskyblue", "green", "gold", "darkorchid"]


def _configure_matplotlib() -> None:
    """Set global Matplotlib RC parameters for publication-quality plots."""
    plt.rcParams.update(
        {
            "figure.figsize": (20, 15),  # Increased height for 3 rows
            "figure.autolayout": True,
            "axes.linewidth": 1.25,
            "font.size": 14,
        }
    )


def generate_plots(
    results: List[SampleResult],
    output_prefix: str,
    show_errorbars: bool = False,
    show_counts: bool = False,
) -> None:
    """Render the five-panel AMBER figure and save it as a PDF."""
    _configure_matplotlib()
    
    # Create a 3x2 grid
    fig, axes = plt.subplots(3, 2)
    ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = axes
    
    # We only have 5 plots, so we hide the 6th (bottom right)
    ax6.set_visible(False)

    # ── Tracking variables for dynamic axis scaling ──────────────────────
    min_read_length = float("inf")
    max_read_length = 0.0
    max_damage = 0.0
    max_coverage_x = 0.0

    for idx, res in enumerate(results):
        colour = COLOURS[idx % len(COLOURS)]

        # ── 1. Mismatch rate vs read length ──────────────────────────────
        if hasattr(res, 'mismatch_lengths') and len(res.mismatch_lengths) > 0:
            min_read_length = min(min_read_length, res.mismatch_lengths.min())
            max_read_length = max(max_read_length, res.mismatch_lengths.max())

            if show_errorbars:
                ax1.errorbar(
                    res.mismatch_lengths, res.mismatch_means, yerr=res.mismatch_ci95,
                    fmt="o", markerfacecolor=colour, markersize=8,
                    markeredgecolor="black", label=res.name, color=colour,
                )
            else:
                ax1.scatter(
                    res.mismatch_lengths, res.mismatch_means,
                    s=40, alpha=1, edgecolors="k", label=res.name, color=colour,
                )

        # ── 2. DNA damage by position ────────────────────────────────────
        positions = np.arange(31)
        sample_max_damage = max(
            np.max(res.damage_other), np.max(res.damage_ct), np.max(res.damage_cpg)
        )
        max_damage = max(max_damage, sample_max_damage)

        ax2.plot(positions, res.damage_other, linewidth=2, linestyle="dotted", color=colour, label=f"{res.name}: other")
        ax2.plot(positions, res.damage_ct, linewidth=2, color=colour, label=f"{res.name}: C to T")
        ax2.plot(positions, res.damage_cpg, linewidth=2, linestyle="dashed", color=colour, label=f"{res.name}: CpG to TpG")

        # ── 3. Read-length distribution ──────────────────────────────────
        if hasattr(res, 'readlen_lengths') and len(res.readlen_lengths) > 0:
            min_read_length = min(min_read_length, res.readlen_lengths.min())
            max_read_length = max(max_read_length, res.readlen_lengths.max())

            if show_counts:
                y_vals = res.readlen_counts.astype(float)
            else:
                total = res.readlen_counts.sum()
                y_vals = (res.readlen_counts / total) * 100.0 if total > 0 else res.readlen_counts.astype(float)

            ax3.fill_between(res.readlen_lengths, 0, y_vals, alpha=0.25, label=res.name, color=colour)
            ax3.plot(res.readlen_lengths, y_vals, color=colour, linewidth=2, linestyle="-", marker="o", markersize=3)

        # ── 4. Coverage histogram ────────────────────────────────────────
        capped = np.clip(res.window_depths, 0, res.avg_depth * 5)
        hist, bin_edges = np.histogram(capped, bins=100)
        max_coverage_x = max(max_coverage_x, bin_edges[-1])

        if show_counts:
            y_hist = hist.astype(float)
        else:
            y_hist = hist / hist.sum() if hist.sum() > 0 else hist.astype(float)

        ax4.bar(bin_edges[1:], y_hist, width=bin_edges[1] - bin_edges[0], color=colour, alpha=0.5, label=res.name)
        ax4.axvline(x=res.avg_depth, linestyle="dashed", color=colour, zorder=0)

        # ── 5. Indel rate vs read length ─────────────────────────────────
        if hasattr(res, 'indel_lengths') and len(res.indel_lengths) > 0:
            # We already track min/max read lengths, but ensure it captures indel boundaries too
            min_read_length = min(min_read_length, res.indel_lengths.min())
            max_read_length = max(max_read_length, res.indel_lengths.max())

            if show_errorbars:
                ax5.errorbar(
                    res.indel_lengths, res.indel_means, yerr=res.indel_ci95,
                    fmt="o", markerfacecolor=colour, markersize=8,
                    markeredgecolor="black", label=res.name, color=colour,
                )
            else:
                ax5.scatter(
                    res.indel_lengths, res.indel_means,
                    s=40, alpha=1, edgecolors="k", label=res.name, color=colour,
                )

    # ── Axis formatting ──────────────────────────────────────────────────
    
    # Fallback if no valid read lengths were found
    if min_read_length == float("inf"):
        min_read_length, max_read_length = 0, 100

    # Read length margin
    rl_margin = max(5, (max_read_length - min_read_length) * 0.05)
    rl_xlim = (max(0, min_read_length - rl_margin), max_read_length + rl_margin)
    rl_ticker = ticker.MaxNLocator(nbins=10, integer=True)

    # Ax1 formatting (Mismatches)
    ax1.set_xlabel("read length (bp)", weight="bold", fontsize=14)
    ax1.set_ylabel("mismatches per basepair", weight="bold", fontsize=14)
    ax1.set_xlim(rl_xlim)
    ax1.xaxis.set_major_locator(rl_ticker)
    ax1.set_ylim(bottom=0)
    ax1.legend(loc="upper right", edgecolor="black")

    # Ax2 formatting (Damage)
    ax2.set_xlabel("distance from read end (bp)", weight="bold", fontsize=14)
    ax2.set_ylabel("mismatch frequency", weight="bold", fontsize=14)
    ax2.set_xlim(-0.5, 30.5)
    ax2.set_xticks(list(range(0, 31, 2)))
    ax2.set_ylim(0, max_damage * 1.1 if max_damage > 0 else 0.1)
    ax2.legend(loc="upper right", edgecolor="black")

    # Ax3 formatting (Fragment lengths)
    ax3.set_xlabel("read length (bp)", weight="bold", fontsize=14)
    ax3.set_ylabel("read counts" if show_counts else "% of reads", weight="bold", fontsize=14)
    ax3.set_xlim(rl_xlim)
    ax3.xaxis.set_major_locator(rl_ticker)
    ax3.set_ylim(bottom=0)
    ax3.legend(loc="upper right", edgecolor="black")

    # Ax4 formatting (Coverage)
    ax4.set_xlabel("depth", weight="bold", fontsize=14)
    ax4.set_ylabel("window counts" if show_counts else "fraction of 1 kb windows", weight="bold", fontsize=14)
    ax4.set_xlim(0, max_coverage_x * 1.05 if max_coverage_x > 0 else 10)
    ax4.xaxis.set_major_locator(ticker.MaxNLocator(nbins=10))
    ax4.set_ylim(bottom=0)
    ax4.legend(loc="upper right", edgecolor="black")

    # Ax5 formatting (Indels)
    ax5.set_xlabel("read length (bp)", weight="bold", fontsize=14)
    ax5.set_ylabel("mean indels per read", weight="bold", fontsize=14)
    ax5.set_xlim(rl_xlim)
    ax5.xaxis.set_major_locator(rl_ticker)
    ax5.set_ylim(bottom=0)
    ax5.legend(loc="upper right", edgecolor="black")

    # ── Save ─────────────────────────────────────────────────────────────
    pdf_path = f"{output_prefix}.pdf"
    plt.savefig(pdf_path)
    plt.close()
    print(f"\n  Plot saved → {pdf_path}")