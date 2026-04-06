"""
amber.damage
==============
Mismatch-rate and DNA-damage pattern analysis.

Iterates over every mapped read in the BAM file to compute:

1. **Per-read-length mismatch rate** — the fraction of mismatching bases
   per read, grouped by read length, with 95 % confidence intervals.
2. **Positional DNA-damage profile** — the frequency of C→T, CpG→TpG,
   and "other" mismatches at each of the first 31 positions from the
   read end (after orienting to the original 5′→3′ direction).
3. **Fragment-length distribution** — counts of reads at each length.

Performance notes
-----------------
* Memory tracking for mismatches uses Welford-like running sums of squares
  via NumPy arrays, completely eliminating list-append overhead.
* Sequence reconstruction is delegated to pysam's C-compiled 
  `get_reference_sequence()`.
* String slicing is applied *before* reverse-complementing to avoid 
  processing full 150+ bp read sequences.
* String allocations in the inner loop have been removed.
"""

from __future__ import annotations

from typing import Set

import numpy as np
import pysam

from .utils import SampleResult, print_progress_bar

# ── CIGAR operation codes that disqualify a read ──────────────────────────
# 1=I (ins), 2=D (del), 3=N (ref-skip), 4=S (soft-clip),
# 5=H (hard-clip), 6=P (pad)
_BAD_CIGAR_OPS = frozenset({1, 2, 3, 4, 5, 6})

# Maximum read length tracked (longer reads are binned here).
_MAX_RL = 300

# Number of positions from the read end analysed for damage.
_DAMAGE_POSITIONS = 31

# Fast translation table for reverse complement
_RC_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _has_bad_cigar(read: pysam.AlignedSegment) -> bool:
    """Return True if the read's CIGAR contains indels, clips, or padding."""
    if read.cigartuples is None:
        return True
    # Fast path: if the read is a single continuous match (op 0)
    if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0:
        return False
    return any(op in _BAD_CIGAR_OPS for op, _ in read.cigartuples)


def analyse_damage(
    filepath: str,
    excluded: Set[str],
    total_reads: int,
) -> SampleResult:
    bam = pysam.AlignmentFile(filepath, "rb")

    # ── Pre-allocate Memory-Efficient Arrays ──────────────────────────────
    readlen_counts = np.zeros(_MAX_RL + 1, dtype=np.int64)
    
    # Running sums for Mismatch Rate (avoids storing millions of floats)
    mm_sum = np.zeros(_MAX_RL + 1, dtype=np.float64)
    mm_sum_sq = np.zeros(_MAX_RL + 1, dtype=np.float64)

    # DNA-damage accumulators (numerators)
    dmg_ct = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)
    dmg_cpg = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)
    dmg_other = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)

    # Denominators (site counts)
    sites_c = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)
    sites_cpg = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)
    sites_other = np.zeros(_DAMAGE_POSITIONS, dtype=np.int64)

    # ── Main read loop ────────────────────────────────────────────────────
    print("  Calculating mismatch rate and DNA damage …")
    read_counter = 0

    for read in bam.fetch():
        if read.mapping_quality == 0 or read.reference_name in excluded:
            continue

        read_counter += 1
        if read_counter % 50_000 == 0:
            print_progress_bar(
                min(read_counter, total_reads), total_reads, prefix="  Progress:"
            )

        read_seq = read.query_sequence
        if read_seq is None or "N" in read_seq:
            continue

        if _has_bad_cigar(read):
            continue

        try:
            # Reconstruct reference using pysam's optimized C implementation
            ref_seq = read.get_reference_sequence()
        except ValueError:
            # Raised if MD tag is missing
            continue

        read_length = len(read_seq)
        rl_bin = min(read_length, _MAX_RL)

        # -- Mismatch rate ------------------------------------------------
        try:
            edit_dist = read.get_tag("NM")
        except KeyError:
            continue
            
        mm_val = edit_dist / read_length
        readlen_counts[rl_bin] += 1
        mm_sum[rl_bin] += mm_val
        mm_sum_sq[rl_bin] += mm_val * mm_val

        # -- Positional damage analysis (first 31 bp) --------------------
        seq_len = len(ref_seq)
        limit = min(_DAMAGE_POSITIONS, seq_len - 1)
        if limit <= 0:
            continue

        # Slice string BEFORE reverse-complementing to save CPU cycles
        if read.is_reverse:
            r_sub = ref_seq[-(limit + 1):].upper().translate(_RC_TRANS)[::-1]
            q_sub = read_seq[-(limit + 1):].upper().translate(_RC_TRANS)[::-1]
        else:
            r_sub = ref_seq[:limit + 1].upper()
            q_sub = read_seq[:limit + 1].upper()

        for pos in range(limit):
            r_base = r_sub[pos]
            q_base = q_sub[pos]

            if r_base == "N" or q_base == "N":
                continue

            if r_base == "C":
                if r_sub[pos + 1] == "G":
                    sites_cpg[pos] += 1
                    if q_base == "T":
                        dmg_cpg[pos] += 1
                else:
                    sites_c[pos] += 1
                    if q_base == "T":
                        dmg_ct[pos] += 1
            elif r_base not in "CT" and q_base not in "CT":
                # Evaluates logical 'Other' conditions without string concatenations
                if not (r_base == "G" and q_base == "A"):
                    sites_other[pos] += 1
                    if r_base != q_base:
                        dmg_other[pos] += 1

    print_progress_bar(total_reads, total_reads, prefix="  Progress:")
    bam.close()

    # ── Aggregate mismatch-rate statistics ────────────────────────────────
    mm_lengths = []
    mm_means = []
    mm_ci95 = []

    for rl in range(_MAX_RL + 1):
        n = readlen_counts[rl]
        if n > 10:
            mean = mm_sum[rl] / n
            # Variance derived from sum of squares: (Sum(x^2) - n * mean^2) / (n - 1)
            var = (mm_sum_sq[rl] - n * (mean ** 2)) / (n - 1)
            std = np.sqrt(max(0.0, var))  # max() guards against slight float inaccuracies
            
            mm_lengths.append(rl)
            mm_means.append(mean)
            mm_ci95.append(1.96 * std / np.sqrt(n))

    # ── Compute damage fractions (avoid division by zero) ────────────────
    with np.errstate(divide="ignore", invalid="ignore"):
        frac_ct = np.where(sites_c > 0, dmg_ct / sites_c, 0.0)
        frac_cpg = np.where(sites_cpg > 0, dmg_cpg / sites_cpg, 0.0)
        frac_other = np.where(sites_other > 0, dmg_other / sites_other, 0.0)

    # ── Build read-length arrays (only lengths with reads) ───────────────
    rl_mask = readlen_counts > 0
    rl_lengths = np.where(rl_mask)[0]
    rl_vals = readlen_counts[rl_mask]

    # ── Package results ──────────────────────────────────────────────────
    result = SampleResult(name="")
    result.mismatch_lengths = np.array(mm_lengths, dtype=np.int64)
    result.mismatch_means = np.array(mm_means, dtype=np.float64)
    result.mismatch_ci95 = np.array(mm_ci95, dtype=np.float64)
    result.readlen_lengths = rl_lengths
    result.readlen_counts = rl_vals
    result.damage_ct = frac_ct
    result.damage_cpg = frac_cpg
    result.damage_other = frac_other

    return result