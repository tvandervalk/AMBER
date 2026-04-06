"""
AMBER — Assess Mapping Biases and Evaluate Read Reliability.

A modular toolkit for analysing BAM files from ancient-DNA sequencing
experiments. AMBER computes per-read mismatch rates, cytosine-deamination
damage profiles, fragment-length distributions, and genome-wide coverage
histograms.

Modules
-------
coverage   : Depth-of-coverage analysis in 1 kb windows.
damage     : Mismatch-rate and DNA-damage pattern analysis.
plotting   : Matplotlib-based visualisation of all statistics.
utils      : Shared helper functions (reverse complement, progress bar, I/O).
"""

__version__ = "2.0.0"
