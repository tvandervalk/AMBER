# AMBER v2.0

**A**ssess **M**apping **B**iases and **E**valuate **R**ead **R**eliability

## Description

AMBER computes sequence-read mapping bias, ancient DNA damage patterns,
fragment-length distributions, and genome breadth of coverage directly from
BAM files. No external packages beyond Python 3, NumPy, Matplotlib, and
pysam are required.


## Installation

```bash
git clone https://github.com/tvandervalk/AMBER.git
cd AMBER
```

Install dependencies (if needed):

```bash
python -m pip install -U pysam matplotlib numpy
```

## Usage

```bash
python -m amber --bamfiles bamfile_list.txt --output my_run
```

### Options

| Flag | Default | Description |
|---|---|---|
| `--bamfiles` | *(required)* | Tab-separated file: column 1 = sample name, column 2 = BAM path |
| `--output` | `amber` | Output file-name prefix |
| `--exclude` | none | File listing contig names to exclude (one per line) |
| `--errorbars` | off | Plot 95% confidence intervals on the mismatch-rate panel |
| `--counts` | off | Plot raw counts instead of fractions |

### Input files

**Required:** a tab-separated BAM-file list (max 6 samples):

```
sampleA	/path/to/sampleA.bam
sampleB	/path/to/sampleB.bam
```

**Optional:** a contig-exclusion list (one name per line):

```
contig001
chrX
chrY
```

## Package structure

```
AMBER/
├── amber/
│   ├── __init__.py      # Package metadata and version
│   ├── __main__.py      # CLI entry point and pipeline orchestration
│   ├── coverage.py      # Depth-of-coverage analysis (1 kb windows)
│   ├── damage.py        # Mismatch rate and DNA-damage profiling
│   ├── output.py        # Tab-delimited text report writer
│   ├── plotting.py      # Four-panel PDF figure generation
│   └── utils.py         # Shared helpers (reverse complement, I/O, etc.)
├── AMBER_original       # Original v1 script (preserved for reference)
├── LICENSE.md
└── README.md
```

## Output

AMBER produces two files:

- **`{output}.pdf`** — Four-panel figure: mismatch rate vs read length,
  positional DNA-damage profile, fragment-length distribution, and
  coverage-depth histogram.
- **`{output}.txt`** — Tab-delimited statistics for each sample.

## Citation

TBA
