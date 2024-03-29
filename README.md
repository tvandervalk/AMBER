# AMBER
Asses Mapping Biases and Evaluating Read Reliability

## Description
AMBER computes sequence read mapping bias, ancient DNA damage patterns, fragment length distribution and genome breadth of coverage directly from a BAM-file. 
No external packages except for python3, numpy, matplotlib and pysam are required.

## Installation

```
git clone https://github.com/tvandervalk/AMBER.git
cd AMBER
chmod +x AMBER
```

If matplotlib,numpy and/or pysam are not yet installed:
```
python -m pip install -U pysam
python -m pip install -U matplotlib
python -m pip install -U numpy
```

## Input files

To assess your ancient samples, a tab seperated file with the bamfiles to be analyzed has to be supplied. The First column of the file should contain sample name, second column should contain path to bamfile. A maximum of 6 bamfiles can be analysed and plotted in the same run (however it is recommended to analyse just one bamfile at the time).

### Required input file
```
sampleA  /path/to/sampleA.bam
sampleB  /path/to/sampleB.bam
sampleC  /path/to/sampleC.bam
....
```

### Optional input file
A text file containing per line the name of a chromosomes/scaffolds or contigs that should be excluded from the analysis. For example:
```
contig001
contig002
chrX
chrY
```

## Run AMBER

```
./AMBER
--bamfiles         (default = empty)            Tab seperated file containing the sample names and path to bamfiles to be included in the analysis, required
--output           (default = amber)            Name of the output files generated by AMBER, optional
--exclude          (default = empty)            A text file containing per line the name of chromosomes/scaffolds or contigs that should be excluded from the analysis, optional
--errorbars        (default = not set)          plot the 95% confidence intervals in the mistmatch rate plot
--counts           (default = not set)          Plot the data as counts instead of fractions
```

## AMBER output example for a ~1 million year old Mammuthus genome (UDG treated)
<img width="1787" alt="Screenshot 2024-01-30 at 08 05 46" src="https://github.com/tvandervalk/AMBER/assets/64150803/f4d7ef24-0c0d-4a75-bf10-60a6d4431f2b">
