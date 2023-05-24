# PyMotif

(Work in progress!)

PyMotif is a bioinformatics tool designed to identify enriched motifs from DNA sequences, similar to the Homer tool. This tool is a Python script that implements a consensus string approach to motif finding.

# Install instructions

Todo

# Basic usage

The basic usage of `pymotif` is:

```
pymotif <peaks file> <genome> <output directory> -size # [options]
```

To run `pymotif` on a small test example (using files in this repo):
```
pymotif ERpeaks.txt hg18 ER_MotifOutput/ -size 200 -mask
```

This should produce the output below:
```
...
```

To compare to output of HOMER run:
```
findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
```

# PyMotif options

PyMotif takes the following command line arguments:

## Required arguments

- `peak_file`: This is the path to your peak or BED file.
- `genome`: This is the path to your genome fasta file.
- `output_directory`: This is the path to the directory where you want PyMotif to output its results.

Example usage with required arguments:

pymotif ERpeaks.txt hg18 ER_MotifOutput/

## Optional arguments

- `-size`: This is the size of the motif to search for. It defaults to 200 if not specified.
- `-mask`: If included in the command, PyMotif will mask the genome. 

Example usage with optional arguments:

pymotif ERpeaks.txt hg18 ER_MotifOutput/ -size 200 -mask

# File format

The output file format is the same as the HOMER knownResults.txt file.  See: http://homer.ucsd.edu/homer/ngs/peakMotifs.html 

# Contributors

This repository was created by Vince Rothenberg and Caitlyn Truong, as part of the final project for CSE 185 at UC San Diego.

Please submit a pull request with any corrections or suggestions.

# Testing

To run tests:
```
# Run command line tests
sh tests/cmdline_tests.sh

# Run unit tests
python -m pytest --cov=.
```



