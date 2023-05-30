# PyMotif

(Work in progress!)

PyMotif is a bioinformatics tool designed to identify enriched motifs from DNA sequences, similar to the Homer tool. This tool is a Python script that implements a consensus string approach to motif finding.

# Install instructions

Installation requires the `pandas`, `pyfaidx`, `pysam`, and `BioPython` libraries to be installed. You can install these with `pip`:

```
pip install pandas pyfaidx pysam BioPython
```

Once required libraries are installed, you can install `pymotif` with the following command:

```
python setup.py install
```

Note: if you do not have root access, you can run the commands above with additional options to install locally:
```
pip install --user pandas BioPython
python setup.py install --user
```

If the install was successful, typing `pymotif --help` should show a useful message.

# Basic usage

The basic usage of `pymotif` is:

```
pymotif <peaks file> <genome> <output directory> -size # [options]
```

The following command assumes you are in the repo directory and have downloaded GRCm38.fa to the example-files folder.

To run `pymotif` on a small test example (using files in this repo):
```
pymotif example-files/peaks.txt example-files/GRCm38.fa ~/testOutput/ -size 200 -mask
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

- `peak_file`: The path to your peak or BED file.
- `genome`: The path to your genome fasta file. For example, you can use the mouse genome GCF_000001635.27_GRCm39_genomic.fna by downloading the assembly from NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27
- `output_directory`: The path to the directory where you want PyMotif to output its results.

Example usage with required arguments:

pymotif ERpeaks.txt hg18 ER_MotifOutput/

## Optional arguments

- `-size`: This is the size of the motif to search for. It defaults to 200 if not specified.
- `-mask`: If included in the command, PyMotif will mask the genome. 
- `-motifs`: The path to a known motifs file.  Default is 

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



