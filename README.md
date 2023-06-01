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
pip install --user pandas pyfaidx pysam BioPython
python setup.py install --user
```

If the installation was successful, typing `pymotif --help` should show a useful message.

# Basic usage

The basic usage of `pymotif` is:

```
pymotif <peaks file> <genome> <output directory> -size # [options]
```

The following command assumes you are in the repo directory and have downloaded GRCm38.fa to the example-files folder.

To run `pymotif` on a small test example (using files in this repo):
```
cp ~/public/genomes/GRCm38.fa ~/pymotif/example-files
pymotif example-files/peaks.txt example-files/GRCm38.fa ~/testOutput/ -size 200 -mask
```

Note, the first command assumes you are using a UCSD JupyterHub terminal for CSE 185. Please visit https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/ to download `GRCm38.fa` otherwise. You may also need to modify file paths to match your directories.

When the program begins, it should output statistics about usage, including version, number of sequences, and number of background sequences. As the program runs, it will output the motif it is evaluating along with its elapsed time. The final product will be in the `knownResults.txt` file which includes information about the motifs passing the p-value threshold, including the motif name, consensus sequence, p-value, and other statistics. Please note that the program currently takes quite some time to complete (about an hour); we are working on fixing this in a future version of our tool.

To compare to output of HOMER run:
```
findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
```

# PyMotif options

PyMotif takes the following command line arguments:

## Required arguments

- `peak_file`: The path to your peak or BED file. This file may be generated using Homer's `findPeaks` function; from a fastqc file of reads, generate an aligned bam file and create tag directories through Homer (`makeTagDirectory`).
- `genome`: The path to your genome fasta file. For example, you can use the mouse genome GCF_000001635.27_GRCm39_genomic.fna by downloading the assembly from NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27
- `output_directory`: The path to the directory where you want PyMotif to output its results.

Example usage with required arguments:

```
pymotif ERpeaks.txt hg18 ER_MotifOutput/
```

## Optional arguments

- `-size`: This is the size of the motif to search for. It defaults to 200 if not specified.
- `-mask`: If included in the command, PyMotif will mask the genome. 
- `-motifs`: The path to a known motifs file.

Example usage with optional arguments:

```
pymotif ERpeaks.txt hg18 ER_MotifOutput/ -size 200 -mask
```

# File format

The output file format is the same as the HOMER knownResults.txt file.  It includes column-formatted information about the motif name, consensus sequence,	p-value,	log p-value,	number of target sequences with the motif, percentage of target sequences with the motif, number of background sequences with the motif, and percentage of background sequences with the motif.

See: http://homer.ucsd.edu/homer/ngs/peakMotifs.html

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



