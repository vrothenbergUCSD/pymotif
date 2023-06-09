# PyMotif

PyMotif is a bioinformatics tool designed to identify enriched motifs from DNA sequences, similar to the HOMER tool. This tool is a Python script that implements a consensus string approach to motif finding.

# Install instructions

Installation requires the `BioPython`, `numpy`, `pandas`, `pyfaidx`, `pysam`, and `scipy` libraries to be installed. You can install these with `pip`:

```
pip install BioPython numpy pandas pyfaidx pysam scipy
```

Once required libraries are installed, you can install `pymotif` with the following command:

```
python setup.py install
```

Note: if you do not have root access, you can run the commands above with additional options to install locally:
```
pip install --user BioPython numpy pandas pyfaidx pysam scipy
python setup.py install --user
```
You may also need to run
```
PATH+=':<directory location>'
```

If the installation was successful, typing `pymotif --help` should show a useful message.

# Basic usage

The basic usage of `pymotif` is:

```
pymotif <peaks file> <genome> <output directory> -size # [options]
```

The following command assumes you are in the repo directory and have downloaded `GRCm38.fa` to the `example-files` folder. To download `GRCm38.fa`, you can visit https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.20/. If you are using a UCSD JupyterHub terminal for CSE 185 SP23, you can run
```
cp ~/public/genomes/GRCm38.fa ~/pymotif/example-files
```
You may also need to modify file paths to match your directories.

To run `pymotif` on a small test example (using files in this repo):
```
pymotif example-files/oct4_peaks.txt example-files/GRCm38.fa ~/testOutput/ -size 200 -mask
```

When the program begins, it should output statistics about usage, including version, number of sequences, and number of background sequences. As the program runs, it will output the motif it is evaluating along with its elapsed time. The final product will be in the `knownResults.txt` file which includes information about the motifs passing the p-value threshold, including the motif name, consensus sequence, p-value, and other statistics. Please note that specifying number of cores (see options below) is strongly recommended for faster results.

To compare to output of HOMER run:
```
findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
```

For example, to compare the small test example to HOMER run:
````
prefix=Oct4
findMotifsGenome.pl \
 ~/lab5/tagdirs/${prefix}/peaks.txt \
 ~/public/genomes/GRCm38.fa \
 ~/lab5/motifs/${prefix} \
 -mask -size 100
 ````
 
 To run `pymotif` on an external dataset, please see the `experimental-data` folder.

# PyMotif options

PyMotif takes the following command line arguments:

## Required arguments

- `peak_file`: The path to your peak or BED file.
- `genome`: The path to your genome fasta file. For example, you can use the mouse genome GCF_000001635.27_GRCm39_genomic.fna by downloading the assembly from NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27
- `output_directory`: The path to the directory where you want PyMotif to output its results.

Example usage with required arguments:

```
pymotif example-files/oct4_peaks.txt example-files/GRCm38.fa ~/motif_output/
```

## Optional arguments

- `-size`: This is the size of the motif to search for. It defaults to 100 if not specified.
- `-motifs`: The path to a known motifs file.
- `-cores`: Number of cores to use. It defaults to 1 if not specified, but it is strongly recommended that you increase the number of cores to maximize efficiency and runtime.
- `-enrichment_pval_thresh`: Enrichment threshold p-value for motif finding. It defaults to 0.00001 if not specified.
- `-num_peaks`: Number of peaks to use out of peaks file, default all.
- `--version`: Print the version and quit.

Example usage with optional arguments:

```
pymotif example-files/oct4_peaks.txt example-files/GRCm38.fa ~/motif_output/ -size 200 -cores 4
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
