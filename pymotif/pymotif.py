#!/usr/bin/env python

"""
Command-line script to perform motif finding

Similar to HOMER
"""

import argparse
from . import utils
from pymotif import __version__
import os
import sys

def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        prog="pymotif",
        description="Command-line script to perform motif finding ...",
    )

    # Required positional arguments
    parser.add_argument("peak_file", help="Peak or BED file", type=str)
    parser.add_argument("genome", help="Genome fasta file", type=str)
    parser.add_argument("output_directory", help="Directory to output results", type=str)

    # Optional arguments
    parser.add_argument("-size", help="Size of motif to search for", type=int, default=200)
    # 'store_true' makes this False by default but can be made True by including '-mask' in the command
    parser.add_argument("-mask", help="Use to mask the genome", action='store_true') 

    # Parse the arguments
    args = parser.parse_args()

    print('PyMotif!')

    # You can then use these arguments in your script like this:
    print(f"Peak file: {args.peak_file}")
    print(f"Genome: {args.genome}")
    print(f"Output directory: {args.output_directory}")
    print(f"Motif size: {args.size}")
    print(f"Mask genome: {args.mask}")
    
    
    sys.exit(0)

    # # Output
    # parser.add_argument("-o", "--out", help="Write output to file. " \
    # 	"Default: stdout", metavar="FILE", type=str, required=False)

    # # Other options
    # parser.add_argument("-f", "--fasta-ref", \
    # 	help="faidx indexed reference sequence file", \
    # 	type=str, metavar="FILE", required=False)
    # parser.add_argument("-r", "--region", help="region in which pileup is " \
    # 	"generated. Format chr:start-end", \
    # 	type=str, metavar="REG", required=False)
    # parser.add_argument("--version", help="Print the version and quit", \
    # 	action="version", version = '{version}'.format(version=__version__))

    # # Parse args
    # args = parser.parse_args()

    # # Set up output file
    # if args.out is None:
    # 	outf = sys.stdout
    # else: outf = open(args.out, "w")

    # # Load FASTA
    # if args.fasta_ref is not None:
    # 	if not os.path.exists(args.fasta_ref):
    # 		myutils.ERROR("{fasta} does not exist".format(fasta=args.fasta_ref))
    # 	reffasta = pyfaidx.Fasta(args.fasta_ref)
    # else:
    # 	reffasta = None

    # # Load BAM
    # bamfile = pysam.AlignmentFile(args.bam, "rb")

    # if args.region is not None:
    # 	region = args.region
    # else:
    # 	region = None

    # # Peform pileup
    # for pileupcolumn in bamfile.pileup(region=region):
    # 	chrom = pileupcolumn.reference_name
    # 	position = pileupcolumn.reference_pos
    # 	numreads = 0
    # 	if reffasta is not None:
    # 		refbase = str(reffasta[chrom][position])
    # 	else:
    # 		refbase = "N"
    # 	read_bases = []
    # 	read_quals = []
    # 	for pileupread in pileupcolumn.pileups:
    # 		if not pileupread.is_del and not pileupread.is_refskip:
    # 		# query position is None if is_del or is_refskip is set.
    # 			read_bases.append(myutils.GetReadBase(pileupread, refbase))
    # 			read_quals.append(myutils.GetReadQual(pileupread))
    # 			numreads += 1
    # 	outf.write("\t".join([chrom, str(position+1), refbase, str(numreads), \
    # 			"".join(read_bases), "".join(read_quals)])+"\n")
    # bamfile.close()
    # outf.close()
    

if __name__ == "__main__":
    main()