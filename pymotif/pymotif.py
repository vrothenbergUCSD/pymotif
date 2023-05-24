#!/usr/bin/env python

"""
Command-line script to perform motif finding

Similar to HOMER
"""

import os
import argparse
from pymotif import __version__
from . import utils
import sys

from Bio import motifs
from Bio.Seq import Seq
# from Bio.Alphabet import IUPAC
from pyfaidx import Fasta

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
    # mask the genome False by default
    parser.add_argument("-mask", help="Use to mask the genome", action='store_true')
    parser.add_argument("--version", help="Print the version and quit", \
        action="version", version = f'{__version__}')

    # Parse the arguments
    args = parser.parse_args()

    print('PyMotif!')

    # Print the arguments
    print(f"Peak file: {args.peak_file}")
    print(f"Genome: {args.genome}")
    print(f"Output directory: {args.output_directory}")
    print(f"Motif size: {args.size}")
    print(f"Mask genome: {args.mask}")

    # motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_meme.txt"

    motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
    motif_db_file = utils.convert_to_absolute(motif_db_file)
    print(f"Motif database file: {motif_db_file}")

    peak_file = utils.convert_to_absolute(args.peak_file)
    genome_file = utils.convert_to_absolute(args.genome)
    output_directory = utils.convert_to_absolute(args.output_directory)

    peaks_df = utils.read_peaks_file(peak_file)

    # print(peaks_df.head())
    peaks = [tuple(x) for x in peaks_df[['chr', 'start', 'end']].values]

    print(peaks)

    with open(motif_db_file, encoding='utf-8') as f:
        motifs_db = motifs.parse(f, 'jaspar')


    # # Read your motifs from a MEME file (replace with the path to your .meme file)
    # with open(motif_db_file) as f:
    #     for format_type in ['pfm', 'transfac', 'jaspar','meme', 'sites']:
    #     # for format_type in ['minimal', 'pfm', 'transfac', 'jaspar', 'meme', 'sites']:
    #         try:
    #             f.seek(0)
    #             motifs_db = motifs.parse(f, format_type)
    #             print(f"Successfully parsed with format '{format_type}'")
    #         except Exception as e:
    #             print(f"Failed to parse with format '{format_type}': {str(e)}")

    # for motif in motifs_db:
    #     print(motif)
    #     print(dir(motif))
        
    #     print("Motif ID:", motif.base_id)
    #     print("Motif Name:", motif.name)
    #     print("Length of Motif:", motif.length)
    #     num_instances = 0
    #     if motif.instances:
    #         num_instances = len(motif.instances)
    #     print("Number of instances:", num_instances)
    #     print("Consensus sequence:", motif.consensus)
    #     break
    
    # print('motifs_db')
    # print(motifs_db)

    # Load the genome
    genome = Fasta(genome_file)
    print('genome')
    print(genome)

    sequences = utils.get_peak_sequences(peaks, genome)
    print('Number of sequences:', len(sequences))
    # print(sequences)

    for sequence in sequences:
        print('sequence:', sequence)
        found_motifs = utils.find_motifs_in_sequence(sequence, motifs_db)
        print(found_motifs)



    # # Then for each peak sequence in your genome, call the function
    # peak_sequence = Seq("GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA")
    # found_motifs = utils.find_motifs_in_sequence(peak_sequence, motifs_db)
    
    
    
    sys.exit(0)


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