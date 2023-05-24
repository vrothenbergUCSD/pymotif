#!/usr/bin/env python

"""
Command-line script to perform motif finding

Similar to HOMER
"""

import os
import sys
import math
import argparse

from pymotif import __version__
from . import utils

from Bio import motifs
from Bio.Seq import Seq
from pyfaidx import Fasta

from collections import defaultdict


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

    # print(peaks)

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

    for motif in motifs_db:
        print(motif)
        print(dir(motif))
        
        print("Motif ID:", motif.base_id)
        print("Motif Name:", motif.name)
        print("Length of Motif:", motif.length)
        num_instances = 0
        if motif.instances:
            num_instances = len(motif.instances)
        print("Number of instances:", num_instances)
        print("Consensus sequence:", motif.consensus, type(motif.consensus))
        break
    
    # print('motifs_db')
    # print(motifs_db)


    # Load the genome
    genome = Fasta(genome_file)
    print('genome')
    print(genome)

    sequences = utils.get_peak_sequences(peaks, genome)
    print('Number of sequences:', len(sequences))

    all_motifs = defaultdict(list)

    import time 
    start_time = time.time()

    i = 0
    for tup, sequence in sequences.items():
        if i % 100 == 0:
            print(i, 'elapsed time:', time.time() - start_time)
        i += 1
        if i == 100:
            break
        found_motifs = utils.find_motifs_in_sequence(str(sequence), motifs_db)
        # print(found_motifs)
        for motif in found_motifs:
            all_motifs[motif].append([tup, sequence])
            # utils.display_motif(motif)
    
    end_time = time.time()
    print('elapsed time:', end_time - start_time)

    for motif in all_motifs:
        print(motif)
        print(all_motifs[motif])
        break

    sorted_dict = dict(sorted(all_motifs.items(), key=lambda x: len(x[1]), reverse=True))

    from itertools import islice
    top_5_keys = [(motif.base_id, len(all_motifs[motif])) for motif in list(islice(sorted_dict.keys(), 5))]

    print(top_5_keys)   

    

    
    # # Use the function like this:
    # utils.write_known_results_file(found_motifs, 'knownResults.txt', 
    #                                len(sequences), len(background_sequences))
    
    sys.exit(0)
    

if __name__ == "__main__":
    main()