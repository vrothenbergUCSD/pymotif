#!/usr/bin/env python

"""
Command-line script to perform motif finding

Similar to HOMER
"""

import os
import sys
# import math
import argparse

# Saving and loading data
import pickle

# Logging
import logging
import datetime

from collections import defaultdict

from Bio import motifs
from pyfaidx import Fasta


from pymotif import __version__
from . import utils

def print_log(logger, message):
    """
    Print a message to the console and log file
    """
    print(message)
    logger.info(message)

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
    parser.add_argument("-motifs", help="Motifs database to use, e.g. homer or jaspar", type=str, default='jaspar')
    parser.add_argument("-cores", help="Number of cores to use", type=int, default=1)
    parser.add_argument("-enrichment_pval_thresh", help="Enrichment threshold p-value for motif finding", type=float, default=0.00001)
    parser.add_argument("--version", help="Print the version and quit", \
        action="version", version = f'{__version__}')

    logger = logging.getLogger('pymotif_logger')
    # Parse the arguments
    args = parser.parse_args()
    
    # Print the arguments
    print_log(logger, f"PyMotif started")
    print_log(logger, f"PyMotif version: {__version__}")
    print_log(logger, f"Python version: {sys.version}")

    peak_file = utils.convert_to_absolute(args.peak_file)
    print_log(logger, f"Peak file: {peak_file}")
    genome_file = utils.convert_to_absolute(args.genome)
    print_log(logger, f"Genome: {genome_file}")
    output_directory = utils.convert_to_absolute(args.output_directory)
    print_log(logger, f"Output directory: {output_directory}")
    print_log(logger, f"Motif size: {args.size}")
    print_log(logger, f"Mask genome: {args.mask}")
    print_log(logger, f"Motifs: {args.motifs}")
    print_log(logger, f"Number of cores: {args.cores}")
    print_log(logger, f"Enrichment p-value threshold: {args.enrichment_pval_thresh}")

    # Create a log folder if it doesn't exist
    log_folder = os.path.join(output_directory, 'logs')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    # Configure the logger settings
    log_file = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.log'
    log_path = os.path.join(log_folder, log_file)
    print_log(logger, f"Log file path: {log_path}")
    # Hack to remove the default handler
    for handler in logger.root.handlers[:]:
        logger.root.removeHandler(handler)
    logging.basicConfig(
        filename=log_path,  # Specify the filename for the log file
        level=logging.INFO,     # Set the logging level (e.g., INFO, DEBUG)
        format='%(asctime)s %(levelname)s: %(message)s',  # Define the log message format
        datefmt='%Y-%m-%d %H:%M:%S'  # Define the date and time format
    )

    # Check which motifs database to use
    motif_db_file = None
    if args.motifs == 'jaspar':
        motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
    elif args.motifs == 'homer':
        motif_db_file = 'example-files/HOMER_known.motifs'
    motif_db_file = utils.convert_to_absolute(motif_db_file)
    print_log(logger, f"Motif database file: {motif_db_file}")

    with open(motif_db_file, encoding='utf-8') as f:
        if args.motifs == 'homer':
            motifs_db = motifs.parse(f, 'pfm-four-columns')
        elif args.motifs == 'jaspar':
            motifs_db = motifs.parse(f, 'jaspar')
        else:
            # Not tested
            for file_format in ['pfm-four-columns', 'minimal', 'pfm', 'transfac', 'jaspar', 'meme', 'sites']:
                try:
                    f.seek(0)
                    motifs_db = motifs.parse(f, file_format)
                    print(f"Successfully parsed motifs file with format '{file_format}'")
                    logger.info(f"Successfully parsed motifs file with format '{file_format}'")
                    break
                except Exception as exception:
                    print_log(logger, f"Failed to parse with format '{file_format}': {str(exception)}")
                    continue

    peaks_df, genome = utils.read_files(peak_file)
    # print(peaks_df.head())
    # sys.exit(0)
    peaks = [tuple(x) for x in peaks_df[['chr', 'start', 'end']].values]

    # Reducing the number of peaks and motifs for testing
    n_peaks = -1 # All peaks use -1
    n_motifs = -1 # All motifs use -1
    df_first_n = peaks_df.head(n_peaks)
    peaks = [tuple(x) for x in df_first_n[['chr', 'start', 'end']].values]
    motifs_db = motifs_db[:n_motifs]

    # for motif in motifs_db:
    #     print('Motif:')
    #     print(motif)
    #     print(dir(motif))
        
    #     # print("Motif ID:", motif.base_id)
    #     print("Motif Name:", motif.name)
    #     # motif_data = motif.name.split('\t')
    #     # print('motif_data', motif_data)
    #     # motif.name = motif_data[1].split('/')[0].split('(')[0]
    #     # motif.log_odds = float(motif_data[2])
    #     print("Length of Motif:", motif.length)
    #     num_instances = 0
    #     if motif.instances:
    #         num_instances = len(motif.instances)
    #     print("Number of instances:", num_instances)
    #     print("Consensus sequence:", motif.consensus, type(motif.consensus))
    #     print("PWM:", motif.pwm, type(motif.pwm), dir(motif.pwm))
    #     break



    # Load the genome
    genome = Fasta(genome_file)
    print('genome keys:')
    print(genome.keys())

    sequences = utils.get_peak_sequences(peaks, genome)
    print_log(logger, f"Number of sequences: {len(sequences)}")
    peak_seqs = [str(sequence) for sequence in sequences.values()]

    print('peak_seqs', len(peak_seqs), len(peak_seqs[0]))
    # logger.info(f"Length of sequences: {len(peak_seqs[0])}")
    print('peak_seqs[0]:', peak_seqs[0])

    bg_seqs = utils.generate_background_sequences(genome_file, len(peak_seqs), len(peak_seqs[0]))
    print('bg_seqs', len(bg_seqs), len(bg_seqs[0]))
    print('bg_seqs[0]:', bg_seqs[0])

    # DEBUG: Load found motifs from pickle file if it exists
    filename = 'example-files/found_motifs.pkl'
    # Set to True if you want to skip processing and load the pickle file
    load_pickle = False
    if os.path.exists(filename) and load_pickle:
        # Reload the dictionary from the file
        with open(filename, 'rb') as file:
            found_motifs = pickle.load(file)
    else:
        found_motifs = utils.compute_motif_pvals(motifs_db, 
                                                 peak_seqs, 
                                                 bg_seqs, 
                                                 num_sim=3000,
                                                 num_cores=args.cores,
                                                 enrichment_pval_thresh=args.enrichment_pval_thresh,
                                                 )
        # Save the dictionary to a file
        with open(filename, 'wb') as file:
            pickle.dump(found_motifs, file)

    # # Sort the defaultdict based on pval
    # sorted_found_motifs = defaultdict(dict, sorted(found_motifs.items(), key=lambda x: x[1]['pval']))

    # # Print the sorted dictionary
    # print_log(logger, "Top motifs:")
    # for key, value in sorted_found_motifs.items():
    #     print_log(logger, f"Motif: {value['motif'].name}, pval: {value['pval']}, consensus: {value['motif'].consensus}")
    
    # Write the results to a file
    out_file = os.path.join(output_directory, 'knownResults.txt')
    utils.write_known_results_file(found_motifs, 
                                   out_file, 
                                   len(sequences), 
                                   len(bg_seqs), 
                                   top_n=25)

    sys.exit(0)
    

if __name__ == "__main__":
    main()