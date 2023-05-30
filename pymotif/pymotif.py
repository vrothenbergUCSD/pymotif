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

# Performance testing
import cProfile
import pstats

from collections import defaultdict

# Debugging
from itertools import islice

from Bio import motifs
# from Bio.Seq import Seq
from pyfaidx import Fasta


from pymotif import __version__
from . import utils


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        prog="pymotif",
        description="Command-line script to perform motif finding ...",
    )
    print('PyMotif started')

    logger = logging.getLogger('pymotif_logger')

    logger.info("PyMotif started")
    logger.info(f"PyMotif version: {__version__}")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Python executable: {sys.executable}")


    # Create a log folder if it doesn't exist
    log_folder = utils.convert_to_absolute('logs')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

    # Configure the logger settings
    log_file = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.log'
    print('log_file:', log_file)
    log_path = os.path.join(log_folder, log_file)
    print('log_path:', log_path)
    # Hack to remove the default handler
    for handler in logger.root.handlers[:]:
        logger.root.removeHandler(handler)
    logging.basicConfig(
        filename=log_path,  # Specify the filename for the log file
        level=logging.INFO,     # Set the logging level (e.g., INFO, DEBUG)
        format='%(asctime)s %(levelname)s: %(message)s',  # Define the log message format
        datefmt='%Y-%m-%d %H:%M:%S'  # Define the date and time format
    )

    # Get the directory path of the pymotif.py file
    curdir = os.path.abspath(os.path.dirname(__file__))
    print('curdir:', curdir)
    logger.info(f"Current directory: {curdir}")

    # Get the directory path from command-line argument
    tool_dir = os.path.dirname(sys.argv[0])
    print('tool_dir:', tool_dir)
    logger.info(f"Tool directory: {tool_dir}")

    # Get the parent directory
    root_dir = os.path.dirname(tool_dir)
    print('root_dir:', root_dir)
    logger.info(f"Root directory: {root_dir}")

    # Required positional arguments
    parser.add_argument("peak_file", help="Peak or BED file", type=str)
    parser.add_argument("genome", help="Genome fasta file", type=str)
    parser.add_argument("output_directory", help="Directory to output results", type=str)

    # Optional arguments
    parser.add_argument("-size", help="Size of motif to search for", type=int, default=200)
    # mask the genome False by default
    parser.add_argument("-mask", help="Use to mask the genome", action='store_true')
    parser.add_argument("-motifs", help="Motifs database to use, e.g. homer or jaspar", type=str, default='homer')
    parser.add_argument("-cores", help="Number of cores to use", type=int, default=1)
    parser.add_argument("--version", help="Print the version and quit", \
        action="version", version = f'{__version__}')

    # Parse the arguments
    args = parser.parse_args()

    # Print the arguments
    print(f"Peak file: {args.peak_file}")
    logger.info(f"Peak file: {args.peak_file}")
    print(f"Genome: {args.genome}")
    logger.info(f"Genome: {args.genome}")
    print(f"Output directory: {args.output_directory}")
    logger.info(f"Output directory: {args.output_directory}")
    print(f"Motif size: {args.size}")
    logger.info(f"Motif size: {args.size}")
    print(f"Mask genome: {args.mask}")
    logger.info(f"Mask genome: {args.mask}")
    print(f"Motifs: {args.motifs}")
    logger.info(f"Motifs: {args.motifs}")
    print(f"Number of cores: {args.cores}")
    logger.info(f"Number of cores: {args.cores}")

    # motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_meme.txt"
    motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
    if args.motifs == 'jaspar':
        motif_db_file = "example-files/JASPAR2022_CORE_non-redundant_pfms_jaspar.txt"
    if args.motifs == 'homer':
        motif_db_file = 'example-files/HOMER_known.motifs'
    logger.info(f"Motif database file: {motif_db_file}")
    


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
                    logger.info(f"Failed to parse with format '{file_format}': {str(exception)}")
                    continue

    # motif_db_file = utils.convert_to_absolute(motif_db_file)
    # print(f"Motif database file: {motif_db_file}")
    # logger.info(f"Motif database file: {motif_db_file}")

    peak_file = utils.convert_to_absolute(args.peak_file)
    logger.info(f"Peak file: {peak_file}")
    genome_file = utils.convert_to_absolute(args.genome)
    logger.info(f"Genome file: {genome_file}")
    output_directory = utils.convert_to_absolute(args.output_directory)
    logger.info(f"Output directory: {output_directory}")

    peaks_df = utils.read_peaks_file(peak_file)
    peaks = [tuple(x) for x in peaks_df[['chr', 'start', 'end']].values]

    
    # Reducing the number of peaks and motifs for testing
    n = 500
    df_first_n = peaks_df.head(n)
    peaks = [tuple(x) for x in df_first_n[['chr', 'start', 'end']].values]
    motifs_db = motifs_db[:n]

    for motif in motifs_db:
        print('Motif:')
        print(motif)
        print(dir(motif))
        
        # print("Motif ID:", motif.base_id)
        print("Motif Name:", motif.name)
        # motif_data = motif.name.split('\t')
        # print('motif_data', motif_data)
        # motif.name = motif_data[1].split('/')[0].split('(')[0]
        # motif.log_odds = float(motif_data[2])
        print("Length of Motif:", motif.length)
        num_instances = 0
        if motif.instances:
            num_instances = len(motif.instances)
        print("Number of instances:", num_instances)
        print("Consensus sequence:", motif.consensus, type(motif.consensus))
        print("PWM:", motif.pwm, type(motif.pwm), dir(motif.pwm))
        break

    
    # print('motifs_db')
    # print(motifs_db)


    # Load the genome
    genome = Fasta(genome_file)

    sequences = utils.get_peak_sequences(peaks, genome)
    print('Number of sequences:', len(sequences))
    logger.info(f"Number of sequences: {len(sequences)}")
    peak_seqs = [str(sequence) for sequence in sequences.values()]
    print('peak_seqs', len(peak_seqs), len(peak_seqs[0]))
    # logger.info(f"Length of sequences: {len(peak_seqs[0])}")
    print('peak_seqs[0]:', peak_seqs[0])

    bg_seqs = utils.generate_background_sequences(genome_file, len(peak_seqs), len(peak_seqs[0]))
    print('bg_seqs', len(bg_seqs), len(bg_seqs[0]))
    print('bg_seqs[0]:', bg_seqs[0])

    # # Define a dictionary with your variables
    # variables = {'utils': utils, 'motifs_db': motifs_db, 'peak_seqs': peak_seqs, 'bg_seqs': bg_seqs}

    # # Run the profiler and save results to a file
    # found_motifs = cProfile.runctx('utils.compute_motif_pvals(motifs_db, peak_seqs, bg_seqs)', variables, {}, 'stats_outputfile')

    # # Create a pstats.Stats object
    # p = pstats.Stats('stats_outputfile')

    # # Print all statistics
    # p.print_stats()

    # # Sort statistics by cumulative time and print them
    # p.sort_stats('cumulative').print_stats()

    # # Sort statistics by function time and print them
    # p.sort_stats('time').print_stats()

    # Load found motifs from pickle file if it exists
    filename = 'example-files/found_motifs.pkl'
    if os.path.exists(filename) and False:
        # Reload the dictionary from the file
        with open(filename, 'rb') as file:
            found_motifs = pickle.load(file)
    else:
        found_motifs = utils.compute_motif_pvals(motifs_db, 
                                                 peak_seqs, 
                                                 bg_seqs, 
                                                 num_sim=3000,
                                                 num_cores=args.cores
                                                 )
        # Save the dictionary to a file
        with open(filename, 'wb') as file:
            pickle.dump(found_motifs, file)
    print('found_motifs', len(found_motifs))
    # print(found_motifs)

    # Sort the defaultdict based on pval
    sorted_dict = defaultdict(dict, sorted(found_motifs.items(), key=lambda x: x[1]['pval']))

    # # Print the sorted dictionary
    for key, value in sorted_dict.items():
        print(key, value['pval'], value['motif'].consensus)
        print(value['motif'].name, value['motif'].tf_class, value['motif'].tf_family, value['motif'].species)
        # pri
        # nt(dir(value['motif']))
        break

    
    top_25_keys = [key for key in list(islice(sorted_dict.keys(), 25))]

    print('top_25_keys')
    for motif_id in top_25_keys:
        motif = sorted_dict[motif_id]['motif']
        print(motif_id, motif.name, sorted_dict[motif_id]['pval'], motif.consensus)
        





    # sorted_dict = dict(sorted(all_motifs.items(), key=lambda x: len(x[1]), reverse=True))





    
    # # Use the function like this:
    # utils.write_known_results_file(found_motifs, 'knownResults.txt', 
    #                                len(sequences), len(background_sequences))


    
    sys.exit(0)
    

if __name__ == "__main__":
    main()