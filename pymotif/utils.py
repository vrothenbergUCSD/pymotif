"""
Utility functions for pymotif.
"""

import os
import sys
import pandas as pd
from Bio import motifs
from Bio.Seq import Seq
import math
from scipy.stats import bernoulli


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
    """
    Print an error message and die

    Parameters
    ----------
    msg : str
       Error message to print
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
    sys.exit(1)

def convert_to_absolute(path):
    """
    Check if path is relative and if true then convert to absolute path
    """
    # Expand ~ if present in the path
    path = os.path.expanduser(path)

    if not os.path.isabs(path):
        # Convert relative path to absolute
        path = os.path.abspath(path)
    return path

def read_peaks_file(file_path):
    """
    Reads the peaks file and returns a DataFrame
    """

    # Count the number of header lines
    header_lines = count_header_lines(file_path)

    # Read the file, assuming whitespace-separated values
    # Skip the header lines of metadata
    df = pd.read_csv(file_path, sep="\t", skiprows=header_lines)

    return df

def count_header_lines(file_path):
    """
    Counts the number of lines in the file before the line starting with "#PeakID"
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if line.startswith("#PeakID"):
                return i

    # If no line starts with "#PeakID", return None
    return None

def find_motifs_in_sequence(sequence, motifs_db):
    """
    This function takes a sequence and a motifs database,
    and returns motifs found in the sequence.
    """
    found_motifs = []

    for motif in motifs_db:
        if str(motif.consensus) in sequence:
            found_motifs.append(motif)

    return found_motifs

def get_peak_sequences(peaks, genome):
    """
    This function extracts the sequences corresponding to each peak from the genome.
    :param peaks: List of tuples. Each tuple should have three elements: (chr, start, end).
    :param genome: Fasta object of the genome.
    :return: Dictionary where keys are the peak tuples and values are the sequences.
    """
    sequences = {}
    for peak in peaks:
        chr, start, end = peak
        chr = str(chr)  # convert integer to string
        sequences[peak] = str(genome[chr][start:end].seq)
    return sequences

def write_known_results_file(found_motifs, filename, total_sequences, total_background):
    with open(filename, 'w') as f:
        # Write the header
        f.write("Motif Name\tConsensus\tP-value\tLog P-value\tq-value (Benjamini)\t# of Target Sequences with Motif(of {0})\t% of Target Sequences with Motif\t# of Background Sequences with Motif(of {1})\t% of Background Sequences with Motif\n".format(total_sequences, total_background))
        
        # Write the values for each motif
        for motif in found_motifs:
            name = motif.name
            consensus = motif.consensus
            # Assuming that p_value, q_value are attributes of the motif object (you'll need to calculate them somehow)
            p_value = motif.p_value  
            q_value = motif.q_value 
            log_p_value = math.log10(p_value) * -1  # convert p-value to log scale and make it positive

            # Calculate the number and percentage of sequences with the motif
            num_target_sequences = len(motif.instances)
            percent_target_sequences = num_target_sequences / total_sequences * 100

            # Calculate the number and percentage of background sequences with the motif
            # (This requires knowing how the motif was found in the background sequences)
            num_background_sequences = len(motif.background_instances)
            percent_background_sequences = num_background_sequences / total_background * 100

            # Write the line for this motif
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6:.2f}%\t{7}\t{8:.2f}%\n".format(
                name, consensus, p_value, log_p_value, q_value, 
                num_target_sequences, percent_target_sequences, 
                num_background_sequences, percent_background_sequences))

