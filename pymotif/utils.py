"""
Utility functions for pymotif.
"""

import os
import sys
import math
import random
import pandas as pd
import numpy as np
import scipy.stats as stats
from Bio import motifs
from Bio.Seq import Seq


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


def display_motif(motif):
    print("Motif:", motif, type(motif))
    print("Motif ID:", motif.base_id)
    print("Motif Name:", motif.name)
    print("Length of Motif:", motif.length)
    num_instances = 0
    if motif.instances:
        num_instances = len(motif.instances)
    print("Number of instances:", num_instances)
    print("Consensus sequence:", motif.consensus, type(motif.consensus))
    print("PWM:", motif.pwm, type(motif.pwm))
    print(dir(motif))


nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
rev_nucs = {"A": "T", "C": "G", "G": "C", "T": "A"}

def score_seq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    for i,b in enumerate(sequence):
        # Check if sequence nucleotide index is out of bounds
        if i >= pwm.shape[1]:
            break
        ind = nucs[b]
        score += pwm[ind,i]     
    return score

def reverse_complement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    comp = [rev_nucs[s] for s in sequence]
    revcomp = ''.join(comp)[::-1]
    return revcomp

def find_max_score(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
    max_score = -1*np.inf
    rev_sequence = reverse_complement(sequence)
    M = pwm.shape[1]
    L = len(sequence)
    for i in range(L-M+1):      
        fwd_seq = sequence[i:i+M]
        fwd_score = score_seq(pwm, fwd_seq)
        rev_seq = rev_sequence[i:i+M]
        rev_score = score_seq(pwm, rev_seq)      
        max_score = max(max_score, fwd_score, rev_score)

    return max_score

def compute_nuc_freqs(sequences):
    """ Compute nucleotide frequencies of a list of sequences
    
    Parameters
    ----------
    sequences : list of str
       List of sequences
       
    Returns
    -------
    freqs : list of float
       Frequencies of A, C, G, T in the sequences
    """
    freqs = [0.25, 0.25, 0.25, 0.25] # compute frequency of A, C, G, T
    nuc_counts = [0,0,0,0]
    for seq in sequences:
        for nuc in seq:
            ind = nucs[nuc]
            nuc_counts[ind] += 1      
    nuc_total = sum(nuc_counts)
    freqs = nuc_counts / nuc_total
    return freqs

def random_sequence(n, freqs):
    """ Generate a random string of nucleotides of length n
    
    Use the given nucleotide frequences
    
    Parameters
    ----------
    n : int
       Length of random string to generate
    freqs : list of float
       List of frequencies of A, C, G, T
       
    Returns
    -------
    seq : str
       random sequence of length n with the specified allele frequencies
    """
    nucleotides = ['A', 'C', 'G', 'T']
    seq = ''.join(np.random.choice(nucleotides, size=n, p=freqs))
    return seq

def get_threshold(null_dist, pval):
    """ Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       pval% of null_dist should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 # set this  below to be the score threshold to obtain a p-value <0.01
    thresh = np.percentile(null_dist, (1 - pval) * 100)
    return thresh

def compute_enrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    pval = -1
    A = peak_motif
    B = peak_total - peak_motif
    C = bg_motif
    D = bg_total - bg_motif
    table = [[A,B],[C,D]]
    odds, pval = stats.fisher_exact(table)
    return pval

def compute_motif_pvals(PWMList, peak_seqs, bg_seqs, pwm_thresholds, pwm_names):
    for i in range(len(PWMList)):
        pwm = PWMList[i]
        thresh = pwm_thresholds[i]
        num_peak_pass = np.sum([int(find_max_score(pwm, seq)>thresh) for seq in peak_seqs])
        num_bg_pass = np.sum([int(find_max_score(pwm, seq)>thresh) for seq in bg_seqs])
        pval = compute_enrichment(len(peak_seqs), num_peak_pass, len(bg_seqs), num_bg_pass)
        print("PWM: %s, %s/%s peaks, %s/%s background; p-val: %s"%(pwm_names[i], \
            num_peak_pass, len(peak_seqs), num_bg_pass, len(bg_seqs), pval))
        