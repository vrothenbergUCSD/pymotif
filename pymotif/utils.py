"""
Utility functions for pymotif.
"""

import sys
import pandas as pd

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



