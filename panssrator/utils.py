# panssrator/utils.py
import os
import sys
import logging
import subprocess
from functools import wraps

# Set up a logger for the application
logger = logging.getLogger("PanSSRAtor")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

def setup_logging(level=logging.INFO):
    """Configure the global logging level."""
    logger.setLevel(level)

def run_command(cmd, cwd=None):
    """Run an external command and log its execution."""
    logger.info("Executing command: %s", " ".join(cmd))
    try:
        subprocess.check_call(cmd, cwd=cwd)
    except subprocess.CalledProcessError as e:
        logger.error("Command failed: %s", e)
        raise

def do_error(message):
    """Log an error and exit."""
    logger.error(message)
    sys.exit(message)

def reverse_complement(seq):
    """Return the reverse complement of the DNA sequence."""
    # Using str.translate for speed
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]

def replace_ambiguity_codes(seq):
    """
    Convert IUPAC ambiguity codes into regex character classes.
    For example, R -> [AG], Y -> [CT], etc.
    """
    codes = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]",
        "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]",
        "N": "[ACGT]", "X": "[ACGT]"
    }
    try:
        return "".join(codes[base] for base in seq.upper())
    except KeyError as e:
        do_error(f"Unrecognized nucleotide code: {e}")

# A decorator for timing functions (for performance debugging)
def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        from time import time
        start = time()
        result = func(*args, **kwargs)
        logger.info("Function %s took %.3f seconds", func.__name__, time() - start)
        return result
    return wrapper

