# panssrator/config.py
import os

# ---------------------------
# SSR Detection Parameters
# ---------------------------
# Minimum repeat numbers for each motif type
DEFAULT_MIN_REPEATS = {
    "mono": 12,
    "di": 7,
    "tri": 5,
    "tetra": 4,
    "penta": 4,
    "hexa": 4,
}

# Maximum SSR length (to avoid very long repeats)
MAX_SSR_LENGTH = 80  # in bp

# Minimum distance between adjacent SSRs (to avoid compound SSRs)
MIN_FLANK_BETWEEN_SSR = 200  # in bp

# ---------------------------
# Primer Design Parameters (for primer3)
# ---------------------------
PRIMER_PARAMS = {
    "PRIMER_PRODUCT_SIZE_RANGE": "100-400",
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_OPT_SIZE": 19,
    "PRIMER_MAX_SIZE": 23,
    "PRIMER_MIN_TM": 52,
    "PRIMER_OPT_TM": 55,
    "PRIMER_MAX_TM": 60,
    "PRIMER_MIN_GC": 40,
    "PRIMER_MAX_GC": 70,
    "PRIMER_MAX_POLY_X": 4,  # maximum mononucleotide repeat in primer
}

# Flanking region to extract for primer design (in bp)
FLANK_SIZE = 100

# ---------------------------
# In Silico PCR (ePCR) Parameters
# ---------------------------
# Maximum cost (i.e. allowed mismatches/edits) when matching primers
MAX_EPCR_COST = 3

# Maximum allowed product length for ePCR simulation (in bp)
MAX_EPCR_PRODUCT = 1500

# ---------------------------
# Genotyping Parameters (BAM processing)
# ---------------------------
# Minimum mapping quality to consider a read
MIN_MAPQ = 45

# Minimum read support for genotype call
MIN_READ_SUPPORT = 3

# ---------------------------
# General Settings
# ---------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATABASE_FILE = os.path.join(BASE_DIR, "panssrator.db")

# You can add more parameters here as needed.

