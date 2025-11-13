# panssrator/epcr.py
import re
import sys
from panssrator import config, utils

# Try to import tre module for fuzzy matching
TRE_AVAILABLE = False
try:
    import tre
    TRE_AVAILABLE = True
except ImportError:
    utils.logger.warning("The 'tre' module is not available. Falling back to exact matching for ePCR simulation.")
    utils.logger.warning("For fuzzy matching support, install tre from: https://github.com/laurikari/tre/")

def find_primer_positions_exact(seq: str, primer: str) -> list:
    """
    Find all positions of a primer in a sequence using exact matching.

    Parameters:
      seq: The sequence to search in
      primer: The primer sequence

    Returns:
      List of start positions (0-indexed)
    """
    positions = []
    primer_pattern = replace_ambiguity_codes(primer)
    pattern = re.compile(primer_pattern, re.IGNORECASE)

    pos = 0
    while pos < len(seq):
        match = pattern.search(seq, pos)
        if not match:
            break
        positions.append(match.start())
        pos = match.start() + 1

    return positions

def find_primer_positions_fuzzy(seq: str, primer: str, max_cost: int) -> list:
    """
    Find all positions of a primer in a sequence using fuzzy matching (requires tre module).

    Parameters:
      seq: The sequence to search in
      primer: The primer sequence
      max_cost: Maximum allowed mismatches

    Returns:
      List of start positions (0-indexed)
    """
    if not TRE_AVAILABLE:
        return find_primer_positions_exact(seq, primer)

    positions = []
    primer_pattern = replace_ambiguity_codes(primer)

    try:
        pattern = tre.compile(primer_pattern, tre.EXTENDED)
        pos = 0
        while pos < len(seq):
            match = pattern.search(seq, pos, fuzzyness=max_cost)
            if not match:
                break
            positions.append(match.start())
            pos = match.start() + 1
    except Exception as e:
        utils.logger.error(f"Error in fuzzy matching: {e}")
        return []

    return positions

def simulate_epcr(genome_seq: str, primer_pair: dict, max_cost: int = config.MAX_EPCR_COST) -> list:
    """
    Simulate in silico PCR by searching for primer binding sites in genome_seq.

    The function searches for the forward primer as-is, and searches for the
    reverse complement of the reverse primer (as it would bind to the opposite strand).

    Parameters:
      genome_seq: The target genome sequence.
      primer_pair: Dictionary with keys 'forward' and 'reverse' containing primer sequences.
      max_cost: Maximum allowed fuzzy matching cost (only used if tre module is available).

    Returns:
      A list of predicted amplicon sizes (in bp) where both primers are found in correct orientation.
    """
    forward_primer = primer_pair.get("forward", "")
    reverse_primer = primer_pair.get("reverse", "")

    if not forward_primer or not reverse_primer:
        utils.logger.error("Invalid primer pair: both forward and reverse primers are required")
        return []

    # Convert reverse primer to reverse complement for searching
    # (the reverse primer binds to the opposite strand)
    reverse_primer_rc = utils.reverse_complement(reverse_primer)

    # Find primer binding positions
    if TRE_AVAILABLE and max_cost > 0:
        f_positions = find_primer_positions_fuzzy(genome_seq, forward_primer, max_cost)
        r_positions = find_primer_positions_fuzzy(genome_seq, reverse_primer_rc, max_cost)
    else:
        f_positions = find_primer_positions_exact(genome_seq, forward_primer)
        r_positions = find_primer_positions_exact(genome_seq, reverse_primer_rc)

    if not f_positions or not r_positions:
        return []

    amplicon_sizes = []
    # For each forward match, find reverse matches that are downstream
    for f_pos in f_positions:
        for r_pos in r_positions:
            # Reverse primer position should be after forward primer
            if r_pos > f_pos:
                # Calculate amplicon size: distance + length of reverse primer
                size = r_pos - f_pos + len(reverse_primer)
                if size <= config.MAX_EPCR_PRODUCT:
                    amplicon_sizes.append(size)

    return amplicon_sizes

def replace_ambiguity_codes(seq: str) -> str:
    """Wrapper to call utils.replace_ambiguity_codes."""
    return utils.replace_ambiguity_codes(seq)

if __name__ == '__main__':
    # Example test for ePCR simulation
    test_genome = "N" * 100 + "ATGCGT" + "N" * 50 + "CATGCA" + "N" * 100
    primers = {"forward": "ATGCGT", "reverse": "CATGCA"}
    products = simulate_epcr(test_genome, primers)
    utils.logger.info("Predicted amplicon sizes: %s", products)

