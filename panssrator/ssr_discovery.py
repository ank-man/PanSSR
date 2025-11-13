# panssrator/ssr_discovery.py
import re
from typing import List, Dict, Any
from panssrator import config, utils

@utils.timeit
def detect_ssrs(sequence: str, min_repeats: Dict[str, int] = None) -> List[Dict[str, Any]]:
    """
    Scan a DNA sequence for perfect SSRs.
    
    Parameters:
      sequence: The DNA sequence to scan.
      min_repeats: Dictionary specifying the minimum number of repeats for each motif size.
                   Defaults to config.DEFAULT_MIN_REPEATS.
    
    Returns:
      A list of SSR records. Each record is a dictionary with keys:
        - 'motif': the repeated unit (e.g., 'AT')
        - 'start': start position (1-indexed)
        - 'end': end position
        - 'repeat_count': number of repeats
        - 'raw_sequence': the matched repeat region
    """
    if min_repeats is None:
        min_repeats = config.DEFAULT_MIN_REPEATS
    ssr_records = []
    sequence = sequence.upper()

    # Loop over motif lengths 1 to 6
    for motif_length in range(1, 7):
        # Build a regex pattern to capture perfect repeats.
        # e.g. for motif_length=2, pattern becomes: (..)\1{min_repeats-1,}
        min_rep = min_repeats.get({1:"mono", 2:"di", 3:"tri", 4:"tetra", 5:"penta", 6:"hexa"}[motif_length], 0)
        if min_rep == 0:
            continue  # skip if no threshold defined
        pattern = f"((.{{{motif_length}}}))\\1{{{min_rep - 1},}}"
        for match in re.finditer(pattern, sequence):
            raw = match.group(0)
            repeat_unit = match.group(2)
            repeat_count = len(raw) // motif_length
            # Optionally filter if the repeat region is too long
            if len(raw) > config.MAX_SSR_LENGTH:
                continue
            ssr_records.append({
                "motif": repeat_unit,
                "start": match.start() + 1,  # convert to 1-indexed
                "end": match.end(),
                "repeat_count": repeat_count,
                "raw_sequence": raw
            })
    return ssr_records

if __name__ == '__main__':
    # Example test
    test_seq = "ATATATATATCGCGCGCGATATAT"
    records = detect_ssrs(test_seq)
    for rec in records:
        utils.logger.info("%s", rec)

