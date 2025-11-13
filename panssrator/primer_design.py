# panssrator/primer_design.py
import primer3
from panssrator import config, utils

def design_primers_for_ssr(ssr_record: dict, genome_seq: str, flank: int = config.FLANK_SIZE,
                           custom_params: dict = None) -> dict:
    """
    Extract flanking sequence around an SSR locus and design primers using Primer3.

    Parameters:
      ssr_record: Dictionary with SSR information (must include 'start' and 'end').
                  Positions should be 1-indexed (biological convention).
      genome_seq: The entire sequence (e.g., contig or chromosome) in which the SSR is located.
      flank: Number of bases to extract upstream and downstream.
      custom_params: Optional dictionary of primer design parameters (overrides defaults).

    Returns:
      A dictionary of primer design results from Primer3.
    """
    # Validate input
    if "start" not in ssr_record or "end" not in ssr_record:
        utils.logger.error("SSR record missing 'start' or 'end' coordinates")
        return {}

    # Convert 1-indexed SSR positions to 0-indexed for Python slicing
    ssr_start_0 = ssr_record["start"] - 1
    ssr_end_0 = ssr_record["end"]  # end is inclusive in 1-indexed, so this becomes exclusive in 0-indexed

    # Calculate extraction boundaries with bounds checking
    extract_start = max(0, ssr_start_0 - flank)
    extract_end = min(len(genome_seq), ssr_end_0 + flank)

    # Extract template sequence including flanking regions
    template_seq = genome_seq[extract_start:extract_end]

    if not template_seq:
        utils.logger.error(f"Failed to extract template sequence for SSR at {ssr_record['start']}-{ssr_record['end']}")
        return {}

    # Calculate the target region position within the extracted template
    # target_start is where the SSR begins in the template (0-indexed relative to template)
    target_start = ssr_start_0 - extract_start
    target_length = ssr_end_0 - ssr_start_0

    # Validate target region
    if target_start < 0 or target_start + target_length > len(template_seq):
        utils.logger.error(f"Invalid target region calculation for SSR at {ssr_record['start']}-{ssr_record['end']}")
        return {}

    seq_args = {
        'SEQUENCE_ID': f"SSR_{ssr_record['start']}_{ssr_record['end']}",
        'SEQUENCE_TEMPLATE': template_seq,
        'SEQUENCE_TARGET': [target_start, target_length]
    }

    params = custom_params if custom_params else config.PRIMER_PARAMS

    try:
        result = primer3.designPrimers(seq_args, params)
        if not result:
            utils.logger.warning(f"Primer3 returned empty result for SSR at {ssr_record['start']}-{ssr_record['end']}")
    except Exception as e:
        utils.logger.error(f"Primer3 design failed for SSR at {ssr_record['start']}-{ssr_record['end']}: {e}")
        result = {}

    return result

if __name__ == '__main__':
    # Example: design primers for a dummy SSR in a synthetic genome sequence.
    dummy_seq = "N" * config.FLANK_SIZE + "AT" * 10 + "N" * config.FLANK_SIZE
    ssr = {"start": config.FLANK_SIZE + 1, "end": config.FLANK_SIZE + 20, "motif": "AT"}
    primers = design_primers_for_ssr(ssr, dummy_seq)
    utils.logger.info("Primer design results: %s", primers)

