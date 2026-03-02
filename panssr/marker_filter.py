# panssr/marker_filter.py
from typing import List, Dict, Any
from panssr import utils

def is_valid_ssr_polymorphism(sizes: List[int], motif: str, max_alleles: int = 10) -> bool:
    """
    Check if amplicon size differences are consistent with SSR repeat variation.

    Parameters:
      sizes: List of amplicon sizes
      motif: SSR motif sequence
      max_alleles: Maximum number of different alleles to allow (to filter non-specific amplification)

    Returns:
      True if the size variation is consistent with SSR polymorphism, False otherwise
    """
    if not sizes or not motif:
        return False

    unique_sizes = list(set(sizes))

    # Check if polymorphic (at least 2 different sizes)
    if len(unique_sizes) < 2:
        return False

    # Filter out markers with too many alleles (likely non-specific)
    if len(unique_sizes) > max_alleles:
        utils.logger.debug(f"Marker has too many alleles ({len(unique_sizes)}), likely non-specific")
        return False

    # Check if size differences are multiples of motif length
    motif_len = len(motif)
    unique_sizes_sorted = sorted(unique_sizes)
    min_size = unique_sizes_sorted[0]

    for size in unique_sizes_sorted[1:]:
        size_diff = size - min_size
        if size_diff % motif_len != 0:
            utils.logger.debug(f"Size difference {size_diff} is not a multiple of motif length {motif_len}")
            return False

    return True

def has_successful_primers(marker: Dict[str, Any]) -> bool:
    """
    Check if the marker has successfully designed primers.

    Parameters:
      marker: Marker dictionary

    Returns:
      True if primers were successfully designed
    """
    primers = marker.get("primers", {})
    if not primers:
        return False

    # Check for at least one primer pair (forward and reverse)
    has_forward = "PRIMER_LEFT_0_SEQUENCE" in primers
    has_reverse = "PRIMER_RIGHT_0_SEQUENCE" in primers

    return has_forward and has_reverse

def is_preferred_annotation(annotation: Any, prefer_intergenic: bool = False) -> bool:
    """
    Check if the marker has preferred annotation.

    Parameters:
      annotation: Annotation information (dict or None)
      prefer_intergenic: If True, prefer intergenic markers

    Returns:
      True if annotation is preferred
    """
    if annotation is None:
        # No annotation = intergenic
        return prefer_intergenic

    if isinstance(annotation, dict):
        feature_type = annotation.get("type", "").lower()
        # Deprioritize repeat regions and low-complexity regions
        if feature_type in ["repeat_region", "low_complexity", "tandem_repeat"]:
            return False
        # If prefer_intergenic is True, only accept intergenic
        if prefer_intergenic and feature_type not in ["intergenic", "intergenic_region"]:
            return False

    return True

def filter_markers(markers: List[Dict[str, Any]],
                   min_amplicon_size: int = 100,
                   max_amplicon_size: int = 500,
                   prefer_intergenic: bool = False,
                   require_primers: bool = True) -> List[Dict[str, Any]]:
    """
    Filter markers based on multiple quality criteria:
      - Must have successful primer design (if require_primers=True)
      - Must be polymorphic (at least 2 different amplicon sizes)
      - Amplicon length differences must be multiples of motif length
      - Amplicon sizes within specified range
      - Not too many alleles (filters non-specific amplification)
      - Optional: prefer intergenic markers

    Parameters:
      markers: List of marker dictionaries
      min_amplicon_size: Minimum acceptable amplicon size (bp)
      max_amplicon_size: Maximum acceptable amplicon size (bp)
      prefer_intergenic: If True, prioritize intergenic markers
      require_primers: If True, only keep markers with successful primer design

    Returns:
      A list of filtered markers.
    """
    filtered = []
    total = len(markers)
    no_sizes = 0
    no_primers = 0
    non_polymorphic = 0
    invalid_polymorphism = 0
    size_out_of_range = 0
    annotation_filtered = 0

    for marker in markers:
        # Check for amplicon sizes from ePCR
        sizes = marker.get("amplicon_sizes", [])
        if not sizes:
            no_sizes += 1
            continue

        # Check primer design success
        if require_primers and not has_successful_primers(marker):
            no_primers += 1
            continue

        # Check amplicon size range
        if any(s < min_amplicon_size or s > max_amplicon_size for s in sizes):
            size_out_of_range += 1
            continue

        # Validate SSR polymorphism
        motif = marker.get("motif", "")
        if not is_valid_ssr_polymorphism(sizes, motif):
            if len(set(sizes)) < 2:
                non_polymorphic += 1
            else:
                invalid_polymorphism += 1
            continue

        # Optional: annotation-based filtering
        if prefer_intergenic:
            annotation = marker.get("annotation")
            if not is_preferred_annotation(annotation, prefer_intergenic=True):
                annotation_filtered += 1
                continue

        filtered.append(marker)

    # Log filtering statistics
    utils.logger.info(f"Marker filtering statistics:")
    utils.logger.info(f"  Total markers: {total}")
    utils.logger.info(f"  Filtered out - no amplicon sizes: {no_sizes}")
    if require_primers:
        utils.logger.info(f"  Filtered out - no primers: {no_primers}")
    utils.logger.info(f"  Filtered out - non-polymorphic: {non_polymorphic}")
    utils.logger.info(f"  Filtered out - invalid polymorphism: {invalid_polymorphism}")
    utils.logger.info(f"  Filtered out - size out of range: {size_out_of_range}")
    if prefer_intergenic:
        utils.logger.info(f"  Filtered out - annotation preference: {annotation_filtered}")
    utils.logger.info(f"  Markers passing filters: {len(filtered)}")

    return filtered

if __name__ == '__main__':
    # Test with dummy markers.
    dummy_markers = [
        {"start": 100, "end": 140, "motif": "AT", "amplicon_sizes": [200, 200, 205]},
        {"start": 150, "end": 190, "motif": "CG", "amplicon_sizes": [300, 300, 300]},
    ]
    filtered = filter_markers(dummy_markers)
    utils.logger.info("Filtered markers: %s", filtered)

