# panssr/genotyper.py
import pysam
import re
from collections import Counter
from panssr import config, utils

def count_repeat_units(seq: str, motif: str) -> int:
    """
    Count the number of times a given motif is repeated consecutively in a sequence.
    This approach uses regex matching to find the longest consecutive repeat stretch.
    """
    if not seq or not motif:
        return 0
    # Build a regex pattern for the motif repeated consecutively
    pattern = f"(?:{re.escape(motif)})+"
    match = re.search(pattern, seq, re.IGNORECASE)
    if match:
        repeated_seq = match.group(0)
        return len(repeated_seq) // len(motif)
    return 0

def extract_ssr_region_from_read(read, ssr_start: int, ssr_end: int) -> str:
    """
    Extract the portion of a read that overlaps with the SSR region using CIGAR-aware alignment.

    Parameters:
      read: pysam AlignedSegment object
      ssr_start: Start position of SSR in reference (1-indexed)
      ssr_end: End position of SSR in reference (1-indexed, inclusive)

    Returns:
      The extracted sequence from the read overlapping the SSR region, or empty string if no overlap.
    """
    if read.is_unmapped or not read.cigartuples:
        return ""

    # Convert to 0-indexed for internal calculations
    ssr_start_0 = ssr_start - 1
    ssr_end_0 = ssr_end

    # Get alignment positions
    ref_pos = read.reference_start  # 0-indexed position in reference
    query_pos = 0  # Position in read sequence

    extracted_seq = []
    in_ssr_region = False

    for op, length in read.cigartuples:
        # CIGAR operations:
        # 0 = M (match/mismatch)
        # 1 = I (insertion to reference)
        # 2 = D (deletion from reference)
        # 3 = N (skipped region)
        # 4 = S (soft clipping)
        # 5 = H (hard clipping)

        if op == 0:  # Match/mismatch - consumes both query and reference
            for i in range(length):
                if ssr_start_0 <= ref_pos < ssr_end_0:
                    if query_pos < len(read.query_sequence):
                        extracted_seq.append(read.query_sequence[query_pos])
                        in_ssr_region = True
                elif in_ssr_region:
                    # We've passed the SSR region
                    return ''.join(extracted_seq)
                ref_pos += 1
                query_pos += 1

        elif op == 1:  # Insertion - consumes query only
            # Include insertions if we're in the SSR region
            if in_ssr_region or (ssr_start_0 <= ref_pos < ssr_end_0):
                for i in range(length):
                    if query_pos < len(read.query_sequence):
                        extracted_seq.append(read.query_sequence[query_pos])
                    query_pos += 1
            else:
                query_pos += length

        elif op == 2:  # Deletion - consumes reference only
            ref_pos += length

        elif op == 4:  # Soft clipping - consumes query only
            query_pos += length

        elif op == 5:  # Hard clipping - consumes nothing
            pass

        elif op == 3:  # Skipped region (e.g., intron) - consumes reference only
            ref_pos += length

    return ''.join(extracted_seq)



def genotype_marker_with_handle(bam, ssr_record: dict) -> dict:
    """
    Genotype a single SSR marker using an already opened BAM handle.

    Parameters:
      bam: Open pysam.AlignmentFile in read-binary mode.
      ssr_record: Dictionary with keys including 'chrom', 'start', 'end', 'motif'.

    Returns:
      A dictionary containing allele counts and a genotype call.
    """
    allele_counts = Counter()
    chrom = ssr_record.get("chrom", None)
    if not chrom:
        utils.logger.error("SSR record does not contain chromosome information.")
        return {"allele_counts": {}, "genotype": None}

    start = ssr_record["start"]
    end = ssr_record["end"]
    motif = ssr_record.get("motif", "")

    if not motif:
        utils.logger.error("SSR record does not contain motif information.")
        return {"allele_counts": {}, "genotype": None}

    # Fetch reads overlapping the SSR region
    try:
        for read in bam.fetch(chrom, start - 1, end + 1):
            # Filter by mapping quality
            if read.mapping_quality < config.MIN_MAPQ:
                continue

            # Skip unmapped, secondary, and supplementary reads
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            # Extract the SSR region from the read using CIGAR-aware method
            ssr_seq = extract_ssr_region_from_read(read, start, end)

            if ssr_seq:
                # Count repeat units in the extracted sequence
                repeat_count = count_repeat_units(ssr_seq, motif)
                if repeat_count > 0:
                    allele_counts[repeat_count] += 1
    except Exception as e:
        utils.logger.error(f"Error fetching reads for {chrom}:{start}-{end}: {e}")

    # Determine genotype based on allele counts
    total = sum(allele_counts.values())
    if total < config.MIN_READ_SUPPORT:
        genotype = None  # Insufficient read support
    else:
        # For diploid species, expect at most two alleles.
        most_common = allele_counts.most_common(2)
        if len(most_common) == 1:
            genotype = {"alleles": [most_common[0][0]], "support": most_common[0][1]}
        else:
            genotype = {
                "alleles": [most_common[0][0], most_common[1][0]],
                "support": {
                    most_common[0][0]: most_common[0][1],
                    most_common[1][0]: most_common[1][1],
                },
            }
    return {"allele_counts": dict(allele_counts), "genotype": genotype}


def genotype_marker(bam_file: str, ssr_record: dict) -> dict:
    """Compatibility wrapper that opens the BAM file and calls genotype_marker_with_handle."""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        utils.logger.error(f"Failed to open BAM file {bam_file}: {e}")
        return {"allele_counts": {}, "genotype": None}

    try:
        return genotype_marker_with_handle(bam, ssr_record)
    finally:
        bam.close()

if __name__ == '__main__':
    # Example test: Replace 'example.bam' with an actual BAM file to run this test.
    ssr_example = {"chrom": "chr1", "start": 100, "end": 140, "motif": "AT"}
    gt = genotype_marker("example.bam", ssr_example)
    utils.logger.info("Genotype call: %s", gt)

