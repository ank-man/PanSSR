# panssr/annotator.py
import os
from typing import List, Dict, Optional
from panssr import utils

try:
    from intervaltree import Interval, IntervalTree
except ImportError:  # fallback for restricted/offline environments
    from panssr._compat.intervaltree import Interval, IntervalTree

def load_annotation(annotation_file: str) -> Dict[str, IntervalTree]:
    """
    Load a GFF/GTF file and return a dictionary mapping chromosomes to an IntervalTree of features.
    
    Each feature is stored as a dictionary with keys: 'type', 'start', 'end', 'strand', 'attributes'.
    """
    trees = {}
    with open(annotation_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            chrom, source, ftype, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)
            feature = {
                "type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "attributes": attributes
            }
            if chrom not in trees:
                trees[chrom] = IntervalTree()
            trees[chrom][start:end+1] = feature  # end+1 because intervaltree is half-open [start, end)
    return trees

# Priority ranking for genomic features (higher priority = more specific annotation)
FEATURE_PRIORITY = {
    # Coding regions (highest priority)
    "cds": 100,
    "CDS": 100,
    "exon": 90,
    "five_prime_utr": 85,
    "five_prime_UTR": 85,
    "5'UTR": 85,
    "three_prime_utr": 85,
    "three_prime_UTR": 85,
    "3'UTR": 85,
    # RNA features
    "mrna": 80,
    "mRNA": 80,
    "transcript": 80,
    "lnc_rna": 75,
    "lncRNA": 75,
    "ncrna": 70,
    "ncRNA": 70,
    "rrna": 70,
    "rRNA": 70,
    "trna": 70,
    "tRNA": 70,
    # Intronic regions
    "intron": 60,
    # Gene regions
    "gene": 50,
    # Regulatory regions
    "promoter": 45,
    "enhancer": 45,
    "regulatory_region": 40,
    # Repeat regions (lower priority - often not informative for SSR markers)
    "repeat_region": 20,
    "tandem_repeat": 15,
    "low_complexity": 10,
    # Intergenic (default)
    "intergenic": 5,
    "intergenic_region": 5,
}

def get_feature_priority(feature_type: str) -> int:
    """
    Get the priority score for a feature type.

    Parameters:
      feature_type: The feature type string

    Returns:
      Priority score (higher = more specific/important)
    """
    # Try exact match first
    if feature_type in FEATURE_PRIORITY:
        return FEATURE_PRIORITY[feature_type]

    # Try case-insensitive match
    feature_lower = feature_type.lower()
    for key, priority in FEATURE_PRIORITY.items():
        if key.lower() == feature_lower:
            return priority

    # Default priority for unknown feature types
    return 30

def annotate_ssr(ssr_record: dict, annot_trees: Dict[str, IntervalTree]) -> Optional[dict]:
    """
    Given an SSR record and annotation interval trees, find the best overlapping feature
    based on a priority ranking system.

    Priority order (highest to lowest):
      CDS/exon > UTR > mRNA/transcript > ncRNA > intron > gene > regulatory > repeats > intergenic

    Parameters:
      ssr_record: SSR dictionary with 'chrom', 'start', 'end'
      annot_trees: Dictionary of IntervalTree objects per chromosome

    Returns:
      The highest-priority overlapping feature dictionary if found; otherwise, None.
    """
    chrom = ssr_record.get("chrom", None)
    if not chrom or chrom not in annot_trees:
        return None

    overlaps = annot_trees[chrom].overlap(ssr_record["start"], ssr_record["end"]+1)

    if not overlaps:
        return None

    # If only one feature overlaps, return it
    if len(overlaps) == 1:
        return list(overlaps)[0].data

    # Multiple features overlap - select the highest priority one
    best_feature = None
    best_priority = -1

    for interval in overlaps:
        feature = interval.data
        feature_type = feature.get("type", "unknown")
        priority = get_feature_priority(feature_type)

        if priority > best_priority:
            best_priority = priority
            best_feature = feature

    return best_feature

if __name__ == '__main__':
    # Example: assume an annotation file "example.gff" exists.
    try:
        trees = load_annotation("example.gff")
        test_ssr = {"chrom": "chr1", "start": 100, "end": 140, "motif": "AT"}
        feature = annotate_ssr(test_ssr, trees)
        if feature:
            utils.logger.info("SSR is located in feature: %s", feature)
        else:
            utils.logger.info("SSR not located in any annotated feature.")
    except Exception as e:
        utils.do_error(str(e))

