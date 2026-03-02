"""Helper utilities for genotype-mode orchestration and output formatting."""

from panssr import utils


def filter_markers_by_reference(markers, ref_seqs):
    """Keep markers whose chromosome exists in reference; log missing chromosomes once."""
    missing = set()
    filtered = []
    for marker in markers:
        chrom = marker.get("chrom")
        if chrom in ref_seqs:
            filtered.append(marker)
        else:
            missing.add(chrom)
    for chrom in sorted(missing):
        utils.logger.warning(f"Chromosome {chrom} not found in reference genome; skipping related markers")
    return filtered


def format_genotype_call(call: dict):
    """Format genotype call dictionary into CSV-friendly genotype and allele-count strings."""
    genotype = call.get("genotype", None)
    allele_counts = call.get("allele_counts", {})

    if genotype:
        gt_str = "/".join(map(str, genotype.get("alleles", [])))
    else:
        gt_str = "NA"

    allele_str = ";".join(f"{allele}:{count}" for allele, count in sorted(allele_counts.items()))
    if not allele_str:
        allele_str = "NA"

    return gt_str, allele_str
