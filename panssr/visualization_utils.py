"""Utilities for summarizing PanSSR output files for GUI visualizations."""

from __future__ import annotations

import csv
from collections import Counter
from statistics import mean


def _to_int(value, default=0):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def load_marker_records(path: str):
    """Load marker TSV rows as dictionaries."""
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def summarize_markers(records):
    """Build marker summary stats and distributions."""
    motif_counter = Counter()
    motif_len_counter = Counter()
    chrom_counter = Counter()
    repeats = []

    for rec in records:
        motif = (rec.get("motif") or "").upper()
        if motif:
            motif_counter[motif] += 1
            motif_len_counter[f"{len(motif)}-mer"] += 1

        chrom = rec.get("chrom") or "unknown"
        chrom_counter[chrom] += 1

        repeats.append(_to_int(rec.get("repeat_count"), 0))

    repeats_nonzero = [r for r in repeats if r > 0]
    return {
        "total_markers": len(records),
        "unique_motifs": len(motif_counter),
        "chromosomes": len(chrom_counter),
        "mean_repeat_count": round(mean(repeats_nonzero), 2) if repeats_nonzero else 0.0,
        "motif_counts": dict(motif_counter.most_common(20)),
        "motif_length_counts": dict(motif_len_counter),
        "chromosome_counts": dict(chrom_counter),
    }


def load_genotype_records(path: str):
    """Load genotype CSV rows as dictionaries."""
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader)


def summarize_genotypes(records):
    """Build genotype summary stats and distributions."""
    bam_counter = Counter()
    genotype_counter = Counter()
    missing = 0

    for rec in records:
        bam_counter[rec.get("BAM_file") or "unknown"] += 1
        genotype = rec.get("Genotype") or "NA"
        genotype_counter[genotype] += 1
        if genotype == "NA":
            missing += 1

    total = len(records)
    call_rate = round((total - missing) / total * 100, 2) if total else 0.0

    return {
        "total_calls": total,
        "samples": len(bam_counter),
        "call_rate": call_rate,
        "missing_calls": missing,
        "sample_counts": dict(bam_counter),
        "genotype_counts": dict(genotype_counter.most_common(20)),
    }
