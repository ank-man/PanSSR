"""Utilities to build cross-genome marker presence/polymorphism matrices."""

from __future__ import annotations

from collections import defaultdict


def build_marker_matrix(markers):
    """
    Build per-locus matrix across genomes.

    Marker key is (chrom, start, end, motif). Each row tracks:
    - genome_support: number of genomes where marker observed
    - polymorphic: repeat_count varies across supporting genomes
    - per-genome repeat_count values
    """
    loci = {}
    genome_ids = sorted({m.get("genome_id", "unknown") for m in markers})

    for rec in markers:
        key = (rec.get("chrom"), rec.get("start"), rec.get("end"), rec.get("motif"))
        locus = loci.setdefault(
            key,
            {
                "chrom": rec.get("chrom"),
                "start": rec.get("start"),
                "end": rec.get("end"),
                "motif": rec.get("motif"),
                "annotation": rec.get("annotation"),
                "genome_values": defaultdict(list),
            },
        )
        gid = rec.get("genome_id", "unknown")
        locus["genome_values"][gid].append(rec.get("repeat_count"))

    rows = []
    for locus in loci.values():
        genome_values = {}
        repeats_flat = []
        for gid in genome_ids:
            values = [v for v in locus["genome_values"].get(gid, []) if isinstance(v, int)]
            if values:
                # representative value per genome (first found)
                genome_values[gid] = values[0]
                repeats_flat.extend(values)
            else:
                genome_values[gid] = "NA"

        unique_repeat_counts = sorted({v for v in repeats_flat})
        rows.append(
            {
                "chrom": locus["chrom"],
                "start": locus["start"],
                "end": locus["end"],
                "motif": locus["motif"],
                "annotation": locus["annotation"],
                "genome_support": sum(1 for v in genome_values.values() if v != "NA"),
                "unique_repeat_counts": unique_repeat_counts,
                "polymorphic": len(unique_repeat_counts) >= 2,
                "genome_values": genome_values,
            }
        )

    rows.sort(key=lambda r: (str(r["chrom"]), int(r["start"]), int(r["end"]), str(r["motif"])))
    return genome_ids, rows


def select_common_polymorphic(rows, min_genome_support):
    """Select loci that are polymorphic and present in at least min_genome_support genomes."""
    return [r for r in rows if r["polymorphic"] and r["genome_support"] >= min_genome_support]
