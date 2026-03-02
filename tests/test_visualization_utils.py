"""Tests for visualization utility summaries."""

from panssr.visualization_utils import (
    load_genotype_records,
    load_marker_records,
    summarize_genotypes,
    summarize_markers,
)


def test_marker_summary(tmp_path):
    marker_file = tmp_path / "markers.tsv"
    marker_file.write_text(
        "chrom\tstart\tend\tmotif\trepeat_count\n"
        "chr1\t1\t10\tAT\t5\n"
        "chr1\t12\t18\tA\t12\n"
        "chr2\t30\t40\tGATA\t4\n"
    )
    records = load_marker_records(str(marker_file))
    summary = summarize_markers(records)
    assert summary["total_markers"] == 3
    assert summary["unique_motifs"] == 3
    assert summary["chromosome_counts"]["chr1"] == 2
    assert summary["motif_length_counts"]["2-mer"] == 1


def test_genotype_summary(tmp_path):
    gt_file = tmp_path / "genotypes.csv"
    gt_file.write_text(
        "BAM_file,Marker_chrom,Marker_start,Marker_end,Genotype,Allele_counts\n"
        "s1.bam,chr1,1,10,10/12,10:3;12:2\n"
        "s1.bam,chr1,20,30,NA,NA\n"
        "s2.bam,chr2,5,14,8/8,8:5\n"
    )
    records = load_genotype_records(str(gt_file))
    summary = summarize_genotypes(records)
    assert summary["total_calls"] == 3
    assert summary["samples"] == 2
    assert summary["missing_calls"] == 1
    assert summary["call_rate"] == 66.67
