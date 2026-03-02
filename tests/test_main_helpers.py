"""Tests for helper functions in main module."""

from panssr.genotype_utils import filter_markers_by_reference, format_genotype_call


def test_filter_markers_by_reference_keeps_existing_chroms():
    markers = [
        {"chrom": "chr1", "start": 1, "end": 2},
        {"chrom": "chr2", "start": 3, "end": 4},
    ]
    ref = {"chr1": "AAAA"}
    filtered = filter_markers_by_reference(markers, ref)
    assert filtered == [{"chrom": "chr1", "start": 1, "end": 2}]


def test_format_genotype_call_handles_missing_and_sorts_alleles():
    gt, alleles = format_genotype_call({"genotype": None, "allele_counts": {}})
    assert gt == "NA"
    assert alleles == "NA"

    gt2, alleles2 = format_genotype_call(
        {
            "genotype": {"alleles": [12, 10]},
            "allele_counts": {12: 4, 10: 8},
        }
    )
    assert gt2 == "12/10"
    assert alleles2 == "10:8;12:4"
