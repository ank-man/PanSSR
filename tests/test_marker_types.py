"""Tests for marker type annotations (SSR/VNTR, cSSR, iSSR)."""

from panssr.marker_types import annotate_compound_and_interrupted, classify_primary_type


def test_classify_primary_type_vntr_threshold():
    assert classify_primary_type({"repeat_count": 3}) == "SSR"
    assert classify_primary_type({"repeat_count": 8}) == "VNTR"


def test_annotate_compound_and_interrupted():
    markers = [
        {"chrom": "chr1", "start": 10, "end": 20, "motif": "AT", "repeat_count": 5},
        {"chrom": "chr1", "start": 25, "end": 35, "motif": "AT", "repeat_count": 6},
        {"chrom": "chr1", "start": 100, "end": 110, "motif": "A", "repeat_count": 12},
    ]
    out = annotate_compound_and_interrupted(markers, max_compound_gap=10, max_interrupt_gap=30)
    assert out[0]["is_cSSR"] is True
    assert out[1]["is_cSSR"] is True
    assert out[0]["is_iSSR"] is True
    assert out[1]["is_iSSR"] is True
    assert out[2]["marker_type"] == "VNTR"
