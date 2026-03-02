"""Tests for cross-genome marker matrix utilities."""

from panssr.marker_matrix import build_marker_matrix, select_common_polymorphic


def test_build_marker_matrix_and_common_selection():
    markers = [
        {"genome_id": "g1", "chrom": "chr1", "start": 10, "end": 20, "motif": "AT", "repeat_count": 5, "annotation": None},
        {"genome_id": "g2", "chrom": "chr1", "start": 10, "end": 20, "motif": "AT", "repeat_count": 7, "annotation": None},
        {"genome_id": "g3", "chrom": "chr1", "start": 10, "end": 20, "motif": "AT", "repeat_count": 7, "annotation": None},
        {"genome_id": "g1", "chrom": "chr2", "start": 30, "end": 40, "motif": "A", "repeat_count": 12, "annotation": None},
    ]

    genome_ids, rows = build_marker_matrix(markers)
    assert genome_ids == ["g1", "g2", "g3"]
    assert len(rows) == 2

    common = select_common_polymorphic(rows, min_genome_support=3)
    assert len(common) == 1
    row = common[0]
    assert row["chrom"] == "chr1"
    assert row["polymorphic"] is True
    assert row["genome_support"] == 3
    assert row["unique_repeat_counts"] == [5, 7]
