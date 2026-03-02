"""Unit tests for io_tools module."""

from pathlib import Path

from panssr import io_tools


def touch(path: Path):
    path.write_text("dummy\n")


def test_get_genome_annotation_pairs_exact_stem_match(tmp_path):
    genome_dir = tmp_path / "genomes"
    annot_dir = tmp_path / "annotations"
    genome_dir.mkdir()
    annot_dir.mkdir()

    g1 = genome_dir / "sampleA.fa"
    g2 = genome_dir / "sampleB.fna"
    a1 = annot_dir / "sampleA.gff"
    a2 = annot_dir / "sampleB.gtf"
    decoy = annot_dir / "sampleA_extra.gff"

    for p in [g1, g2, a1, a2, decoy]:
        touch(p)

    pairs = io_tools.get_genome_annotation_pairs(str(genome_dir), str(annot_dir))
    pair_map = {Path(g).stem: Path(a).stem for g, a in pairs}

    assert pair_map["sampleA"] == "sampleA"
    assert pair_map["sampleB"] == "sampleB"
    assert "sampleA_extra" not in pair_map.values()


def test_get_genome_annotation_pairs_missing_annotation(tmp_path):
    genome_dir = tmp_path / "genomes"
    annot_dir = tmp_path / "annotations"
    genome_dir.mkdir()
    annot_dir.mkdir()

    g1 = genome_dir / "only.fa"
    touch(g1)

    pairs = io_tools.get_genome_annotation_pairs(str(genome_dir), str(annot_dir))
    assert pairs == []
