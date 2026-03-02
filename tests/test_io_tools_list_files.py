"""Tests for file listing behavior."""

from pathlib import Path

from panssr import io_tools


def test_list_files_in_dir_filters_and_sorts(tmp_path):
    (tmp_path / "b.fa").write_text("x")
    (tmp_path / "a.fa").write_text("x")
    (tmp_path / "c.txt").write_text("x")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "nested.fa").write_text("x")

    files = io_tools.list_files_in_dir(str(tmp_path), extensions=[".fa"])
    names = [Path(f).name for f in files]

    assert names == ["a.fa", "b.fa"]
