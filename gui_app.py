#!/usr/bin/env python3
"""Streamlit GUI for PanSSR pipeline execution and output visualization."""

from __future__ import annotations

import shlex
import subprocess
import sys

import streamlit as st

from panssr.visualization_utils import (
    load_genotype_records,
    load_marker_records,
    summarize_genotypes,
    summarize_markers,
)


st.set_page_config(page_title="PanSSR GUI", layout="wide")


def _run_command(args: list[str]) -> tuple[int, str]:
    """Run a command and stream logs to UI, returning code and collected logs."""
    log_placeholder = st.empty()
    collected: list[str] = []

    process = subprocess.Popen(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None
    for line in process.stdout:
        collected.append(line)
        log_placeholder.code("".join(collected[-300:]), language="bash")

    return process.wait(), "".join(collected)


def _show_common_options() -> tuple[str, str]:
    output_path = st.text_input("Output file path", value="./output.tsv")
    python_executable = st.text_input("Python executable", value=sys.executable)
    return output_path, python_executable


def _render_genome_mode() -> None:
    st.subheader("Genome mode (SSR discovery + marker design)")
    genome_dir = st.text_input("Genome directory", value="./genomes")
    annot_dir = st.text_input("Annotation directory", value="./annotations")
    output_path, python_executable = _show_common_options()
    min_genome_support = st.number_input("Minimum genome support for common polymorphic set", min_value=1, value=1, step=1)

    if st.button("Run genome mode", type="primary"):
        cmd = [
            python_executable,
            "main.py",
            "--mode",
            "genome",
            "--genome_dir",
            genome_dir,
            "--annot_dir",
            annot_dir,
            "--output",
            output_path,
            "--min_genome_support",
            str(int(min_genome_support)),
        ]
        st.info(f"Running: {shlex.join(cmd)}")
        code, logs = _run_command(cmd)
        if code == 0:
            st.success(f"Genome mode finished successfully. Output: {output_path}")
        else:
            st.error(f"Genome mode failed with exit code {code}")
        st.download_button("Download run logs", logs, "panssr_genome_run.log", "text/plain")


def _render_genotype_mode() -> None:
    st.subheader("Genotype mode (BAM genotyping)")
    reference = st.text_input("Reference FASTA", value="./reference.fasta")
    markers = st.text_input("Markers TSV", value="./markers.tsv")
    bam_dir = st.text_input("BAM directory", value="./bams")
    output_path, python_executable = _show_common_options()
    workers = st.number_input("Genotype workers", min_value=1, value=1, step=1)

    if st.button("Run genotype mode", type="primary"):
        cmd = [
            python_executable,
            "main.py",
            "--mode",
            "genotype",
            "--reference",
            reference,
            "--markers",
            markers,
            "--bam_dir",
            bam_dir,
            "--output",
            output_path,
            "--workers",
            str(int(workers)),
        ]
        st.info(f"Running: {shlex.join(cmd)}")
        code, logs = _run_command(cmd)
        if code == 0:
            st.success(f"Genotype mode finished successfully. Output: {output_path}")
        else:
            st.error(f"Genotype mode failed with exit code {code}")
        st.download_button("Download run logs", logs, "panssr_genotype_run.log", "text/plain")


def _safe_bar_chart(data: dict, title: str) -> None:
    st.markdown(f"**{title}**")
    if not data:
        st.info("No data available")
        return
    try:
        st.bar_chart(data)
    except Exception:
        st.write(data)


def _render_visualization_studio() -> None:
    st.subheader("Visualization Studio")
    st.caption("Inspired by marker-centric dashboards in MegaSSR/Krait2 style workflows")

    marker_file = st.text_input("Marker TSV path", value="./markers.tsv")
    genotype_file = st.text_input("Genotype CSV path", value="./genotypes.csv")

    left, right = st.columns(2)

    with left:
        if st.button("Analyze marker file"):
            try:
                records = load_marker_records(marker_file)
                summary = summarize_markers(records)
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Total markers", summary["total_markers"])
                c2.metric("Unique motifs", summary["unique_motifs"])
                c3.metric("Chromosomes", summary["chromosomes"])
                c4.metric("Mean repeat count", summary["mean_repeat_count"])
                _safe_bar_chart(summary["motif_length_counts"], "Motif class distribution")
                _safe_bar_chart(summary["motif_counts"], "Top motifs")
                _safe_bar_chart(summary["chromosome_counts"], "Marker density by chromosome")
            except Exception as exc:
                st.error(f"Could not analyze marker file: {exc}")

    with right:
        if st.button("Analyze genotype file"):
            try:
                records = load_genotype_records(genotype_file)
                summary = summarize_genotypes(records)
                c1, c2, c3, c4 = st.columns(4)
                c1.metric("Total calls", summary["total_calls"])
                c2.metric("Samples", summary["samples"])
                c3.metric("Call rate (%)", summary["call_rate"])
                c4.metric("Missing calls", summary["missing_calls"])
                _safe_bar_chart(summary["sample_counts"], "Calls per sample")
                _safe_bar_chart(summary["genotype_counts"], "Top genotype states")
            except Exception as exc:
                st.error(f"Could not analyze genotype file: {exc}")


def main() -> None:
    st.title("PanSSR GUI")
    st.caption("Run PanSSR and visualize SSR/genotype output quality in one place")

    tab_run, tab_visual = st.tabs(["Run Pipelines", "Visualize Outputs"])

    with tab_run:
        mode = st.radio("Select workflow", ["Genome mode", "Genotype mode"], horizontal=True)
        if mode == "Genome mode":
            _render_genome_mode()
        else:
            _render_genotype_mode()

    with tab_visual:
        _render_visualization_studio()

    st.divider()
    st.markdown("**Tip:** run from repository root so `main.py` and relative paths resolve correctly.")
    st.markdown("`streamlit run gui_app.py`")


if __name__ == "__main__":
    main()
