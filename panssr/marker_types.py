"""Marker type annotation utilities (SSR/cSSR/iSSR/VNTR-like classes)."""

from __future__ import annotations


def classify_primary_type(marker: dict, vntr_min_repeat_units: int = 8) -> str:
    """Classify marker into SSR or VNTR-like category based on repeat count."""
    repeat_count = marker.get("repeat_count", 0)
    try:
        repeat_count = int(repeat_count)
    except (TypeError, ValueError):
        repeat_count = 0
    return "VNTR" if repeat_count >= vntr_min_repeat_units else "SSR"


def annotate_compound_and_interrupted(markers: list[dict], max_compound_gap: int = 10, max_interrupt_gap: int = 30) -> list[dict]:
    """
    Annotate markers with cSSR/iSSR style flags based on nearby loci on same chromosome.

    Rules:
    - cSSR: adjacent markers where gap <= max_compound_gap.
    - iSSR: same motif appears nearby with 0 < gap <= max_interrupt_gap.
    """
    grouped = {}
    for m in markers:
        grouped.setdefault(m.get("chrom"), []).append(m)

    for chrom_markers in grouped.values():
        chrom_markers.sort(key=lambda x: (int(x.get("start", 0)), int(x.get("end", 0))))
        for m in chrom_markers:
            m.setdefault("marker_type", classify_primary_type(m))
            m["is_cSSR"] = False
            m["is_iSSR"] = False

        for i in range(len(chrom_markers) - 1):
            a = chrom_markers[i]
            b = chrom_markers[i + 1]
            a_end = int(a.get("end", 0))
            b_start = int(b.get("start", 0))
            gap = b_start - a_end - 1

            if gap <= max_compound_gap:
                a["is_cSSR"] = True
                b["is_cSSR"] = True

            if 0 < gap <= max_interrupt_gap and str(a.get("motif", "")).upper() == str(b.get("motif", "")).upper():
                a["is_iSSR"] = True
                b["is_iSSR"] = True

    return markers
