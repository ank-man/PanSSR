#!/usr/bin/env python3
"""
PanSSRAtor – A Pan‑Species SSR Annotator

Usage:
  Genome Mode:
    python main.py --mode genome --genome_dir ./genomes/ --annot_dir ./annotations/ --output markers.tsv
  Genotype Mode:
    python main.py --mode genotype --reference ref.fasta --markers markers.tsv --bam_dir ./bams/ --output genotypes.csv
"""

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import pysam
from panssr import config, utils, io_tools, ssr_discovery, annotator, primer_design, epcr, genotyper, marker_filter
from panssr.genotype_utils import filter_markers_by_reference, format_genotype_call
from panssr.marker_matrix import build_marker_matrix, select_common_polymorphic
from panssr.marker_types import annotate_compound_and_interrupted, classify_primary_type

def validate_input_directory(dir_path: str, dir_type: str) -> bool:
    """
    Validate that a directory exists and is accessible.

    Parameters:
      dir_path: Path to directory
      dir_type: Description of directory type (for error messages)

    Returns:
      True if valid, False otherwise
    """
    if not os.path.exists(dir_path):
        utils.logger.error(f"{dir_type} directory does not exist: {dir_path}")
        return False
    if not os.path.isdir(dir_path):
        utils.logger.error(f"{dir_type} path is not a directory: {dir_path}")
        return False
    if not os.access(dir_path, os.R_OK):
        utils.logger.error(f"{dir_type} directory is not readable: {dir_path}")
        return False
    return True

def validate_input_file(file_path: str, file_type: str) -> bool:
    """
    Validate that a file exists and is accessible.

    Parameters:
      file_path: Path to file
      file_type: Description of file type (for error messages)

    Returns:
      True if valid, False otherwise
    """
    if not os.path.exists(file_path):
        utils.logger.error(f"{file_type} file does not exist: {file_path}")
        return False
    if not os.path.isfile(file_path):
        utils.logger.error(f"{file_type} path is not a file: {file_path}")
        return False
    if not os.access(file_path, os.R_OK):
        utils.logger.error(f"{file_type} file is not readable: {file_path}")
        return False
    return True



def genome_mode(genome_dir: str, annot_dir: str, output: str, min_genome_support: int = 1):
    """
    Genome mode: Discover SSRs, design primers, and create marker database.
    """
    utils.logger.info("Running Genome Mode")

    if not validate_input_directory(genome_dir, "Genome"):
        sys.exit(1)
    if not validate_input_directory(annot_dir, "Annotation"):
        sys.exit(1)

    try:
        pairs = io_tools.get_genome_annotation_pairs(genome_dir, annot_dir)
    except Exception as e:
        utils.logger.error(f"Failed to get genome-annotation pairs: {e}")
        sys.exit(1)

    if not pairs:
        utils.logger.error("No matching genome-annotation pairs found")
        utils.logger.error(f"  Genome directory: {genome_dir}")
        utils.logger.error(f"  Annotation directory: {annot_dir}")
        sys.exit(1)

    utils.logger.info(f"Found {len(pairs)} genome-annotation pair(s)")

    headers = [
        "genome_id", "chrom", "start", "end", "motif", "repeat_count", "annotation", "primers", "amplicon_sizes",
        "marker_type", "is_cSSR", "is_iSSR",
    ]

    filtered_markers_for_matrix = []
    total_detected = 0
    total_filtered = 0

    try:
        with open(output, "w") as f:
            f.write("\t".join(headers) + "\n")

            for genome_file, annot_file in pairs:
                utils.logger.info("Processing genome: %s", genome_file)
                genome_id = os.path.splitext(os.path.basename(genome_file))[0]
                genome_markers = []

                try:
                    annot_trees = annotator.load_annotation(annot_file)
                except Exception as e:
                    utils.logger.error(f"Failed to load annotation from {annot_file}: {e}")
                    continue

                try:
                    for header, seq in io_tools.read_fasta(genome_file):
                        utils.logger.info(f"  Processing sequence: {header} ({len(seq)} bp)")
                        try:
                            ssrs = ssr_discovery.detect_ssrs(seq)
                            utils.logger.info(f"    Detected {len(ssrs)} SSRs")
                        except Exception as e:
                            utils.logger.error(f"    SSR detection failed: {e}")
                            continue

                        for rec in ssrs:
                            rec["chrom"] = header
                            rec["genome_id"] = genome_id
                            rec["marker_type"] = classify_primary_type(rec)

                            try:
                                rec["annotation"] = annotator.annotate_ssr(rec, annot_trees)
                            except Exception as e:
                                utils.logger.warning(f"    Annotation failed for SSR at {rec['start']}: {e}")
                                rec["annotation"] = None

                            try:
                                rec["primers"] = primer_design.design_primers_for_ssr(rec, seq, flank=config.FLANK_SIZE)
                            except Exception as e:
                                utils.logger.warning(f"    Primer design failed for SSR at {rec['start']}: {e}")
                                rec["primers"] = {}

                            if rec.get("primers"):
                                primer_pair = {
                                    "forward": rec["primers"].get("PRIMER_LEFT_0_SEQUENCE", ""),
                                    "reverse": rec["primers"].get("PRIMER_RIGHT_0_SEQUENCE", ""),
                                }
                                try:
                                    rec["amplicon_sizes"] = epcr.simulate_epcr(seq, primer_pair, max_cost=config.MAX_EPCR_COST)
                                except Exception as e:
                                    utils.logger.warning(f"    ePCR simulation failed for SSR at {rec['start']}: {e}")
                                    rec["amplicon_sizes"] = []
                            else:
                                rec["amplicon_sizes"] = []

                        ssrs = annotate_compound_and_interrupted(ssrs)
                        genome_markers.extend(ssrs)

                except Exception as e:
                    utils.logger.error(f"Failed to process genome file {genome_file}: {e}")
                    continue

                total_detected += len(genome_markers)

                try:
                    filtered_markers = marker_filter.filter_markers(genome_markers)
                except Exception as e:
                    utils.logger.error(f"Marker filtering failed for genome {genome_id}: {e}")
                    continue

                total_filtered += len(filtered_markers)
                filtered_markers_for_matrix.extend(filtered_markers)

                for rec in filtered_markers:
                    row = [str(rec.get(col, "")) for col in headers]
                    f.write("\t".join(row) + "\n")

        utils.logger.info(f"Total markers detected across all genomes: {total_detected}")
        utils.logger.info(f"Total filtered markers written: {total_filtered}")
        utils.logger.info(f"Marker database saved to {output}")
    except Exception as e:
        utils.logger.error(f"Failed to write output file {output}: {e}")
        sys.exit(1)

    try:
        genome_ids, matrix_rows = build_marker_matrix(filtered_markers_for_matrix)
        matrix_file = f"{output}.matrix.tsv"
        with open(matrix_file, "w") as f:
            matrix_headers = [
                "chrom", "start", "end", "motif", "annotation", "genome_support", "polymorphic", "unique_repeat_counts",
            ] + genome_ids
            f.write("\t".join(matrix_headers) + "\n")
            for row in matrix_rows:
                base = [
                    str(row.get("chrom", "")),
                    str(row.get("start", "")),
                    str(row.get("end", "")),
                    str(row.get("motif", "")),
                    str(row.get("annotation", "")),
                    str(row.get("genome_support", 0)),
                    str(row.get("polymorphic", False)),
                    ",".join(map(str, row.get("unique_repeat_counts", []))),
                ]
                vals = [str(row["genome_values"].get(gid, "NA")) for gid in genome_ids]
                f.write("\t".join(base + vals) + "\n")

        common_rows = select_common_polymorphic(matrix_rows, max(1, min_genome_support))
        common_file = f"{output}.common_poly.tsv"
        with open(common_file, "w") as f:
            common_headers = [
                "chrom", "start", "end", "motif", "annotation", "genome_support", "unique_repeat_counts",
            ]
            f.write("\t".join(common_headers) + "\n")
            for row in common_rows:
                f.write(
                    "\t".join(
                        [
                            str(row.get("chrom", "")),
                            str(row.get("start", "")),
                            str(row.get("end", "")),
                            str(row.get("motif", "")),
                            str(row.get("annotation", "")),
                            str(row.get("genome_support", 0)),
                            ",".join(map(str, row.get("unique_repeat_counts", []))),
                        ]
                    )
                    + "\n"
                )
        utils.logger.info(
            f"Marker matrix saved to {matrix_file}; common polymorphic markers (>= {max(1, min_genome_support)} genomes) saved to {common_file}"
        )
    except Exception as e:
        utils.logger.error(f"Failed to build marker support matrix: {e}")


def _genotype_one_bam(bam_path: str, markers: list[dict]):
    """Worker helper for per-BAM genotyping."""
    calls = {}
    try:
        bam_handle = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        return bam_path, calls, f"Failed to open BAM file {bam_path}: {e}"

    try:
        for marker in markers:
            marker_key = (marker["chrom"], marker["start"], marker["end"])
            try:
                calls[marker_key] = genotyper.genotype_marker_with_handle(bam_handle, marker)
            except Exception:
                calls[marker_key] = {"allele_counts": {}, "genotype": None}
    finally:
        bam_handle.close()

    return bam_path, calls, None

def genotype_mode(reference: str, markers_file: str, bam_dir: str, output: str, workers: int = 1):
    """
    Genotype mode: Call genotypes at SSR loci from BAM files.

    Parameters:
      reference: Reference genome FASTA file
      markers_file: Marker database file (TSV from genome mode)
      bam_dir: Directory containing BAM files
      output: Output file path for genotype calls
    """
    utils.logger.info("Running Genotype Mode")

    # Validate input files and directories
    if not validate_input_file(reference, "Reference genome"):
        sys.exit(1)
    if not validate_input_file(markers_file, "Marker database"):
        sys.exit(1)
    if not validate_input_directory(bam_dir, "BAM"):
        sys.exit(1)

    # Load reference genome
    ref_seqs = {}
    try:
        utils.logger.info(f"Loading reference genome from {reference}")
        for header, seq in io_tools.read_fasta(reference):
            ref_seqs[header] = seq
            utils.logger.info(f"  Loaded sequence: {header} ({len(seq)} bp)")
    except Exception as e:
        utils.logger.error(f"Failed to load reference genome: {e}")
        sys.exit(1)

    if not ref_seqs:
        utils.logger.error("No sequences loaded from reference genome")
        sys.exit(1)

    # Load marker file (TSV)
    markers = []
    try:
        utils.logger.info(f"Loading markers from {markers_file}")
        with open(markers_file, "r") as f:
            header_line = next(f)  # skip header
            line_num = 1
            for line in f:
                line_num += 1
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    utils.logger.warning(f"  Skipping malformed line {line_num}: insufficient columns")
                    continue
                try:
                    marker = {
                        "chrom": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "motif": parts[3],
                        "repeat_count": int(parts[4])
                    }
                    markers.append(marker)
                except (ValueError, IndexError) as e:
                    utils.logger.warning(f"  Skipping malformed line {line_num}: {e}")
                    continue
        utils.logger.info(f"Loaded {len(markers)} markers")
    except Exception as e:
        utils.logger.error(f"Failed to load marker file: {e}")
        sys.exit(1)

    if not markers:
        utils.logger.error("No valid markers found in marker file")
        sys.exit(1)

    markers = filter_markers_by_reference(markers, ref_seqs)
    if not markers:
        utils.logger.error("No markers matched chromosomes in reference genome")
        sys.exit(1)

    # List BAM files
    try:
        bam_files = io_tools.list_files_in_dir(bam_dir, extensions=[".bam"])
    except Exception as e:
        utils.logger.error(f"Failed to list BAM files: {e}")
        sys.exit(1)

    if not bam_files:
        utils.logger.error(f"No BAM files found in {bam_dir}")
        sys.exit(1)

    utils.logger.info(f"Found {len(bam_files)} BAM file(s) to process")

    # Process each BAM file (optionally parallelized for HPC usage)
    genotype_calls = {}
    max_workers = max(1, int(workers))

    if max_workers == 1:
        for bam in bam_files:
            utils.logger.info(f"Processing BAM file: {bam}")
            bam_path, calls, error = _genotype_one_bam(bam, markers)
            if error:
                utils.logger.error(error)
            genotype_calls[bam_path] = calls
    else:
        utils.logger.info(f"Running genotyping in parallel with {max_workers} workers")
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_genotype_one_bam, bam, markers) for bam in bam_files]
            for fut in as_completed(futures):
                bam_path, calls, error = fut.result()
                if error:
                    utils.logger.error(error)
                genotype_calls[bam_path] = calls

    # Write genotype table as CSV
    try:
        with open(output, "w") as f:
            header_line = "BAM_file,Marker_chrom,Marker_start,Marker_end,Genotype,Allele_counts\n"
            f.write(header_line)
            for bam, calls in genotype_calls.items():
                bam_name = os.path.basename(bam)
                for (marker_chrom, marker_start, marker_end), call in sorted(calls.items()):
                    gt_str, allele_str = format_genotype_call(call)
                    f.write(f"{bam_name},{marker_chrom},{marker_start},{marker_end},{gt_str},{allele_str}\n")

        utils.logger.info(f"Genotype calls saved to {output}")
    except Exception as e:
        utils.logger.error(f"Failed to write output file {output}: {e}")
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(description="PanSSRAtor – Pan‑Species SSR Annotator")
    parser.add_argument("--mode", choices=["genome", "genotype"], required=True,
                        help="Select the mode of operation: genome (SSR discovery) or genotype (genotyping from BAM)")
    parser.add_argument("--genome_dir", help="Directory of genome FASTA files (for genome mode)")
    parser.add_argument("--annot_dir", help="Directory of annotation (GFF/GTF) files (for genome mode)")
    parser.add_argument("--reference", help="Reference genome FASTA (for genotype mode)")
    parser.add_argument("--markers", help="Marker file (from genome mode, for genotype mode)")
    parser.add_argument("--bam_dir", help="Directory of BAM files (for genotype mode)")
    parser.add_argument("--output", required=True, help="Output file (or prefix) for results")
    parser.add_argument("--min_genome_support", type=int, default=1,
                        help="For genome mode: minimum number of genomes for common polymorphic marker set")
    parser.add_argument("--workers", type=int, default=1,
                        help="Number of worker processes for genotype mode (parallel across BAM files)")
    return parser.parse_args()

def main():
    args = parse_args()
    start = time.time()
    if args.mode == "genome":
        if not args.genome_dir or not args.annot_dir:
            utils.do_error("Genome mode requires --genome_dir and --annot_dir.")
        genome_mode(args.genome_dir, args.annot_dir, args.output, args.min_genome_support)
    elif args.mode == "genotype":
        if not args.reference or not args.markers or not args.bam_dir:
            utils.do_error("Genotype mode requires --reference, --markers, and --bam_dir.")
        genotype_mode(args.reference, args.markers, args.bam_dir, args.output, args.workers)
    end = time.time()
    utils.logger.info("PanSSRAtor run time: %.2f minutes", (end - start) / 60)

if __name__ == '__main__':
    main()

