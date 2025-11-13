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
from panssrator import config, utils, io_tools, ssr_discovery, annotator, primer_design, epcr, genotyper, marker_filter

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

def genome_mode(genome_dir: str, annot_dir: str, output: str):
    """
    Genome mode: Discover SSRs, design primers, and create marker database.

    Parameters:
      genome_dir: Directory containing genome FASTA files
      annot_dir: Directory containing annotation GFF/GTF files
      output: Output file path for marker database
    """
    utils.logger.info("Running Genome Mode")

    # Validate input directories
    if not validate_input_directory(genome_dir, "Genome"):
        sys.exit(1)
    if not validate_input_directory(annot_dir, "Annotation"):
        sys.exit(1)

    # Get genome-annotation pairs
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

    all_markers = []

    for genome_file, annot_file in pairs:
        utils.logger.info("Processing genome: %s", genome_file)

        try:
            # Load annotation once per genome file
            annot_trees = annotator.load_annotation(annot_file)
        except Exception as e:
            utils.logger.error(f"Failed to load annotation from {annot_file}: {e}")
            continue

        # Process each sequence in the genome file
        try:
            for header, seq in io_tools.read_fasta(genome_file):
                utils.logger.info(f"  Processing sequence: {header} ({len(seq)} bp)")

                # Detect SSRs
                try:
                    ssrs = ssr_discovery.detect_ssrs(seq)
                    utils.logger.info(f"    Detected {len(ssrs)} SSRs")
                except Exception as e:
                    utils.logger.error(f"    SSR detection failed: {e}")
                    continue

                # Process each SSR
                for rec in ssrs:
                    rec["chrom"] = header

                    # Annotate SSR
                    try:
                        rec["annotation"] = annotator.annotate_ssr(rec, annot_trees)
                    except Exception as e:
                        utils.logger.warning(f"    Annotation failed for SSR at {rec['start']}: {e}")
                        rec["annotation"] = None

                    # Design primers
                    try:
                        rec["primers"] = primer_design.design_primers_for_ssr(rec, seq, flank=config.FLANK_SIZE)
                    except Exception as e:
                        utils.logger.warning(f"    Primer design failed for SSR at {rec['start']}: {e}")
                        rec["primers"] = {}

                    # Run ePCR simulation
                    if rec.get("primers"):
                        primer_pair = {
                            "forward": rec["primers"].get("PRIMER_LEFT_0_SEQUENCE", ""),
                            "reverse": rec["primers"].get("PRIMER_RIGHT_0_SEQUENCE", "")
                        }
                        try:
                            rec["amplicon_sizes"] = epcr.simulate_epcr(seq, primer_pair, max_cost=config.MAX_EPCR_COST)
                        except Exception as e:
                            utils.logger.warning(f"    ePCR simulation failed for SSR at {rec['start']}: {e}")
                            rec["amplicon_sizes"] = []
                    else:
                        rec["amplicon_sizes"] = []

                all_markers.extend(ssrs)

        except Exception as e:
            utils.logger.error(f"Failed to process genome file {genome_file}: {e}")
            continue

    utils.logger.info(f"Total markers detected across all genomes: {len(all_markers)}")

    # Filter markers
    try:
        filtered_markers = marker_filter.filter_markers(all_markers)
    except Exception as e:
        utils.logger.error(f"Marker filtering failed: {e}")
        sys.exit(1)

    if not filtered_markers:
        utils.logger.warning("No markers passed filtering criteria")

    # Write results to output
    try:
        with open(output, "w") as f:
            headers = ["chrom", "start", "end", "motif", "repeat_count", "annotation", "primers", "amplicon_sizes"]
            f.write("\t".join(headers) + "\n")
            for rec in filtered_markers:
                line = [str(rec.get(col, "")) for col in headers]
                f.write("\t".join(line) + "\n")
        utils.logger.info(f"Marker database saved to {output} ({len(filtered_markers)} markers)")
    except Exception as e:
        utils.logger.error(f"Failed to write output file {output}: {e}")
        sys.exit(1)

def genotype_mode(reference: str, markers_file: str, bam_dir: str, output: str):
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

    # Process each BAM file
    genotype_calls = {}
    for bam in bam_files:
        utils.logger.info(f"Processing BAM file: {bam}")
        genotype_calls[bam] = {}

        for marker in markers:
            # Check if reference sequence is available for this marker
            seq = ref_seqs.get(marker["chrom"], "")
            if not seq:
                utils.logger.warning(f"  Chromosome {marker['chrom']} not found in reference genome")
                continue

            # Genotype the marker
            try:
                gt = genotyper.genotype_marker(bam, marker)
                genotype_calls[bam][marker["start"]] = gt
            except Exception as e:
                utils.logger.error(f"  Failed to genotype marker at {marker['chrom']}:{marker['start']}: {e}")
                genotype_calls[bam][marker["start"]] = {"allele_counts": {}, "genotype": None}

    # Write genotype table as CSV
    try:
        with open(output, "w") as f:
            header_line = "BAM_file,Marker_chrom,Marker_start,Genotype,Allele_counts\n"
            f.write(header_line)
            for bam, calls in genotype_calls.items():
                bam_name = os.path.basename(bam)
                for marker_start, call in calls.items():
                    # Find corresponding marker for chromosome info
                    marker_chrom = "unknown"
                    for m in markers:
                        if m["start"] == marker_start:
                            marker_chrom = m["chrom"]
                            break

                    genotype = call.get("genotype", None)
                    allele_counts = call.get("allele_counts", {})

                    if genotype:
                        gt_str = "/".join(map(str, genotype.get("alleles", [])))
                    else:
                        gt_str = "NA"

                    allele_str = ";".join(f"{allele}:{count}" for allele, count in allele_counts.items())
                    if not allele_str:
                        allele_str = "NA"

                    f.write(f"{bam_name},{marker_chrom},{marker_start},{gt_str},{allele_str}\n")

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
    return parser.parse_args()

def main():
    args = parse_args()
    start = time.time()
    if args.mode == "genome":
        if not args.genome_dir or not args.annot_dir:
            utils.do_error("Genome mode requires --genome_dir and --annot_dir.")
        genome_mode(args.genome_dir, args.annot_dir, args.output)
    elif args.mode == "genotype":
        if not args.reference or not args.markers or not args.bam_dir:
            utils.do_error("Genotype mode requires --reference, --markers, and --bam_dir.")
        genotype_mode(args.reference, args.markers, args.bam_dir, args.output)
    end = time.time()
    utils.logger.info("PanSSRAtor run time: %.2f minutes", (end - start) / 60)

if __name__ == '__main__':
    main()

