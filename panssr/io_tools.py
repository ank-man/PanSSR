# panssr/io_tools.py
import os
from panssr import utils

def list_files_in_dir(directory, extensions=None):
    """
    List regular files in a directory; if extensions is provided, filter by extension.
    Returns files in sorted order for deterministic behavior.
    """
    files = []
    for entry in os.scandir(directory):
        if not entry.is_file():
            continue
        fname = entry.name
        if extensions and not any(fname.lower().endswith(ext.lower()) for ext in extensions):
            continue
        files.append(entry.path)
    return sorted(files)

def get_genome_annotation_pairs(genome_dir, annot_dir):
    """
    Returns a list of tuples (genome_fasta, annotation_file) by matching based on file basename.
    """
    genomes = list_files_in_dir(genome_dir, extensions=[".fa", ".fasta", ".fna"])
    annots = list_files_in_dir(annot_dir, extensions=[".gff", ".gtf"])

    annot_by_stem = {}
    for annot in annots:
        stem = os.path.splitext(os.path.basename(annot))[0]
        annot_by_stem.setdefault(stem, []).append(annot)

    pairs = []
    for genome in genomes:
        stem = os.path.splitext(os.path.basename(genome))[0]
        matching = annot_by_stem.get(stem, [])
        if len(matching) == 1:
            pairs.append((genome, matching[0]))
        elif len(matching) > 1:
            utils.logger.warning(
                "Multiple annotations found for genome %s (stem=%s); using first: %s",
                genome,
                stem,
                matching[0],
            )
            pairs.append((genome, matching[0]))
        else:
            utils.logger.warning("No annotation found for genome %s", genome)
    return pairs

def read_fasta(filepath):
    """
    Generator that yields (header, sequence) from a FASTA file.
    Uses a simple parser – for very large files, consider using pyfastx.
    """
    with open(filepath, "r") as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            yield header, "".join(seq_lines)
