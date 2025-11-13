# panssrator/io_tools.py
import os
from panssrator import utils

def list_files_in_dir(directory, extensions=None):
    """
    List files in a directory; if extensions is provided, filter files by extension.
    """
    files = []
    for fname in os.listdir(directory):
        if extensions:
            if any(fname.lower().endswith(ext.lower()) for ext in extensions):
                files.append(os.path.join(directory, fname))
        else:
            files.append(os.path.join(directory, fname))
    return files

def get_genome_annotation_pairs(genome_dir, annot_dir):
    """
    Returns a list of tuples (genome_fasta, annotation_file) by matching based on file basename.
    """
    genomes = list_files_in_dir(genome_dir, extensions=[".fa", ".fasta", ".fna"])
    annots = list_files_in_dir(annot_dir, extensions=[".gff", ".gtf"])
    pairs = []
    for genome in genomes:
        base = os.path.splitext(os.path.basename(genome))[0]
        matching = [a for a in annots if base in os.path.basename(a)]
        if matching:
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

