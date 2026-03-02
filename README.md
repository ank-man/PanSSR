# PanSSR - Pangenomics SSR Discovery and Genotyping Tool

A comprehensive bioinformatics tool for discovering, annotating, and genotyping Simple Sequence Repeats (SSRs) across multiple genomic sequences at pangenome scale.

## Features

- **SSR Discovery**: Detect perfect tandem repeats (mono- to hexanucleotide motifs)
- **Marker Typing**: Annotate SSR/VNTR-like classes with cSSR/iSSR-style flags
- **Genomic Annotation**: Map SSRs to genomic features (exons, introns, genes, etc.)
- **Primer Design**: Automated PCR primer design using Primer3
- **In Silico PCR**: Validate primers with ePCR simulation (fuzzy matching support)
- **Genotyping**: CIGAR-aware genotype calling from BAM files
- **Quality Filtering**: Multi-criteria marker filtering for high-quality polymorphic markers

## Installation

### Using Conda (Recommended)

```bash
bash install_panssr.sh
```

### Manual Installation

```bash
# Create conda environment
conda create -n panssr python=3.9 -y
conda activate panssr

# Install dependencies
conda install -c conda-forge -c bioconda pysam primer3-py intervaltree numpy pandas tqdm pyfastx regex -y
```

## Usage


## GUI (Streamlit)

PanSSR includes a simple GUI for users who prefer not to run CLI commands directly.

```bash
streamlit run gui_app.py
```

The GUI supports:
- Genome mode inputs (`genome_dir`, `annot_dir`, `output`)
- Genotype mode inputs (`reference`, `markers`, `bam_dir`, `output`)
- Live run logs and downloadable log files
- A visualization studio inspired by MegaSSR/Krait2-style marker dashboards, including motif class plots, top motif charts, chromosome marker density, call-rate summaries, and genotype state distributions

> Run the GUI from the repository root so relative paths resolve against `main.py`.


### Mode 1: Genome Mode (SSR Discovery and Marker Design)

Discover SSRs, design primers, and create a polymorphic marker database.
PanSSR now also produces a **cross-genome marker matrix** and a **common polymorphic marker set** filtered by minimum genome support.


```bash
python main.py --mode genome \
  --genome_dir ./genomes/ \
  --annot_dir ./annotations/ \
  --output markers.tsv
```

**Input:**
- `--genome_dir`: Directory with genome FASTA files
- `--annot_dir`: Directory with annotation GFF/GTF files
- `--output`: Output TSV file for marker database
- `--min_genome_support`: Minimum number of genomes a marker must be present in to be included in `*.common_poly.tsv`

**Output:** `markers.tsv` containing:
- Chromosome, start, end positions
- SSR motif and repeat count
- Genomic annotation
- Primer sequences
- Expected amplicon sizes
- Marker type columns (`marker_type`, `is_cSSR`, `is_iSSR`)

Additional outputs:
- `<output>.matrix.tsv`: marker-by-genome matrix with repeat counts and polymorphism status
- `<output>.common_poly.tsv`: markers polymorphic and present in at least `--min_genome_support` genomes

### Mode 2: Genotype Mode (Population Genotyping)

Call genotypes at SSR loci from sequencing data:

```bash
python main.py --mode genotype \
  --reference ref.fasta \
  --markers markers.tsv \
  --bam_dir ./bams/ \
  --output genotypes.csv
```

**Input:**
- `--reference`: Reference genome FASTA file
- `--markers`: Marker database from genome mode
- `--bam_dir`: Directory with BAM files (indexed)
- `--output`: Output CSV file for genotype calls
- `--workers`: Number of worker processes for parallel BAM genotyping (HPC-friendly)

**Output:** `genotypes.csv` containing:
- Sample name
- Marker position
- Genotype calls (e.g., 10/12)
- Allele counts and read support

## Key Improvements in This Version

### 1. **Package Structure Fixed**
- Organized code into proper `panssr` package
- Added `__init__.py` for module imports
- Now properly importable and installable

### 2. **CIGAR-Aware Genotyping**
- Improved `genotyper.py` with alignment-aware SSR extraction
- Handles insertions/deletions correctly
- More accurate repeat counting from BAM files

### 3. **Enhanced ePCR Simulation**
- Fixed reverse complement handling for reverse primers
- Graceful fallback when `tre` module unavailable
- Supports both fuzzy and exact matching

### 4. **Robust Primer Design**
- Fixed coordinate calculation bugs
- Proper handling of edge cases (SSRs near sequence boundaries)
- Better error handling and validation

### 5. **Intelligent Marker Filtering**
- Validates that size differences match motif length multiples
- Filters non-specific amplification
- Configurable amplicon size ranges
- Optional intergenic marker preference

### 6. **Annotation Priority Ranking**
- Intelligent feature prioritization (CDS > exon > intron > gene)
- Handles overlapping features correctly
- Case-insensitive matching

### 7. **Comprehensive Error Handling**
- Input validation for all file paths
- Graceful error recovery
- Detailed logging throughout pipeline
- Clear error messages

### 8. **Test Suite**
- Unit tests for core functionality
- Test coverage for SSR discovery, filtering, and utilities
- Easy to run: `python -m pytest tests/`

## Configuration

Edit `panssr/config.py` to customize parameters:

### SSR Detection
```python
DEFAULT_MIN_REPEATS = {
    "mono": 12,
    "di": 7,
    "tri": 5,
    "tetra": 4,
    "penta": 4,
    "hexa": 4
}
MAX_SSR_LENGTH = 80  # bp
```

### Primer Design
```python
PRIMER_OPT_SIZE = 19  # bp
PRIMER_OPT_TM = 55    # °C
PRIMER_MIN_GC = 40    # %
PRIMER_MAX_GC = 70    # %
FLANK_SIZE = 100      # bp
```

### Genotyping
```python
MIN_MAPQ = 45              # Minimum mapping quality
MIN_READ_SUPPORT = 3       # Minimum reads for genotype call
```

## Dependencies

- **Python 3.9+**
- **Core libraries:**
  - `pysam` - BAM file handling
  - `primer3-py` - Primer design
  - `intervaltree` - Genomic interval operations
  - `numpy`, `pandas` - Data manipulation
  - `pyfastx` - Fast FASTA parsing
  - `regex` - Advanced pattern matching

- **Optional:**
  - `tre` - Fuzzy matching for ePCR (falls back to exact matching if unavailable)

## Architecture

```
panssr/
├── __init__.py           # Package initialization
├── config.py            # Configuration parameters
├── utils.py             # Utility functions
├── ssr_discovery.py     # SSR detection engine
├── annotator.py         # Genomic feature annotation
├── primer_design.py     # Primer3 interface
├── epcr.py             # In silico PCR simulation
├── genotyper.py         # BAM-based genotyping
├── marker_filter.py     # Marker quality filtering
├── database.py          # SQLite operations
├── report_generator.py  # HTML report generation
└── io_tools.py         # File I/O utilities

tests/
├── test_ssr_discovery.py
├── test_utils.py
├── test_marker_filter.py
└── test_annotator.py

main.py                  # Command-line interface
```

## Citation

If you use PanSSR in your research, please cite:

```
[Citation information to be added]
```

## License

[License information to be added]

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Support

For issues and questions:
- GitHub Issues: [https://github.com/ank-man/PanSSR/issues](https://github.com/ank-man/PanSSR/issues)

## Authors

PanSSR Development Team

## Acknowledgments

- Primer3 for primer design algorithms
- pysam for BAM file handling
- The bioinformatics community for valuable feedback
