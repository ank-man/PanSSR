"""
PanSSRAtor - A Pan-Species SSR Annotator for Pangenomics

This package provides tools for discovering, annotating, and genotyping
Simple Sequence Repeats (SSRs) across multiple genomic sequences at the
pangenome scale.

Main modules:
- ssr_discovery: SSR detection from DNA sequences
- annotator: Genomic feature annotation
- primer_design: PCR primer design using Primer3
- epcr: In silico PCR simulation
- genotyper: Genotype calling from BAM files
- marker_filter: Marker quality filtering
- database: SQLite database operations
- report_generator: HTML report generation
- io_tools: File I/O utilities
- utils: General utility functions
- config: Configuration parameters
"""

__version__ = "1.0.0"
__author__ = "PanSSR Development Team"

# Import main components for easier access
from . import config
from . import utils
from . import ssr_discovery
from . import annotator
from . import primer_design
from . import epcr
from . import genotyper
from . import marker_filter
from . import database
from . import report_generator
from . import io_tools

__all__ = [
    'config',
    'utils',
    'ssr_discovery',
    'annotator',
    'primer_design',
    'epcr',
    'genotyper',
    'marker_filter',
    'database',
    'report_generator',
    'io_tools',
]
