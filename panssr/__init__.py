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

from importlib import import_module

__version__ = "1.0.0"
__author__ = "PanSSR Development Team"

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
    'genotype_utils',
    'visualization_utils',
    'marker_matrix',
    'marker_types',
]


def __getattr__(name):
    """Lazily load submodules to avoid hard failures on optional dependencies."""
    if name in __all__:
        module = import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
