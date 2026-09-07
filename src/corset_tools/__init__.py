"""
Companion tools for working with transcript clusters produced by Corset.

The package exposes a single console entry point, ``corset-tools``, which
groups four subcommands:

``fetch-seqs``
    Extract and cluster-tag transcript sequences from a multi-FASTA.
``cross-count``
    Count cluster members originating from each of two transcriptomes.
``dist-calc``
    Estimate the alignment penalty required for reads to cross-map between
    two transcriptomes.
``edger-summary``
    Summarise and annotate differentially expressed clusters reported by edgeR.
"""

try:
    # Written at build time by hatch-vcs; absent when running from a source
    # checkout that has not been installed.
    from ._version import __version__
except ImportError:  # pragma: no cover - only hit in an uninstalled checkout
    __version__ = '0.0.0'

__all__ = ['__version__']
