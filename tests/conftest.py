"""Shared pytest fixtures for the corset-tools test suite."""

from pathlib import Path

import pytest


@pytest.fixture(scope='session')
def data_dir() -> Path:
    """
    Directory holding the checked-in test fixtures.

    Returns
    -------
    pathlib.Path
        Path to ``tests/data``.
    """
    return Path(__file__).parent / 'data'


@pytest.fixture
def cluster_map(data_dir: Path) -> Path:
    """
    Corset cluster map fixture.

    Parameters
    ----------
    data_dir : pathlib.Path
        Test data directory.

    Returns
    -------
    pathlib.Path
        Path to ``clusters.txt``.
    """
    return data_dir / 'clusters.txt'


@pytest.fixture
def transcripts(data_dir: Path) -> Path:
    """
    Multi-fasta of transcript sequences.

    Parameters
    ----------
    data_dir : pathlib.Path
        Test data directory.

    Returns
    -------
    pathlib.Path
        Path to ``transcripts.fa``.
    """
    return data_dir / 'transcripts.fa'


@pytest.fixture
def target_clusters(data_dir: Path) -> Path:
    """
    List of target cluster names, including one that does not exist.

    Parameters
    ----------
    data_dir : pathlib.Path
        Test data directory.

    Returns
    -------
    pathlib.Path
        Path to ``significantClusters.txt``.
    """
    return data_dir / 'significantClusters.txt'
