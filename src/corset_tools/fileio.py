"""
Shared input/output helpers for corset-tools.

The original standalone scripts each carried their own copy of the Corset
cluster-map reader, the FASTA loader and the FASTA line-wrapper. Those
implementations disagreed in small ways (delimiters, wrap width, error
handling); this module is the single implementation used by every subcommand.
"""

from collections.abc import Iterable, Iterator, Sequence
import gzip
from importlib.resources import files
import logging
import os
from pathlib import Path
from typing import IO, Optional, Union

from Bio import SeqIO

from .exceptions import DuplicateTranscriptError, InputFormatError

logger = logging.getLogger(__name__)

# Type alias for the many "path-like" parameters in this module.
PathLike = Union[str, os.PathLike]


def package_data_path(name: str) -> Path:
    """
    Resolve the on-disk path of a file shipped inside ``corset_tools.data``.

    Parameters
    ----------
    name : str
        File name relative to the package data directory, e.g. ``EDNAFULL.txt``.

    Returns
    -------
    pathlib.Path
        Absolute path to the requested data file.

    Raises
    ------
    FileNotFoundError
        If the named file is not present in the package data directory.
    """
    resource = files('corset_tools.data').joinpath(name)
    path = Path(str(resource))
    if not path.is_file():
        raise FileNotFoundError(f'Packaged data file not found: {name}')
    return path


def _open_text(path: PathLike) -> IO[str]:
    """
    Open a plain or gzipped text file for reading.

    Parameters
    ----------
    path : str or os.PathLike
        File to open. A ``.gz`` suffix selects transparent decompression.

    Returns
    -------
    typing.IO[str]
        A file handle yielding ``str`` lines.
    """
    # The original scripts called gzip.open() without a mode, which yields
    # bytes on Python 3 and silently broke every downstream str comparison.
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, encoding='utf-8')


def read_cluster_map(path: PathLike) -> dict[str, list[str]]:
    """
    Read a Corset transcript-to-cluster map.

    The map is a two-column table of ``transcript`` and ``cluster``. Columns may
    be separated by tabs or spaces; blank lines and lines beginning with ``#``
    are ignored.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the Corset ``clusters.txt`` file.

    Returns
    -------
    dict of str to list of str
        Mapping of cluster identifier to the list of member transcript names,
        in the order encountered.

    Raises
    ------
    InputFormatError
        If no usable records are found in the file.
    """
    clusters: dict[str, list[str]] = {}
    malformed = 0

    with _open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Whitespace splitting keeps the deliberately delimiter-agnostic
            # behaviour of the original fetchClusterSeqs.py.
            fields = line.split()
            if len(fields) < 2:
                malformed += 1
                logger.warning(
                    'Skipping malformed cluster map line %d in %s: %r',
                    line_number,
                    path,
                    line,
                )
                continue
            transcript, cluster = fields[0], fields[1]
            clusters.setdefault(cluster, []).append(transcript)

    if not clusters:
        raise InputFormatError(
            f'No transcript/cluster records could be read from {path} '
            f'({malformed} malformed lines).'
        )

    logger.info(
        'Read %d clusters covering %d transcripts from %s',
        len(clusters),
        sum(len(members) for members in clusters.values()),
        path,
    )
    return clusters


def invert_cluster_map(clusters: dict[str, list[str]]) -> dict[str, str]:
    """
    Build a transcript-to-cluster lookup from a cluster-to-members map.

    Parameters
    ----------
    clusters : dict of str to list of str
        Mapping of cluster identifier to member transcript names.

    Returns
    -------
    dict of str to str
        Mapping of transcript name to its cluster identifier.
    """
    return {
        transcript: cluster
        for cluster, members in clusters.items()
        for transcript in members
    }


def read_fasta(path: PathLike) -> dict[str, str]:
    """
    Load a (optionally gzipped) FASTA file into a dictionary of sequences.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the FASTA file.

    Returns
    -------
    dict of str to str
        Mapping of record identifier to sequence string.

    Raises
    ------
    DuplicateTranscriptError
        If the same record identifier appears more than once.
    """
    sequences: dict[str, str] = {}
    with _open_text(path) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in sequences:
                raise DuplicateTranscriptError(
                    f'Duplicate sequence name {record.id!r} in {path}'
                )
            sequences[record.id] = str(record.seq)

    logger.info('Read %d sequences from %s', len(sequences), path)
    return sequences


def read_multi_fasta(
    paths: Sequence[PathLike], labels: Optional[Sequence[str]] = None
) -> dict[str, dict[str, str]]:
    """
    Load several FASTA files, keyed by a label for each file.

    Parameters
    ----------
    paths : sequence of str or os.PathLike
        FASTA files to load.
    labels : sequence of str, optional
        Label for each file. Defaults to each file's base name.

    Returns
    -------
    dict of str to (dict of str to str)
        Mapping of label to that file's sequence dictionary.

    Raises
    ------
    ValueError
        If ``labels`` is given and its length differs from ``paths``.
    """
    if labels is None:
        labels = [os.path.basename(str(path)) for path in paths]
    if len(labels) != len(paths):
        raise ValueError(
            f'Got {len(labels)} labels for {len(paths)} sequence files; '
            'these must match.'
        )
    return {label: read_fasta(path) for label, path in zip(labels, paths)}


def label_transcripts(sequence_sets: dict[str, dict[str, str]]) -> dict[str, str]:
    """
    Map every transcript name to the label of the set it came from.

    Parameters
    ----------
    sequence_sets : dict of str to (dict of str to str)
        Mapping of set label to that set's sequence dictionary, as returned by
        :func:`read_multi_fasta`.

    Returns
    -------
    dict of str to str
        Mapping of transcript name to its set label.

    Raises
    ------
    DuplicateTranscriptError
        If a transcript name occurs in more than one set. Cluster membership
        counts would be ambiguous in that case, so this is fatal.
    """
    labelled: dict[str, str] = {}
    for label, sequences in sequence_sets.items():
        for name in sequences:
            if name in labelled:
                raise DuplicateTranscriptError(
                    f'Transcript {name!r} occurs in both {labelled[name]!r} and '
                    f'{label!r}. Transcript names must be unique across input sets.'
                )
            labelled[name] = label
    return labelled


def read_target_list(path: PathLike) -> list[str]:
    """
    Read a list of target cluster names, one per line.

    The first whitespace-delimited field of each line is used, so additional
    annotation columns are tolerated. Blank lines and ``#`` comments are
    skipped and duplicates are removed while preserving first-seen order.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the target cluster list.

    Returns
    -------
    list of str
        Deduplicated cluster names in file order.
    """
    names: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            names.append(line.split()[0])

    # dict.fromkeys is an order-preserving dedupe, and replaces the original
    # O(n^2) `if not list.count(x)` check.
    unique = list(dict.fromkeys(names))
    logger.info('Read %d unique target clusters from %s', len(unique), path)
    return unique


def wrap_sequence(sequence: str, width: int = 60) -> Iterator[str]:
    """
    Split a sequence into fixed-width lines.

    Parameters
    ----------
    sequence : str
        Sequence to wrap.
    width : int, optional
        Line width in characters, by default 60.

    Yields
    ------
    str
        Successive lines of at most ``width`` characters. An empty sequence
        yields a single empty line so the FASTA record stays well formed.

    Raises
    ------
    ValueError
        If ``width`` is not a positive integer.
    """
    if width < 1:
        raise ValueError(f'Line width must be >= 1, got {width}')
    if not sequence:
        yield ''
        return
    # Slicing avoids the original `len(seq) / width` float that crashed
    # range() on Python 3, and needs no separate remainder handling.
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def write_fasta(handle: IO[str], name: str, sequence: str, width: int = 60) -> None:
    """
    Write one FASTA record to an open text handle.

    Parameters
    ----------
    handle : typing.IO[str]
        Open, writable text handle.
    name : str
        Record identifier, written after the ``>``.
    sequence : str
        Sequence to write.
    width : int, optional
        Line wrap width, by default 60.

    Returns
    -------
    None
        The record is written to ``handle``.
    """
    handle.write(f'>{name}\n')
    for line in wrap_sequence(sequence, width):
        handle.write(f'{line}\n')


def write_fasta_records(
    handle: IO[str], records: Iterable[tuple[str, str]], width: int = 60
) -> int:
    """
    Write several FASTA records to an open text handle.

    Parameters
    ----------
    handle : typing.IO[str]
        Open, writable text handle.
    records : iterable of (str, str)
        ``(name, sequence)`` pairs to write.
    width : int, optional
        Line wrap width, by default 60.

    Returns
    -------
    int
        Number of records written.
    """
    count = 0
    for name, sequence in records:
        write_fasta(handle, name, sequence, width)
        count += 1
    return count
