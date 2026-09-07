"""
Build the list of matched transcript pairs to align.

Pairs come either from an explicit two-column table supplied by the user, or
are derived from reciprocal best BLAST hits between the two transcriptomes.
"""

from collections import Counter
import logging
import os
from typing import Optional

from ..exceptions import InputFormatError
from ..fileio import PathLike

logger = logging.getLogger(__name__)

# Column indices in NCBI BLAST tabular output (-outfmt 6).
BLAST_QUERY = 0
BLAST_SUBJECT = 1
BLAST_ALIGN_LEN = 3
BLAST_EVALUE = 10
BLAST_MIN_COLUMNS = 12


def read_pairs(path: PathLike) -> list[tuple[str, str]]:
    """
    Read an explicit table of transcript pairs.

    The first two whitespace-delimited fields of each line are taken as the
    pair. Trailing columns are ignored, and lines beginning with ``#`` are
    treated as comments.

    Parameters
    ----------
    path : str or os.PathLike
        Path to the pair table.

    Returns
    -------
    list of (str, str)
        Transcript pairs in file order.

    Raises
    ------
    InputFormatError
        If no usable pairs are found.
    """
    pairs: list[tuple[str, str]] = []
    with open(path, encoding='utf-8') as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            # The original code called map(str.strip, ...) and then subscripted
            # the result, which is a lazy iterator on Python 3.
            fields = [field.strip() for field in line.split()]
            if len(fields) < 2:
                logger.warning(
                    'Skipping line %d of %s, fewer than two fields: %r',
                    line_number,
                    path,
                    line,
                )
                continue
            pairs.append((fields[0], fields[1]))

    if not pairs:
        raise InputFormatError(f'No transcript pairs could be read from {path}')

    logger.info('Read %d transcript pairs from %s', len(pairs), path)
    return pairs


def _best_hits(
    path: PathLike, min_len: int, e_value: float, flip: bool = False
) -> list[tuple[str, str]]:
    """
    Extract the best hit per query from a BLAST tabular file.

    The best hit is the first surviving line for each query, which assumes the
    standard BLAST output ordering of descending hit quality within a query.

    Parameters
    ----------
    path : str or os.PathLike
        BLAST tabular (``-outfmt 6``) file.
    min_len : int
        Minimum alignment length to accept.
    e_value : float
        Exclusive e-value ceiling; hits at or above this are discarded.
    flip : bool, optional
        If True, emit ``(subject, query)`` rather than ``(query, subject)`` so
        that B-versus-A hits are expressed in A-versus-B orientation.

    Returns
    -------
    list of (str, str)
        Best-hit pairs.
    """
    pairs: list[tuple[str, str]] = []
    last_query: Optional[str] = None

    with open(path, encoding='utf-8') as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith('#'):
                continue
            fields = [field.strip() for field in line.split('\t')]
            if len(fields) < BLAST_MIN_COLUMNS:
                logger.warning(
                    'Skipping line %d of %s, expected %d BLAST columns, found %d',
                    line_number,
                    path,
                    BLAST_MIN_COLUMNS,
                    len(fields),
                )
                continue

            # The original compared the raw string column against an int, which
            # is a TypeError on Python 3; cast both numeric columns explicitly.
            try:
                hit_evalue = float(fields[BLAST_EVALUE])
                hit_length = int(fields[BLAST_ALIGN_LEN])
            except ValueError:
                logger.warning(
                    'Skipping line %d of %s, unparsable numeric columns',
                    line_number,
                    path,
                )
                continue

            if hit_evalue >= e_value or hit_length <= min_len:
                continue

            query, subject = fields[BLAST_QUERY], fields[BLAST_SUBJECT]
            if query != last_query:
                pairs.append((subject, query) if flip else (query, subject))
            last_query = query

    return pairs


def read_blast_reciprocal_pairs(
    blast_a_vs_b: PathLike,
    blast_b_vs_a: PathLike,
    min_len: int = 0,
    e_value: float = 0.001,
    write_pairs_file: bool = False,
) -> list[tuple[str, str]]:
    """
    Derive reciprocal best-hit pairs from two BLAST tabular files.

    Parameters
    ----------
    blast_a_vs_b : str or os.PathLike
        BLAST tabular output of transcriptome A queried against B.
    blast_b_vs_a : str or os.PathLike
        BLAST tabular output of transcriptome B queried against A.
    min_len : int, optional
        Minimum alignment length to accept, by default 0.
    e_value : float, optional
        Exclusive e-value ceiling, by default 0.001.
    write_pairs_file : bool, optional
        If True, write the reciprocal pairs to
        ``<A>_<B>_reciprocal_pairs.tab`` in the working directory.

    Returns
    -------
    list of (str, str)
        Pairs that were each other's best hit in both directions.
    """
    pairs = _best_hits(blast_a_vs_b, min_len, e_value, flip=False)
    pairs += _best_hits(blast_b_vs_a, min_len, e_value, flip=True)

    # A pair seen from both directions is a reciprocal best hit.
    reciprocal = [pair for pair, count in Counter(pairs).items() if count > 1]
    logger.info('Found %d reciprocal best-hit pairs', len(reciprocal))

    if write_pairs_file:
        base_a = os.path.splitext(os.path.basename(str(blast_a_vs_b)))[0]
        base_b = os.path.splitext(os.path.basename(str(blast_b_vs_a)))[0]
        out_name = f'{base_a}_{base_b}_reciprocal_pairs.tab'
        with open(out_name, 'w', encoding='utf-8') as handle:
            handle.write('#SetA\tSetB\n')
            for name_a, name_b in reciprocal:
                handle.write(f'{name_a}\t{name_b}\n')
        logger.info('Wrote reciprocal pairs to %s', out_name)

    return reciprocal
