"""
Global pairwise alignment of matched transcripts.

Uses :class:`Bio.Align.PairwiseAligner` in global (Needleman-Wunsch) mode.
The original implementation used ``Bio.pairwise2``, which was deprecated and
then removed from Biopython.
"""

from dataclasses import dataclass
import logging

from Bio import Align

logger = logging.getLogger(__name__)

# Alignment parameters retained from the original script: identical characters
# score 2, mismatches -1, gap opening -5 and extension -3. These shape the
# alignment only; the reported read penalties use the separate Bowtie2-style
# penalties supplied on the command line.
MATCH_SCORE = 2
MISMATCH_SCORE = -1
OPEN_GAP_SCORE = -5
EXTEND_GAP_SCORE = -3

GAP = '-'


@dataclass(frozen=True)
class AlignedPair:
    """
    A trimmed global alignment of two transcripts.

    Attributes
    ----------
    name_a : str
        Name of the transcript from set A.
    name_b : str
        Name of the transcript from set B.
    aligned_a : str
        Gapped alignment string for transcript A.
    aligned_b : str
        Gapped alignment string for transcript B.
    """

    name_a: str
    name_b: str
    aligned_a: str
    aligned_b: str

    @property
    def length(self) -> int:
        """
        Length of the trimmed alignment.

        Returns
        -------
        int
            Number of alignment columns.
        """
        return len(self.aligned_a)


def make_aligner() -> Align.PairwiseAligner:
    """
    Build the global aligner used for every transcript pair.

    Returns
    -------
    Bio.Align.PairwiseAligner
        Aligner configured with the package's global alignment penalties.
    """
    aligner = Align.PairwiseAligner(
        mode='global',
        match_score=MATCH_SCORE,
        mismatch_score=MISMATCH_SCORE,
        open_gap_score=OPEN_GAP_SCORE,
        extend_gap_score=EXTEND_GAP_SCORE,
    )
    return aligner


def trim_alignment(aligned_a: str, aligned_b: str) -> tuple[str, str]:
    """
    Remove terminal gap-only overhangs from an alignment.

    Leading and trailing columns are dropped until the first column in which
    both sequences have a residue, so that the reported alignment length and
    mismatch density are not diluted by the unaligned ends that global
    alignment necessarily produces. Columns that are a gap in both sequences
    are dropped anywhere they occur.

    Parameters
    ----------
    aligned_a : str
        Gapped alignment string for sequence A.
    aligned_b : str
        Gapped alignment string for sequence B.

    Returns
    -------
    tuple of (str, str)
        The trimmed alignment strings.
    """
    columns = [
        (a, b) for a, b in zip(aligned_a, aligned_b) if not (a == GAP and b == GAP)
    ]

    # Find the first and last columns where both sequences carry a residue.
    first = next((i for i, (a, b) in enumerate(columns) if a != GAP and b != GAP), None)
    if first is None:
        # No column aligns two residues; there is nothing meaningful to keep.
        return '', ''
    last = next(
        i
        for i in range(len(columns) - 1, -1, -1)
        if columns[i][0] != GAP and columns[i][1] != GAP
    )

    kept = columns[first : last + 1]
    return ''.join(a for a, _ in kept), ''.join(b for _, b in kept)


def count_gaps(aligned_a: str, aligned_b: str) -> int:
    """
    Count alignment columns containing a gap in either sequence.

    Parameters
    ----------
    aligned_a : str
        Gapped alignment string for sequence A.
    aligned_b : str
        Gapped alignment string for sequence B.

    Returns
    -------
    int
        Number of gapped columns.
    """
    return sum(a == GAP or b == GAP for a, b in zip(aligned_a, aligned_b))


def count_mismatches(aligned_a: str, aligned_b: str) -> int:
    """
    Count alignment columns where two differing residues are aligned.

    Gapped columns are not counted as mismatches; they are counted by
    :func:`count_gaps`.

    Parameters
    ----------
    aligned_a : str
        Gapped alignment string for sequence A.
    aligned_b : str
        Gapped alignment string for sequence B.

    Returns
    -------
    int
        Number of mismatched columns.
    """
    return sum(a != b and a != GAP and b != GAP for a, b in zip(aligned_a, aligned_b))


def align_pair(task: tuple[str, str, str, str]) -> AlignedPair:
    """
    Globally align one transcript pair and trim its terminal gaps.

    Takes a single packed tuple so that it can be dispatched directly through
    :meth:`multiprocessing.Pool.imap`.

    Parameters
    ----------
    task : tuple of (str, str, str, str)
        ``(sequence_a, sequence_b, name_a, name_b)``.

    Returns
    -------
    AlignedPair
        The trimmed alignment for the pair.
    """
    seq_a, seq_b, name_a, name_b = task
    logger.debug('Aligning %s-%s', name_a, name_b)

    # A fresh aligner per call keeps the function picklable and process-safe.
    alignment = make_aligner().align(seq_a, seq_b)[0]
    trimmed_a, trimmed_b = trim_alignment(alignment[0], alignment[1])

    return AlignedPair(
        name_a=name_a, name_b=name_b, aligned_a=trimmed_a, aligned_b=trimmed_b
    )


def format_alignment(pair: AlignedPair, score: float, bitscore: float) -> str:
    """
    Render an alignment as a human-readable three-line block.

    Parameters
    ----------
    pair : AlignedPair
        The alignment to render.
    score : float
        Substitution-matrix score for the alignment.
    bitscore : float
        Score divided by alignment length.

    Returns
    -------
    str
        A printable representation with a match/mismatch ruler between the
        two aligned sequences.
    """
    ruler = []
    for a, b in zip(pair.aligned_a, pair.aligned_b):
        if a == GAP or b == GAP:
            ruler.append(' ')
        elif a == b:
            ruler.append('|')
        else:
            ruler.append('.')

    return (
        f'{pair.aligned_a}\n'
        f'{"".join(ruler)}\n'
        f'{pair.aligned_b}\n'
        f'Score= {score:g}\n'
        f'BitScore= {bitscore:g}\n'
    )
