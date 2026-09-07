"""
Scoring, read-penalty estimation and Bowtie2 ``--score-min`` advice.

Bowtie2 rejects an alignment whose score falls below
``intercept + slope * readLength``. This module converts the divergence
observed between two aligned transcripts into the alignment penalty a read of
a given length would be expected to accrue, and works backwards from the
distribution of those penalties to the threshold a user must tolerate for a
nominated fraction of the transcript population to cross-map.
"""

from dataclasses import dataclass
from functools import lru_cache
import logging
from typing import Sequence

import numpy as np
import pandas as pd

from ..fileio import package_data_path

logger = logging.getLogger(__name__)

# The substitution matrix uses '*' for a gap; alignment strings use '-'.
GAP = '-'
MATRIX_GAP = '*'


@lru_cache(maxsize=1)
def load_ednafull() -> pd.DataFrame:
    """
    Load the packaged EDNAFULL nucleotide substitution matrix.

    Returns
    -------
    pandas.DataFrame
        Square matrix indexed and columned by IUPAC nucleotide codes.
    """
    path = package_data_path('EDNAFULL.txt')
    return pd.read_csv(path, header=0, index_col=0, sep=r'\s+')


def matrix_score(aligned_a: str, aligned_b: str) -> tuple[float, float]:
    """
    Score an alignment against the EDNAFULL substitution matrix.

    Parameters
    ----------
    aligned_a : str
        Gapped alignment string for sequence A.
    aligned_b : str
        Gapped alignment string for sequence B.

    Returns
    -------
    tuple of (float, float)
        The summed matrix score and the score divided by alignment length
        (zero for an empty alignment).
    """
    matrix = load_ednafull()
    score = 0.0
    for a, b in zip(aligned_a, aligned_b):
        row = MATRIX_GAP if a == GAP else a
        column = MATRIX_GAP if b == GAP else b
        # .at replaces the pandas .ix accessor removed in pandas 1.0.
        score += float(matrix.at[row, column])

    bitscore = score / len(aligned_a) if aligned_a else 0.0
    return score, bitscore


def read_penalty(
    gap_count: int,
    mismatch_count: int,
    align_length: int,
    gap_open: float,
    gap_extend: float,
    mismatch: float,
    read_length: int,
) -> float:
    """
    Estimate the mean Bowtie2 alignment penalty for a read of a given length.

    Gap and mismatch densities observed across the whole transcript alignment
    are scaled to the read length, then multiplied by the corresponding
    Bowtie2 penalties.

    Parameters
    ----------
    gap_count : int
        Gapped columns in the alignment.
    mismatch_count : int
        Mismatched columns in the alignment.
    align_length : int
        Total alignment columns.
    gap_open : float
        Bowtie2 gap-open penalty (negative).
    gap_extend : float
        Bowtie2 gap-extension penalty (negative).
    mismatch : float
        Bowtie2 mismatch penalty (negative).
    read_length : int
        Read length to model.

    Returns
    -------
    float
        Expected total penalty for one read, as a negative number.

    Raises
    ------
    ValueError
        If ``align_length`` is not positive.
    """
    if align_length <= 0:
        raise ValueError('Alignment length must be greater than zero.')

    mean_gaps_per_read = (gap_count / align_length) * read_length
    mean_mismatches_per_read = (mismatch_count / align_length) * read_length

    return (
        mismatch * mean_mismatches_per_read
        + (gap_open + gap_extend) * mean_gaps_per_read
    )


def score_min(intercept: float, slope: float, read_length: int) -> float:
    """
    Evaluate a Bowtie2 linear ``--score-min`` function.

    Parameters
    ----------
    intercept : float
        Intercept term of ``L,intercept,slope``.
    slope : float
        Slope term of ``L,intercept,slope``.
    read_length : int
        Read length.

    Returns
    -------
    float
        The minimum alignment score Bowtie2 will accept for that read length.
    """
    return intercept + slope * read_length


@dataclass(frozen=True)
class ScoreMinAdvice:
    """
    Recommended ``--score-min`` setting for a target cross-mapping rate.

    Attributes
    ----------
    target_percentile : float
        Percentage of transcript pairs the recommendation aims to cover.
    n_pairs : int
        Number of aligned pairs the recommendation is based on.
    read_length : int
        Read length the recommendation was computed for.
    slope : float
        Slope term held fixed while solving for the intercept.
    threshold : float
        Minimum alignment score that must be tolerated, i.e. the read penalty
        at the lower tail of the distribution.
    intercept : float
        Intercept that yields ``threshold`` at ``read_length``.
    achieved_fraction : float
        Fraction of pairs that actually meet ``threshold``.
    current_threshold : float, optional
        Score minimum implied by the user's current settings.
    current_fraction : float, optional
        Fraction of pairs that cross-map under the current settings.
    """

    target_percentile: float
    n_pairs: int
    read_length: int
    slope: float
    threshold: float
    intercept: float
    achieved_fraction: float
    current_threshold: float | None = None
    current_fraction: float | None = None

    @property
    def score_min_flag(self) -> str:
        """
        The recommendation as a ready-to-paste Bowtie2 flag.

        Returns
        -------
        str
            The flag string, e.g. ``--score-min L,-1.42,-0.6``.
        """
        return f'--score-min L,{self.intercept:.4g},{self.slope:g}'

    @property
    def current_is_sufficient(self) -> bool:
        """
        Whether the current settings already reach the target percentile.

        Returns
        -------
        bool
            True when the current threshold is at least as permissive as the
            recommendation. False when no current threshold was supplied.
        """
        if self.current_threshold is None:
            return False
        return self.current_threshold <= self.threshold


def recommended_score_min(
    read_penalties: Sequence[float],
    percentile: float,
    read_length: int,
    slope: float,
    current_intercept: float | None = None,
) -> ScoreMinAdvice:
    """
    Find the minimum alignment score needed for a target cross-mapping rate.

    A read pair cross-maps when the score minimum is no greater than the pair's
    expected read penalty. To let ``percentile`` percent of the transcript
    population cross-map, the score minimum must sit at or below the
    ``100 - percentile`` quantile of the penalty distribution.

    Parameters
    ----------
    read_penalties : sequence of float
        Expected per-read penalties, one per aligned transcript pair. These are
        negative numbers; less negative means more similar transcripts.
    percentile : float
        Target percentage of transcript pairs to allow to cross-map, e.g. 95.
    read_length : int
        Read length to express the recommendation for.
    slope : float
        Slope term of the Bowtie2 ``--score-min`` function to hold fixed while
        solving for the intercept.
    current_intercept : float, optional
        The user's current intercept. If given, the advice also reports what
        that setting achieves.

    Returns
    -------
    ScoreMinAdvice
        The recommended threshold and intercept, and how the current setting
        compares.

    Raises
    ------
    ValueError
        If ``read_penalties`` is empty or ``percentile`` is outside (0, 100].
    """
    if not len(read_penalties):
        raise ValueError('Cannot recommend a score minimum from zero alignments.')
    if not 0 < percentile <= 100:
        raise ValueError(f'Percentile must be in (0, 100], got {percentile}')

    penalties = np.asarray(read_penalties, dtype=float)

    # The lower tail of the penalty distribution is what must be tolerated:
    # covering 95% of pairs means accepting down to the 5th percentile.
    # method='lower' snaps to an observed penalty rather than interpolating
    # between two, so the recommendation always achieves at least the target
    # coverage instead of falling just short of it.
    threshold = float(np.percentile(penalties, 100.0 - percentile, method='lower'))
    achieved = float(np.mean(penalties >= threshold))

    # score_min(intercept, slope, L) == threshold  =>  intercept
    intercept = threshold - slope * read_length

    current_threshold = None
    current_fraction = None
    if current_intercept is not None:
        current_threshold = score_min(current_intercept, slope, read_length)
        current_fraction = float(np.mean(penalties >= current_threshold))

    advice = ScoreMinAdvice(
        target_percentile=percentile,
        n_pairs=int(penalties.size),
        read_length=read_length,
        slope=slope,
        threshold=threshold,
        intercept=intercept,
        achieved_fraction=achieved,
        current_threshold=current_threshold,
        current_fraction=current_fraction,
    )
    logger.info(
        'Recommended score-min for %.4g%% cross-mapping: %s',
        percentile,
        advice.score_min_flag,
    )
    return advice


def format_advice(advice: ScoreMinAdvice) -> str:
    """
    Render a :class:`ScoreMinAdvice` as a short human-readable report.

    Parameters
    ----------
    advice : ScoreMinAdvice
        The recommendation to render.

    Returns
    -------
    str
        A multi-line block suitable for printing to the terminal.
    """
    lines = [
        f'To allow cross-mapping for {advice.target_percentile:g}% of transcript '
        f'pairs (n={advice.n_pairs}) at readLength={advice.read_length}:',
        f'  minimum tolerable alignment score : {advice.threshold:.4g}',
        f'  suggested bowtie2 setting         : {advice.score_min_flag}',
        f'  covers                            : '
        f'{advice.achieved_fraction * 100:.1f}% of pairs',
    ]
    if advice.current_threshold is not None:
        verdict = 'sufficient' if advice.current_is_sufficient else 'TOO STRICT'
        lines.append(
            f'  your current setting              : '
            f'{advice.current_threshold:.4g} -> covers '
            f'{advice.current_fraction * 100:.1f}% of pairs ({verdict})'
        )
    return '\n'.join(lines)
