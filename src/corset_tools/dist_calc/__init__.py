"""
Transcriptome distance calculator.

Aligns matched transcript pairs drawn from two transcriptomes, converts the
observed gaps and mismatches into an expected per-read Bowtie2 alignment
penalty, and reports whether reads would cross-map under a given
``--score-min`` setting.
"""

from .align import AlignedPair, align_pair, count_gaps, count_mismatches, trim_alignment
from .pairs import read_blast_reciprocal_pairs, read_pairs
from .runner import DistCalcResult, run_dist_calc
from .scoring import (
    ScoreMinAdvice,
    load_ednafull,
    matrix_score,
    read_penalty,
    recommended_score_min,
    score_min,
)

__all__ = [
    'AlignedPair',
    'DistCalcResult',
    'ScoreMinAdvice',
    'align_pair',
    'count_gaps',
    'count_mismatches',
    'load_ednafull',
    'matrix_score',
    'read_blast_reciprocal_pairs',
    'read_pairs',
    'read_penalty',
    'recommended_score_min',
    'run_dist_calc',
    'score_min',
    'trim_alignment',
]
