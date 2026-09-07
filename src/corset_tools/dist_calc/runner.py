"""Orchestration for the transcriptome distance calculation."""

from dataclasses import dataclass
import logging
from multiprocessing import Pool
from typing import Optional, Sequence

import numpy as np

from ..fileio import PathLike, read_fasta
from .align import (
    AlignedPair,
    align_pair,
    count_gaps,
    count_mismatches,
    format_alignment,
)
from .pairs import read_blast_reciprocal_pairs, read_pairs
from .plot import plot_penalty_distribution
from .scoring import (
    ScoreMinAdvice,
    matrix_score,
    read_penalty,
    recommended_score_min,
    score_min,
)

logger = logging.getLogger(__name__)

# Pairs whose lengths differ by more than this fraction are not comparable
# enough for a global alignment to be meaningful.
MIN_LENGTH_RATIO = 0.5

OUTPUT_HEADER = [
    'transNameA',
    'transNameB',
    'align_len',
    'gaps',
    'mismatches',
    'align_score',
    'mean_read_penalty',
    'cross_map',
]


@dataclass(frozen=True)
class PairStats:
    """
    Alignment statistics for a single transcript pair.

    Attributes
    ----------
    name_a : str
        Transcript name from set A.
    name_b : str
        Transcript name from set B.
    align_length : int
        Trimmed alignment length.
    gaps : int
        Gapped alignment columns.
    mismatches : int
        Mismatched alignment columns.
    score : float
        EDNAFULL substitution-matrix score.
    bitscore : float
        Matrix score divided by alignment length.
    read_penalty : float
        Expected per-read Bowtie2 penalty.
    cross_maps : bool
        Whether the pair is predicted to cross-map at the requested
        ``--score-min``.
    """

    name_a: str
    name_b: str
    align_length: int
    gaps: int
    mismatches: int
    score: float
    bitscore: float
    read_penalty: float
    cross_maps: bool


@dataclass
class DistCalcResult:
    """
    Outcome of a :func:`run_dist_calc` call.

    Attributes
    ----------
    stats : list of PairStats
        Per-pair alignment statistics.
    medians : dict of str to float
        Median alignment length, gaps, mismatches, score and read penalty.
    min_score : float
        Score minimum implied by the supplied intercept and slope.
    cross_map_pass : int
        Number of pairs predicted to cross-map.
    advice : ScoreMinAdvice
        Recommended score-minimum setting for the target percentile.
    """

    stats: list[PairStats]
    medians: dict[str, float]
    min_score: float
    cross_map_pass: int
    advice: ScoreMinAdvice


def select_comparable_pairs(
    pairs: Sequence[tuple[str, str]], sequences: dict[str, str]
) -> list[tuple[str, str, str, str]]:
    """
    Drop pairs that are missing a sequence or too different in length.

    Parameters
    ----------
    pairs : sequence of (str, str)
        Candidate transcript pairs.
    sequences : dict of str to str
        Sequences for both transcriptomes, keyed by transcript name.

    Returns
    -------
    list of (str, str, str, str)
        ``(sequence_a, sequence_b, name_a, name_b)`` tasks ready to align.
    """
    tasks: list[tuple[str, str, str, str]] = []
    for name_a, name_b in pairs:
        seq_a = sequences.get(name_a)
        seq_b = sequences.get(name_b)
        if seq_a is None or seq_b is None:
            missing = name_a if seq_a is None else name_b
            logger.warning(
                'Skipping pair %s-%s, no sequence for %s', name_a, name_b, missing
            )
            continue
        # Global alignment of wildly different lengths is dominated by the
        # terminal overhang rather than by genuine divergence.
        if len(seq_a) < MIN_LENGTH_RATIO * len(seq_b) or len(
            seq_b
        ) < MIN_LENGTH_RATIO * len(seq_a):
            logger.info(
                'Skipping pair %s-%s, lengths differ by more than %.0f%%',
                name_a,
                name_b,
                (1 - MIN_LENGTH_RATIO) * 100,
            )
            continue
        tasks.append((seq_a, seq_b, name_a, name_b))

    logger.info('Aligning %d of %d candidate pairs', len(tasks), len(pairs))
    return tasks


def align_all(
    tasks: Sequence[tuple[str, str, str, str]], processes: int = 4
) -> list[AlignedPair]:
    """
    Align every task, optionally across several processes.

    Parameters
    ----------
    tasks : sequence of (str, str, str, str)
        Alignment tasks from :func:`select_comparable_pairs`.
    processes : int, optional
        Number of worker processes, by default 4. A value of 1 runs inline,
        which keeps tracebacks readable and makes testing simpler.

    Returns
    -------
    list of AlignedPair
        Trimmed alignments, in task order.
    """
    if not tasks:
        return []
    if processes <= 1:
        return [align_pair(task) for task in tasks]

    # Alignment is CPU bound and each task is independent, so a plain process
    # pool is enough. `with` handles terminate/join, which the original code
    # got wrong by calling pool.join(10).
    with Pool(processes=processes) as pool:
        return list(pool.map(align_pair, tasks, chunksize=1))


def compute_stats(
    alignments: Sequence[AlignedPair],
    gap_open: float,
    gap_extend: float,
    mismatch: float,
    read_length: int,
    min_score: float,
    verbose: bool = False,
) -> list[PairStats]:
    """
    Derive per-pair statistics and cross-mapping predictions.

    Parameters
    ----------
    alignments : sequence of AlignedPair
        Trimmed alignments.
    gap_open : float
        Bowtie2 gap-open penalty.
    gap_extend : float
        Bowtie2 gap-extension penalty.
    mismatch : float
        Bowtie2 mismatch penalty.
    read_length : int
        Read length to model.
    min_score : float
        Score minimum a pair must reach to be predicted to cross-map.
    verbose : bool, optional
        If True, print each formatted alignment.

    Returns
    -------
    list of PairStats
        Statistics for every alignment of non-zero length.
    """
    stats: list[PairStats] = []
    for alignment in alignments:
        if alignment.length == 0:
            logger.warning(
                'Alignment %s-%s has no aligned columns, skipping',
                alignment.name_a,
                alignment.name_b,
            )
            continue

        score, bitscore = matrix_score(alignment.aligned_a, alignment.aligned_b)
        gaps = count_gaps(alignment.aligned_a, alignment.aligned_b)
        mismatches = count_mismatches(alignment.aligned_a, alignment.aligned_b)
        penalty = read_penalty(
            gaps,
            mismatches,
            alignment.length,
            gap_open,
            gap_extend,
            mismatch,
            read_length,
        )

        if verbose:
            print(f'\nAlignment: {alignment.name_a}-{alignment.name_b}')
            print(format_alignment(alignment, score, bitscore))

        stats.append(
            PairStats(
                name_a=alignment.name_a,
                name_b=alignment.name_b,
                align_length=alignment.length,
                gaps=gaps,
                mismatches=mismatches,
                score=score,
                bitscore=bitscore,
                read_penalty=penalty,
                cross_maps=min_score <= penalty,
            )
        )
    return stats


def write_stats_table(path: PathLike, stats: Sequence[PairStats]) -> None:
    """
    Write the per-pair statistics table as a TSV.

    Parameters
    ----------
    path : str or os.PathLike
        Output file path.
    stats : sequence of PairStats
        Statistics to write.

    Returns
    -------
    None
        The table is written to ``path``.
    """
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(OUTPUT_HEADER) + '\n')
        for row in stats:
            handle.write(
                '\t'.join(
                    [
                        row.name_a,
                        row.name_b,
                        str(row.align_length),
                        str(row.gaps),
                        str(row.mismatches),
                        f'{row.score:g}',
                        # `+ 0.0` normalises -0.0 to 0.0 in the output.
                        f'{row.read_penalty + 0.0:g}',
                        'SUCCESS' if row.cross_maps else 'FAIL',
                    ]
                )
                + '\n'
            )


def summarise(stats: Sequence[PairStats]) -> dict[str, float]:
    """
    Compute median statistics across all aligned pairs.

    Parameters
    ----------
    stats : sequence of PairStats
        Per-pair statistics.

    Returns
    -------
    dict of str to float
        Median ``len``, ``gaps``, ``mismatch``, ``score`` and ``readScore``.
    """
    if not stats:
        return {
            key: float('nan')
            for key in ('len', 'gaps', 'mismatch', 'score', 'readScore')
        }
    return {
        'len': float(np.median([s.align_length for s in stats])),
        'gaps': float(np.median([s.gaps for s in stats])),
        'mismatch': float(np.median([s.mismatches for s in stats])),
        'score': float(np.median([s.score for s in stats])),
        'readScore': float(np.median([s.read_penalty for s in stats])),
    }


def run_dist_calc(
    fasta_a: PathLike,
    fasta_b: PathLike,
    pair_names: Optional[PathLike] = None,
    blast_a_vs_b: Optional[PathLike] = None,
    blast_b_vs_a: Optional[PathLike] = None,
    read_length: int = 100,
    score_min_intercept: float = -0.6,
    score_min_slope: float = -0.6,
    gap_open: float = -5,
    gap_extend: float = -3,
    mismatch: float = -6,
    percentile: float = 95,
    out_file: PathLike = 'alignmentStats.txt',
    out_fig: Optional[PathLike] = 'readPenaltyDist.pdf',
    processes: int = 4,
    min_len: int = 0,
    e_value: float = 0.001,
    write_pairs_file: bool = False,
    verbose: bool = False,
) -> DistCalcResult:
    """
    Align matched transcript pairs and report cross-mapping expectations.

    Parameters
    ----------
    fasta_a : str or os.PathLike
        Transcriptome A FASTA.
    fasta_b : str or os.PathLike
        Transcriptome B FASTA.
    pair_names : str or os.PathLike, optional
        Explicit two-column table of transcript pairs. Takes precedence over
        the BLAST inputs.
    blast_a_vs_b : str or os.PathLike, optional
        BLAST tabular output of A queried against B.
    blast_b_vs_a : str or os.PathLike, optional
        BLAST tabular output of B queried against A.
    read_length : int, optional
        Read length to model, by default 100.
    score_min_intercept : float, optional
        Intercept of the Bowtie2 ``--score-min`` function, by default -0.6.
    score_min_slope : float, optional
        Slope of the Bowtie2 ``--score-min`` function, by default -0.6.
    gap_open : float, optional
        Bowtie2 gap-open penalty, by default -5.
    gap_extend : float, optional
        Bowtie2 gap-extension penalty, by default -3.
    mismatch : float, optional
        Bowtie2 mismatch penalty, by default -6.
    percentile : float, optional
        Target percentage of transcript pairs to allow to cross-map when
        recommending a score minimum, by default 95.
    out_file : str or os.PathLike, optional
        Output statistics table, by default ``alignmentStats.txt``.
    out_fig : str or os.PathLike, optional
        Penalty distribution figure. Pass None to skip plotting.
    processes : int, optional
        Worker processes to align with, by default 4.
    min_len : int, optional
        Minimum BLAST hit length, by default 0.
    e_value : float, optional
        BLAST e-value ceiling, by default 0.001.
    write_pairs_file : bool, optional
        If True, write the reciprocal best-hit pairs to file.
    verbose : bool, optional
        If True, print each formatted alignment.

    Returns
    -------
    DistCalcResult
        Per-pair statistics, medians, and the score-minimum recommendation.

    Raises
    ------
    ValueError
        If neither a pair table nor both BLAST tables are supplied, or if no
        pair could be aligned.
    """
    if pair_names:
        pairs = read_pairs(pair_names)
    elif blast_a_vs_b and blast_b_vs_a:
        pairs = read_blast_reciprocal_pairs(
            blast_a_vs_b,
            blast_b_vs_a,
            min_len=min_len,
            e_value=e_value,
            write_pairs_file=write_pairs_file,
        )
    else:
        raise ValueError(
            'Provide a list of transcript pairs (-n) or both BLAST tables '
            '(-x and -y) for reciprocal best-hit analysis.'
        )

    # Both transcriptomes share one lookup; names must therefore be unique
    # across the two files, which read_fasta enforces per file.
    sequences = read_fasta(fasta_a)
    sequences.update(read_fasta(fasta_b))

    tasks = select_comparable_pairs(pairs, sequences)
    alignments = align_all(tasks, processes=processes)

    min_score = score_min(score_min_intercept, score_min_slope, read_length)
    stats = compute_stats(
        alignments,
        gap_open=gap_open,
        gap_extend=gap_extend,
        mismatch=mismatch,
        read_length=read_length,
        min_score=min_score,
        verbose=verbose,
    )

    if not stats:
        raise ValueError(
            'No transcript pair could be aligned. Check that the pair names '
            'match the FASTA record names.'
        )

    write_stats_table(out_file, stats)
    medians = summarise(stats)
    penalties = [s.read_penalty for s in stats]

    advice = recommended_score_min(
        penalties,
        percentile=percentile,
        read_length=read_length,
        slope=score_min_slope,
        current_intercept=score_min_intercept,
    )

    if out_fig:
        plot_penalty_distribution(penalties, out_fig, threshold=advice.threshold)

    return DistCalcResult(
        stats=stats,
        medians=medians,
        min_score=min_score,
        cross_map_pass=sum(1 for s in stats if s.cross_maps),
        advice=advice,
    )
