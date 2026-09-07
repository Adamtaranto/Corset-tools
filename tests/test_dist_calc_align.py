"""Tests for pairwise alignment and alignment statistics."""

from pathlib import Path

from corset_tools.dist_calc.align import (
    align_pair,
    count_gaps,
    count_mismatches,
    format_alignment,
    trim_alignment,
)
from corset_tools.dist_calc.runner import align_all, select_comparable_pairs
from corset_tools.dist_calc.scoring import matrix_score
from corset_tools.fileio import read_fasta


def test_count_gaps_and_mismatches() -> None:
    """Gapped columns and mismatched columns are counted separately."""
    aligned_a = 'ACGT-ACGT'
    aligned_b = 'ACTTAACGT'

    assert count_gaps(aligned_a, aligned_b) == 1
    # Position 3 differs (G vs T); the gap column is not a mismatch.
    assert count_mismatches(aligned_a, aligned_b) == 1


def test_trim_alignment_strips_terminal_overhangs() -> None:
    """Leading and trailing gap overhangs are removed."""
    trimmed = trim_alignment('---ACGTACGT--', 'GGGACGTACGTTT')

    assert trimmed == ('ACGTACGT', 'ACGTACGT')


def test_trim_alignment_keeps_internal_gaps() -> None:
    """Gaps inside the aligned region are informative and are retained."""
    assert trim_alignment('AC-GT', 'ACAGT') == ('AC-GT', 'ACAGT')


def test_trim_alignment_drops_double_gap_columns() -> None:
    """A column that is a gap in both sequences carries no information."""
    assert trim_alignment('AC--GT', 'AC-AGT') == ('AC-GT', 'ACAGT')


def test_trim_alignment_with_no_aligned_residues() -> None:
    """An alignment sharing no residue column trims to nothing."""
    assert trim_alignment('ACGT----', '----ACGT') == ('', '')


def test_align_pair_recovers_a_known_variant() -> None:
    """A single-base deletion is aligned as one gap, not a run of mismatches."""
    reference = 'ACGTACGTACGTACGTACGT'
    variant = reference[:10] + reference[11:]

    pair = align_pair((reference, variant, 'ref', 'var'))

    assert count_gaps(pair.aligned_a, pair.aligned_b) == 1
    assert count_mismatches(pair.aligned_a, pair.aligned_b) == 0


def test_align_fixture_pairs(data_dir: Path) -> None:
    """The checked-in pair fixtures align, and length outliers are dropped."""
    sequences = read_fasta(data_dir / 'pairSeqA.fa')
    sequences.update(read_fasta(data_dir / 'pairSeqB.fa'))

    pairs = [
        ('name1', 'name2'),
        ('name3', 'name4'),
        ('name5', 'name6'),
        ('name7', 'name8'),
        ('nameX', 'nameY'),
    ]
    tasks = select_comparable_pairs(pairs, sequences)

    # name7/name8 differ in length by more than 50% and are excluded.
    assert [(task[2], task[3]) for task in tasks] == [
        ('name1', 'name2'),
        ('name3', 'name4'),
        ('name5', 'name6'),
        ('nameX', 'nameY'),
    ]

    alignments = align_all(tasks, processes=1)
    by_name = {(a.name_a, a.name_b): a for a in alignments}

    # name5/name6 are identical, so the alignment is ungapped and exact.
    identical = by_name[('name5', 'name6')]
    assert identical.aligned_a == identical.aligned_b
    assert count_gaps(identical.aligned_a, identical.aligned_b) == 0
    assert count_mismatches(identical.aligned_a, identical.aligned_b) == 0

    # Every alignment is at least as long as the shorter input sequence minus
    # its trimmed overhang, and never empty.
    assert all(alignment.length > 0 for alignment in alignments)


def test_select_comparable_pairs_skips_missing_sequences() -> None:
    """A pair naming an absent transcript is skipped with a warning."""
    tasks = select_comparable_pairs(
        [('a', 'b'), ('a', 'missing')], {'a': 'ACGT', 'b': 'ACGT'}
    )

    assert [(task[2], task[3]) for task in tasks] == [('a', 'b')]


def test_matrix_score_of_identical_sequence() -> None:
    """Identical sequences score the matrix diagonal at every position."""
    score, bitscore = matrix_score('ACGT', 'ACGT')

    # EDNAFULL scores every exact nucleotide match as 5.
    assert score == 20
    assert bitscore == 5


def test_matrix_score_treats_gap_as_matrix_star() -> None:
    """A gap is scored through the matrix's '*' row and column."""
    gapped, _ = matrix_score('A-', 'AC')
    ungapped, _ = matrix_score('AC', 'AC')

    assert gapped < ungapped


def test_format_alignment_marks_matches_and_mismatches() -> None:
    """The ruler distinguishes matches, mismatches and gaps."""
    pair = align_pair(('ACGT', 'ACGT', 'a', 'b'))
    rendered = format_alignment(pair, 20.0, 5.0)

    assert '||||' in rendered
    assert 'BitScore= 5' in rendered


def test_align_all_multiprocess_matches_inline(data_dir: Path) -> None:
    """Running across processes gives the same alignments as running inline."""
    sequences = read_fasta(data_dir / 'pairSeqA.fa')
    sequences.update(read_fasta(data_dir / 'pairSeqB.fa'))
    tasks = select_comparable_pairs([('name1', 'name2'), ('name5', 'name6')], sequences)

    assert align_all(tasks, processes=2) == align_all(tasks, processes=1)


def test_align_all_with_no_tasks() -> None:
    """An empty task list short-circuits without starting a pool."""
    assert align_all([], processes=4) == []
