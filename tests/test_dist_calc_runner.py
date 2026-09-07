"""End-to-end tests for the transcriptome distance calculation."""

from pathlib import Path

import pytest

from corset_tools.dist_calc.runner import run_dist_calc


def test_run_dist_calc_end_to_end(tmp_path: Path, data_dir: Path) -> None:
    """Aligning the fixture pairs produces a stats table and a recommendation."""
    out_file = tmp_path / 'stats.txt'
    out_fig = tmp_path / 'dist.pdf'

    result = run_dist_calc(
        fasta_a=data_dir / 'pairSeqA.fa',
        fasta_b=data_dir / 'pairSeqB.fa',
        pair_names=data_dir / 'pairs.txt',
        out_file=out_file,
        out_fig=out_fig,
        processes=1,
    )

    # name7/name8 are excluded on the length-ratio filter, leaving four pairs.
    assert len(result.stats) == 4
    assert result.min_score == pytest.approx(-60.6)
    # The fixture pairs are all close variants, so all cross-map at the default.
    assert result.cross_map_pass == 4
    assert result.advice.n_pairs == 4
    assert result.advice.achieved_fraction >= 0.95

    lines = out_file.read_text().splitlines()
    assert lines[0].split('\t') == [
        'transNameA',
        'transNameB',
        'align_len',
        'gaps',
        'mismatches',
        'align_score',
        'mean_read_penalty',
        'cross_map',
    ]
    assert len(lines) == 5
    assert all('SUCCESS' in line for line in lines[1:])

    # name5/name6 are identical, so their penalty is exactly zero.
    identical = next(s for s in result.stats if s.name_a == 'name5')
    assert identical.gaps == 0
    assert identical.mismatches == 0
    assert identical.read_penalty == 0

    assert out_fig.is_file()
    assert out_fig.stat().st_size > 0


def test_run_dist_calc_can_skip_the_figure(tmp_path: Path, data_dir: Path) -> None:
    """Passing out_fig=None skips plotting entirely."""
    result = run_dist_calc(
        fasta_a=data_dir / 'pairSeqA.fa',
        fasta_b=data_dir / 'pairSeqB.fa',
        pair_names=data_dir / 'pairs.txt',
        out_file=tmp_path / 'stats.txt',
        out_fig=None,
        processes=1,
    )

    assert not list(tmp_path.glob('*.pdf'))
    assert result.stats


def test_strict_score_min_makes_pairs_fail(tmp_path: Path, data_dir: Path) -> None:
    """A stricter threshold is reflected in both the table and the advice."""
    out_file = tmp_path / 'stats.txt'

    result = run_dist_calc(
        fasta_a=data_dir / 'pairSeqA.fa',
        fasta_b=data_dir / 'pairSeqB.fa',
        pair_names=data_dir / 'pairs.txt',
        score_min_intercept=-0.6,
        score_min_slope=-0.1,
        out_file=out_file,
        out_fig=None,
        processes=1,
    )

    # -0.6 + (-0.1 * 100) = -10.6, which most of these pairs cannot reach.
    assert result.min_score == pytest.approx(-10.6)
    assert result.cross_map_pass < len(result.stats)
    assert 'FAIL' in out_file.read_text()
    # The advice must therefore be more permissive than the current setting.
    assert not result.advice.current_is_sufficient


def test_run_dist_calc_requires_pairs_or_blast(tmp_path: Path, data_dir: Path) -> None:
    """Neither a pair table nor BLAST tables means there is nothing to do."""
    with pytest.raises(ValueError, match='Provide a list of transcript pairs'):
        run_dist_calc(
            fasta_a=data_dir / 'pairSeqA.fa',
            fasta_b=data_dir / 'pairSeqB.fa',
            out_file=tmp_path / 'stats.txt',
            out_fig=None,
        )


def test_run_dist_calc_with_no_alignable_pairs(tmp_path: Path, data_dir: Path) -> None:
    """Pair names that match no sequence produce a clear error."""
    pairs = tmp_path / 'pairs.txt'
    pairs.write_text('ghostA\tghostB\n')

    with pytest.raises(ValueError, match='No transcript pair could be aligned'):
        run_dist_calc(
            fasta_a=data_dir / 'pairSeqA.fa',
            fasta_b=data_dir / 'pairSeqB.fa',
            pair_names=pairs,
            out_file=tmp_path / 'stats.txt',
            out_fig=None,
        )
