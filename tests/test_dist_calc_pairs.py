"""Tests for transcript pair discovery."""

from pathlib import Path

import pytest

from corset_tools.dist_calc.pairs import read_blast_reciprocal_pairs, read_pairs
from corset_tools.exceptions import InputFormatError


def test_read_pairs_handles_ragged_input(data_dir: Path) -> None:
    """Comments, ragged tabs and extra columns are all tolerated.

    ``pairs.txt`` deliberately mixes tab and space separators, trailing empty
    fields, extra columns and commented-out pairs.
    """
    pairs = read_pairs(data_dir / 'pairs.txt')

    assert pairs == [
        ('name1', 'name2'),
        ('name3', 'name4'),
        ('name5', 'name6'),
        ('name7', 'name8'),
        ('nameX', 'nameY'),
    ]


def test_read_pairs_skips_single_field_lines(tmp_path: Path) -> None:
    """A line with only one name cannot form a pair and is skipped."""
    path = tmp_path / 'pairs.txt'
    path.write_text('a\tb\nlonely\nc\td\n')

    assert read_pairs(path) == [('a', 'b'), ('c', 'd')]


def test_read_pairs_rejects_empty(tmp_path: Path) -> None:
    """A file with no pairs is an error rather than an empty run."""
    path = tmp_path / 'pairs.txt'
    path.write_text('# nothing here\n')

    with pytest.raises(InputFormatError):
        read_pairs(path)


def _blast_line(query: str, subject: str, length: int, evalue: str) -> str:
    """
    Build one line of BLAST tabular (-outfmt 6) output.

    Parameters
    ----------
    query : str
        Query sequence name.
    subject : str
        Subject sequence name.
    length : int
        Alignment length, column 4.
    evalue : str
        E-value, column 11.

    Returns
    -------
    str
        A tab-delimited, newline-terminated BLAST record.
    """
    columns = [
        query,
        subject,
        '99.0',
        str(length),
        '0',
        '0',
        '1',
        '100',
        '1',
        '100',
        evalue,
        '200',
    ]
    return '\t'.join(columns) + '\n'


@pytest.fixture
def blast_tables(tmp_path: Path) -> tuple[Path, Path]:
    """
    Build a pair of BLAST tables with one reciprocal and one one-way best hit.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest temporary directory.

    Returns
    -------
    tuple of (pathlib.Path, pathlib.Path)
        The A-versus-B and B-versus-A table paths.
    """
    a_vs_b = tmp_path / 'AvB.tab'
    b_vs_a = tmp_path / 'BvA.tab'

    a_vs_b.write_text(
        _blast_line('a1', 'b1', 200, '1e-50')
        # Second hit for the same query is not the best hit and is ignored.
        + _blast_line('a1', 'b2', 150, '1e-20')
        + _blast_line('a2', 'b2', 200, '1e-40')
    )
    b_vs_a.write_text(
        _blast_line('b1', 'a1', 200, '1e-50')
        # b2's best hit is a3, so a2-b2 is not reciprocal.
        + _blast_line('b2', 'a3', 200, '1e-45')
    )
    return a_vs_b, b_vs_a


def test_reciprocal_best_hits(blast_tables: tuple[Path, Path]) -> None:
    """Only pairs that are each other's best hit in both directions survive."""
    a_vs_b, b_vs_a = blast_tables

    assert read_blast_reciprocal_pairs(a_vs_b, b_vs_a) == [('a1', 'b1')]


def test_reciprocal_best_hits_respects_filters(blast_tables: tuple[Path, Path]) -> None:
    """Hits failing the e-value or length filters are discarded."""
    a_vs_b, b_vs_a = blast_tables

    # Every hit is 200 long, so a floor of 200 excludes all of them.
    assert read_blast_reciprocal_pairs(a_vs_b, b_vs_a, min_len=200) == []
    # A stricter e-value than any hit achieves also excludes all of them.
    assert read_blast_reciprocal_pairs(a_vs_b, b_vs_a, e_value=1e-99) == []


def test_reciprocal_pairs_file(
    blast_tables: tuple[Path, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The optional pairs file is written next to the working directory."""
    a_vs_b, b_vs_a = blast_tables
    monkeypatch.chdir(tmp_path)

    read_blast_reciprocal_pairs(a_vs_b, b_vs_a, write_pairs_file=True)

    written = tmp_path / 'AvB_BvA_reciprocal_pairs.tab'
    assert written.read_text().splitlines() == ['#SetA\tSetB', 'a1\tb1']


def test_blast_rows_with_too_few_columns_are_skipped(tmp_path: Path) -> None:
    """Truncated BLAST rows are skipped rather than raising IndexError."""
    a_vs_b = tmp_path / 'AvB.tab'
    b_vs_a = tmp_path / 'BvA.tab'
    a_vs_b.write_text('a1\tb1\n' + _blast_line('a1', 'b1', 200, '1e-50'))
    b_vs_a.write_text(_blast_line('b1', 'a1', 200, '1e-50'))

    assert read_blast_reciprocal_pairs(a_vs_b, b_vs_a) == [('a1', 'b1')]
