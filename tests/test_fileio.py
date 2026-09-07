"""Tests for the shared IO helpers."""

import gzip
from pathlib import Path

import pytest

from corset_tools.exceptions import DuplicateTranscriptError, InputFormatError
from corset_tools.fileio import (
    invert_cluster_map,
    label_transcripts,
    package_data_path,
    read_cluster_map,
    read_fasta,
    read_multi_fasta,
    read_target_list,
    wrap_sequence,
    write_fasta,
)


def test_read_cluster_map_groups_members(cluster_map: Path) -> None:
    """Members are grouped under their cluster, preserving file order."""
    clusters = read_cluster_map(cluster_map)

    assert set(clusters) == {
        'Cluster-46274.0',
        'Cluster-31069.0',
        'Cluster-46275.0',
    }
    assert clusters['Cluster-46274.0'][0] == 'q0Qmv2lr6i'
    assert len(clusters['Cluster-46274.0']) == 10
    assert clusters['Cluster-31069.0'] == ['RJ1s42IZgz', 'RJ1s42IZgy']


def test_read_cluster_map_accepts_space_and_tab(tmp_path: Path) -> None:
    """Both tab- and space-delimited maps parse identically."""
    path = tmp_path / 'mixed.txt'
    path.write_text('t1\tcA\nt2 cA\n\n# a comment\nt3   cB\n')

    assert read_cluster_map(path) == {'cA': ['t1', 't2'], 'cB': ['t3']}


def test_read_cluster_map_skips_malformed_lines(tmp_path: Path) -> None:
    """Single-column lines are skipped rather than raising IndexError."""
    path = tmp_path / 'ragged.txt'
    path.write_text('t1\tcA\nlonely\nt2\tcA\n')

    assert read_cluster_map(path) == {'cA': ['t1', 't2']}


def test_read_cluster_map_rejects_empty_input(tmp_path: Path) -> None:
    """An input with no usable records is an error, not an empty result."""
    path = tmp_path / 'empty.txt'
    path.write_text('# only comments\n\n')

    with pytest.raises(InputFormatError):
        read_cluster_map(path)


def test_invert_cluster_map(cluster_map: Path) -> None:
    """Inverting the map yields a transcript to cluster lookup."""
    lookup = invert_cluster_map(read_cluster_map(cluster_map))

    assert lookup['q0Qmv2lr6i'] == 'Cluster-46274.0'
    assert lookup['RJ1s42IZgy'] == 'Cluster-31069.0'


def test_read_fasta(transcripts: Path) -> None:
    """Every record is loaded, keyed by its identifier."""
    sequences = read_fasta(transcripts)

    assert len(sequences) == 10
    assert sequences['q0Qmv2lr6i'].startswith('ATTGGCAGCTCTACAAATCGCC')


def test_read_fasta_handles_gzip(tmp_path: Path, transcripts: Path) -> None:
    """A gzipped fasta is read as text, not bytes."""
    gz_path = tmp_path / 'transcripts.fa.gz'
    with gzip.open(gz_path, 'wt', encoding='utf-8') as handle:
        handle.write(transcripts.read_text())

    assert read_fasta(gz_path) == read_fasta(transcripts)


def test_read_fasta_rejects_duplicate_ids(tmp_path: Path) -> None:
    """A repeated record name is fatal, since it would silently overwrite."""
    path = tmp_path / 'dupes.fa'
    path.write_text('>a\nACGT\n>a\nTGCA\n')

    with pytest.raises(DuplicateTranscriptError):
        read_fasta(path)


def test_label_transcripts(data_dir: Path) -> None:
    """Each transcript is labelled with the set it came from."""
    sets = read_multi_fasta(
        [data_dir / 'transcripts_X.fa', data_dir / 'transcripts_Y.fa'],
        ['SetX', 'SetY'],
    )
    labels = label_transcripts(sets)

    assert labels['q0Qmv2lr6i'] == 'SetX'
    assert labels['6WoPXU2WEA'] == 'SetY'


def test_label_transcripts_rejects_cross_set_duplicates(tmp_path: Path) -> None:
    """A name shared by two sets makes membership counts ambiguous."""
    a = tmp_path / 'a.fa'
    b = tmp_path / 'b.fa'
    a.write_text('>shared\nACGT\n')
    b.write_text('>shared\nTGCA\n')

    with pytest.raises(DuplicateTranscriptError):
        label_transcripts(read_multi_fasta([a, b], ['A', 'B']))


def test_read_multi_fasta_rejects_label_mismatch(tmp_path: Path) -> None:
    """Label and file counts must agree."""
    path = tmp_path / 'a.fa'
    path.write_text('>a\nACGT\n')

    with pytest.raises(ValueError, match='must match'):
        read_multi_fasta([path], ['A', 'B'])


def test_read_target_list_dedupes_in_order(tmp_path: Path) -> None:
    """Duplicates collapse while first-seen order is preserved."""
    path = tmp_path / 'targets.txt'
    path.write_text('cB\ncA\ncB\n\n#skip\ncC extra column\n')

    assert read_target_list(path) == ['cB', 'cA', 'cC']


def test_read_target_list_fixture(target_clusters: Path) -> None:
    """The checked-in target list includes a deliberately missing cluster."""
    assert read_target_list(target_clusters) == [
        'Cluster-46274.0',
        'Cluster-31069.0',
        'Cluster-fakeClust.0',
    ]


@pytest.mark.parametrize(
    ('sequence', 'width', 'expected'),
    [
        # An exact multiple of the width must not emit a trailing empty line.
        ('ABCDEF', 3, ['ABC', 'DEF']),
        # One over a multiple exercises the remainder branch that the original
        # float division broke.
        ('ABCDEFG', 3, ['ABC', 'DEF', 'G']),
        ('AB', 3, ['AB']),
        ('', 3, ['']),
    ],
)
def test_wrap_sequence(sequence: str, width: int, expected: list[str]) -> None:
    """Sequences wrap at the requested width without losing residues."""
    assert list(wrap_sequence(sequence, width)) == expected


def test_wrap_sequence_rejects_bad_width() -> None:
    """A non-positive width is rejected rather than looping forever."""
    with pytest.raises(ValueError, match='width'):
        list(wrap_sequence('ACGT', 0))


def test_write_fasta(tmp_path: Path) -> None:
    """A written record round-trips through the reader unchanged."""
    path = tmp_path / 'out.fa'
    sequence = 'ACGT' * 25

    with open(path, 'w', encoding='utf-8') as handle:
        write_fasta(handle, 'rec1', sequence, width=60)

    assert path.read_text().splitlines()[0] == '>rec1'
    assert read_fasta(path) == {'rec1': sequence}


def test_package_data_path_finds_matrix() -> None:
    """The EDNAFULL matrix ships with the package."""
    assert package_data_path('EDNAFULL.txt').is_file()


def test_package_data_path_missing() -> None:
    """A missing data file raises rather than returning a bad path."""
    with pytest.raises(FileNotFoundError):
        package_data_path('does_not_exist.txt')
