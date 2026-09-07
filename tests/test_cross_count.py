"""Tests for the cross-count subcommand."""

from pathlib import Path

import pytest

from corset_tools.cross_count import (
    count_cluster_members,
    run_cross_count,
    summarise_counts,
)
from corset_tools.fileio import read_cluster_map


def test_run_cross_count_writes_table(
    tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """Counts are split by source transcriptome and written as a TSV."""
    out_file = tmp_path / 'counts.tab'
    counts, summary = run_cross_count(
        fasta_x=data_dir / 'transcripts_X.fa',
        fasta_y=data_dir / 'transcripts_Y.fa',
        cluster_map=cluster_map,
        name_x='SetX',
        name_y='SetY',
        out_file=out_file,
    )

    by_cluster = {row.cluster: row for row in counts}
    # Cluster-46274.0 has 10 listed members: 4 from X, 5 from Y and 'faketrans',
    # which is in neither transcriptome.
    assert by_cluster['Cluster-46274.0'].count_x == 4
    assert by_cluster['Cluster-46274.0'].count_y == 5
    assert by_cluster['Cluster-46274.0'].total == 10
    assert by_cluster['Cluster-46274.0'].unassigned == 1

    assert by_cluster['Cluster-31069.0'].count_x == 2
    assert by_cluster['Cluster-31069.0'].count_y == 0
    assert by_cluster['Cluster-46275.0'].count_x == 0
    assert by_cluster['Cluster-46275.0'].count_y == 2

    assert summary.n_clusters == 3
    assert summary.zero_x == 1
    assert summary.zero_y == 1
    assert summary.zero_both == 0

    lines = out_file.read_text().splitlines()
    assert lines[0] == 'clusterID\tSetX\tSetY\tTotal_Members'
    assert 'Cluster-31069.0\t2\t0\t2' in lines


def test_labels_default_to_fasta_basenames(
    tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """Omitting the labels uses each fasta's base name as the column heading."""
    out_file = tmp_path / 'counts.tab'
    _, summary = run_cross_count(
        fasta_x=data_dir / 'transcripts_X.fa',
        fasta_y=data_dir / 'transcripts_Y.fa',
        cluster_map=cluster_map,
        out_file=out_file,
    )

    # The resolved labels are reported back so the console echo matches the file.
    assert (summary.label_x, summary.label_y) == ('transcripts_X', 'transcripts_Y')

    header = out_file.read_text().splitlines()[0]
    assert header == 'clusterID\ttranscripts_X\ttranscripts_Y\tTotal_Members'


def test_identical_labels_rejected(
    tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """Two identical labels would make the count columns indistinguishable."""
    with pytest.raises(ValueError, match='must differ'):
        run_cross_count(
            fasta_x=data_dir / 'transcripts_X.fa',
            fasta_y=data_dir / 'transcripts_Y.fa',
            cluster_map=cluster_map,
            name_x='same',
            name_y='same',
            out_file=tmp_path / 'counts.tab',
        )


def test_summary_counts_zero_both(cluster_map: Path) -> None:
    """A cluster with no labelled members counts toward every zero total."""
    clusters = read_cluster_map(cluster_map)
    counts = count_cluster_members(clusters, {}, 'SetX', 'SetY')
    summary = summarise_counts(counts, 'SetX', 'SetY')

    assert summary.zero_x == summary.zero_y == summary.zero_both == 3
    assert (summary.label_x, summary.label_y) == ('SetX', 'SetY')
