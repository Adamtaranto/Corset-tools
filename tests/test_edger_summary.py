"""Tests for the edger-summary subcommand."""

from pathlib import Path

import pytest

from corset_tools.edger_summary import (
    add_membership_columns,
    filter_de_table,
    run_edger_summary,
    summarise_membership,
    write_annotation_records,
)
from corset_tools.fileio import read_cluster_map, read_fasta


@pytest.fixture
def toptags(tmp_path: Path) -> Path:
    """
    Write a small edgeR topTags export covering all three fixture clusters.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest temporary directory.

    Returns
    -------
    pathlib.Path
        Path to the written CSV.
    """
    path = tmp_path / 'toptags.csv'
    path.write_text(
        ',logFC,logCPM,LR,PValue,FDR\n'
        # Passes both thresholds.
        'Cluster-46274.0,3.5,8.1,40.2,1e-10,1e-8\n'
        # Passes, but with a weaker FDR so it sorts second.
        'Cluster-31069.0,-2.1,7.0,20.5,1e-5,0.001\n'
        # Fails both the fold-change and FDR thresholds.
        'Cluster-46275.0,0.2,6.5,1.1,0.3,0.4\n'
    )
    return path


@pytest.fixture
def annotations(tmp_path: Path) -> Path:
    """
    Write a free-form annotation file mentioning some cluster members.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Pytest temporary directory.

    Returns
    -------
    pathlib.Path
        Path to the written annotation file.
    """
    path = tmp_path / 'annot.txt'
    path.write_text(
        'q0Qmv2lr6i\tkinase domain\n'
        'RJ1s42IZgz\ttransporter\n'
        '6WoPXU2WEA\tmember of a non-significant cluster\n'
        'unrelated\tno transcript here\n'
    )
    return path


def test_filter_de_table_applies_thresholds(toptags: Path) -> None:
    """Rows failing the FDR or fold-change cut-offs are dropped."""
    table = filter_de_table(toptags, max_fdr=0.05, min_fc=1.0)

    assert list(table['ClusterID']) == ['Cluster-46274.0', 'Cluster-31069.0']
    # Sorted by ascending FDR.
    assert table['FDR'].is_monotonic_increasing


def test_filter_de_table_topx(toptags: Path) -> None:
    """--topx truncates after sorting by FDR."""
    table = filter_de_table(toptags, topx=1)

    assert list(table['ClusterID']) == ['Cluster-46274.0']


def test_filter_de_table_rejects_wrong_shape(tmp_path: Path) -> None:
    """A table without the six expected columns is rejected up front."""
    path = tmp_path / 'bad.csv'
    path.write_text('a,b\n1,2\n')

    with pytest.raises(ValueError, match='columns'):
        filter_de_table(path)


def test_summarise_membership_counts_by_label(
    cluster_map: Path, data_dir: Path
) -> None:
    """Every label is reported for every cluster, including zero counts."""
    clusters = read_cluster_map(cluster_map)
    labels = dict.fromkeys(read_fasta(data_dir / 'transcripts_X.fa'), 'SetX')
    labels.update(dict.fromkeys(read_fasta(data_dir / 'transcripts_Y.fa'), 'SetY'))

    summary = summarise_membership(
        clusters, labels, ['Cluster-46274.0', 'Cluster-31069.0'], ['SetX', 'SetY']
    )

    assert summary['Cluster-46274.0'] == {'SetX': 4, 'SetY': 5}
    assert summary['Cluster-31069.0'] == {'SetX': 2, 'SetY': 0}


def test_summarise_membership_handles_unknown_cluster(cluster_map: Path) -> None:
    """A DE cluster absent from the map reports zeros rather than raising."""
    clusters = read_cluster_map(cluster_map)

    summary = summarise_membership(clusters, {}, ['Cluster-ghost.0'], ['SetX'])

    assert summary == {'Cluster-ghost.0': {'SetX': 0}}


def test_add_membership_columns(toptags: Path) -> None:
    """Count columns and a Totals column are appended in label order."""
    table = filter_de_table(toptags)
    membership = {
        'Cluster-46274.0': {'SetX': 4, 'SetY': 5},
        'Cluster-31069.0': {'SetX': 2, 'SetY': 0},
    }

    merged = add_membership_columns(table, membership, ['SetX', 'SetY'])

    assert list(merged.columns[-3:]) == ['SetX', 'SetY', 'Totals']
    row = merged[merged.ClusterID == 'Cluster-46274.0'].iloc[0]
    assert row['SetX'] == 4
    assert row['Totals'] == 9


def test_write_annotation_records_buckets_by_cluster(
    tmp_path: Path, cluster_map: Path, annotations: Path
) -> None:
    """Matching lines are grouped under their cluster's header."""
    clusters = read_cluster_map(cluster_map)
    out_file = tmp_path / 'records.txt'

    write_annotation_records(
        clusters,
        annotations,
        out_file,
        ['Cluster-46274.0', 'Cluster-31069.0'],
    )

    assert out_file.read_text().splitlines() == [
        '#Cluster-46274.0',
        'q0Qmv2lr6i\tkinase domain',
        '#Cluster-31069.0',
        'RJ1s42IZgz\ttransporter',
    ]


def test_run_edger_summary_end_to_end(
    tmp_path: Path,
    toptags: Path,
    annotations: Path,
    cluster_map: Path,
    data_dir: Path,
) -> None:
    """The full pipeline writes a report, annotation records and cluster fastas."""
    result = run_edger_summary(
        infile=toptags,
        cluster_map=cluster_map,
        seq_files=[data_dir / 'transcripts_X.fa', data_dir / 'transcripts_Y.fa'],
        labels=['SetX', 'SetY'],
        grep_files=[annotations],
        out_dir=tmp_path / 'reports',
        write_clusters=True,
    )

    assert list(result.table['ClusterID']) == ['Cluster-46274.0', 'Cluster-31069.0']
    assert result.report_path.is_file()
    assert result.report_path.name == 'toptags_significant_cluster_report.tab'

    header = result.report_path.read_text().splitlines()[0]
    assert header.endswith('SetX\tSetY\tTotals')

    assert len(result.annotation_paths) == 1
    assert result.annotation_paths[0].name == (
        'annot_toptags_significant_cluster_annotation_records.txt'
    )

    assert len(result.fasta_paths) == 2
    # Both members of this cluster are present in transcripts_X.fa.
    written = read_fasta(tmp_path / 'reports' / 'Cluster-31069.0.fa')
    assert sorted(written) == ['RJ1s42IZgy', 'RJ1s42IZgz']
    # Record names carry transcript, source label and cluster on one line.
    first_header = (
        (tmp_path / 'reports' / 'Cluster-31069.0.fa').read_text().splitlines()[0]
    )
    assert first_header.split('\t')[1:] == ['SetX', 'Cluster-31069.0']


def test_run_edger_summary_skips_fastas_by_default(
    tmp_path: Path, toptags: Path, cluster_map: Path, data_dir: Path
) -> None:
    """Per-cluster fastas are written only when --writeClusters is set."""
    result = run_edger_summary(
        infile=toptags,
        cluster_map=cluster_map,
        seq_files=[data_dir / 'transcripts_X.fa'],
        labels=['SetX'],
        out_dir=tmp_path / 'reports',
    )

    assert result.fasta_paths == []
    assert not list((tmp_path / 'reports').glob('*.fa'))


def test_run_edger_summary_rejects_label_mismatch(
    tmp_path: Path, toptags: Path, cluster_map: Path, data_dir: Path
) -> None:
    """Labels and sequence files must be supplied in equal numbers."""
    with pytest.raises(ValueError, match='must match'):
        run_edger_summary(
            infile=toptags,
            cluster_map=cluster_map,
            seq_files=[data_dir / 'transcripts_X.fa'],
            labels=['SetX', 'SetY'],
            out_dir=tmp_path / 'reports',
        )
