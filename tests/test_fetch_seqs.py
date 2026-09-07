"""Tests for the fetch-seqs subcommand."""

from pathlib import Path

from corset_tools.fetch_seqs import not_found_path, run_fetch_seqs
from corset_tools.fileio import read_cluster_map, read_fasta, read_target_list


def test_fetch_all_clusters(
    tmp_path: Path, transcripts: Path, cluster_map: Path
) -> None:
    """With no target list, every resolvable cluster member is written."""
    out_fasta = tmp_path / 'out.fa'
    result = run_fetch_seqs(
        in_fasta=transcripts,
        cluster_map=cluster_map,
        out_fasta=out_fasta,
        not_found_log=tmp_path / 'notfound.log',
    )

    sequences = read_fasta(out_fasta)
    assert result.written == len(sequences)
    # Records are renamed <cluster>_<transcript>.
    assert 'Cluster-46274.0_q0Qmv2lr6i' in sequences
    # 'faketrans' is in the map but not the fasta, so it is logged not written.
    assert ('Cluster-46274.0', 'faketrans') in result.missing_transcripts


def test_fetch_target_clusters_logs_missing(
    tmp_path: Path, transcripts: Path, cluster_map: Path, target_clusters: Path
) -> None:
    """Unresolvable clusters and transcripts are reported, not fatal."""
    out_fasta = tmp_path / 'out.fa'
    log = tmp_path / 'notfound.log'

    result = run_fetch_seqs(
        in_fasta=transcripts,
        cluster_map=cluster_map,
        targets=read_target_list(target_clusters),
        out_fasta=out_fasta,
        not_found_log=log,
    )

    assert result.missing_clusters == ['Cluster-fakeClust.0']
    assert ('Cluster-46274.0', 'faketrans') in result.missing_transcripts
    assert ('Cluster-31069.0', 'RJ1s42IZgy') in result.missing_transcripts

    log_text = log.read_text()
    assert 'Cluster-fakeClust.0' in log_text
    assert 'Cluster-46274.0\tfaketrans' in log_text


def test_fetch_longest_picks_one_per_cluster(
    tmp_path: Path, transcripts: Path, cluster_map: Path, target_clusters: Path
) -> None:
    """--longest emits exactly one, maximal-length, record per cluster."""
    out_fasta = tmp_path / 'long.fa'
    all_sequences = read_fasta(transcripts)
    clusters = read_cluster_map(cluster_map)
    targets = read_target_list(target_clusters)

    run_fetch_seqs(
        in_fasta=transcripts,
        cluster_map=cluster_map,
        targets=targets,
        out_fasta=out_fasta,
        longest=True,
        not_found_log=tmp_path / 'notfound.log',
    )

    written = read_fasta(out_fasta)
    # One record for each requested cluster that exists in the map.
    assert len(written) == 2

    for name, sequence in written.items():
        cluster, transcript = name.split('_', 1)
        assert sequence == all_sequences[transcript]
        # No other member of that cluster is longer than the one chosen.
        present = [
            all_sequences[member]
            for member in clusters[cluster]
            if member in all_sequences
        ]
        assert len(sequence) == max(len(candidate) for candidate in present)


def test_fetch_longest_survives_cluster_with_no_sequences(
    tmp_path: Path, transcripts: Path
) -> None:
    """A cluster whose members are all absent is skipped, not an error.

    The original implementation raised UnboundLocalError here, or silently
    re-emitted the previous cluster's sequence.
    """
    partial_map = tmp_path / 'map.txt'
    partial_map.write_text('q0Qmv2lr6i\tcGood\nghost1\tcEmpty\nghost2\tcEmpty\n')
    out_fasta = tmp_path / 'out.fa'

    result = run_fetch_seqs(
        in_fasta=transcripts,
        cluster_map=partial_map,
        out_fasta=out_fasta,
        longest=True,
        not_found_log=tmp_path / 'notfound.log',
    )

    assert result.empty_clusters == ['cEmpty']
    assert list(read_fasta(out_fasta)) == ['cGood_q0Qmv2lr6i']


def test_not_found_path_naming(tmp_path: Path) -> None:
    """The log name is derived from the target list, or a fixed default."""
    assert not_found_path(None) == 'NotFound_Clusters_Transcripts.log'
    assert not_found_path(tmp_path / 'myClusters.csv') == 'NotFound_myClusters.csv.log'
