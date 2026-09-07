"""Tests for the click command line interface."""

from pathlib import Path

from click.testing import CliRunner
import pytest

from corset_tools import __version__
from corset_tools.cli import main
from corset_tools.fileio import read_fasta


@pytest.fixture
def runner() -> CliRunner:
    """
    Build a click test runner.

    Returns
    -------
    click.testing.CliRunner
        Runner used to invoke the CLI in-process.
    """
    return CliRunner()


def test_version(runner: CliRunner) -> None:
    """--version reports the installed package version."""
    result = runner.invoke(main, ['--version'])

    assert result.exit_code == 0
    assert __version__ in result.output


def test_group_help_lists_every_subcommand(runner: CliRunner) -> None:
    """All four tools are reachable from the single entry point."""
    result = runner.invoke(main, ['--help'])

    assert result.exit_code == 0
    for command in ('fetch-seqs', 'cross-count', 'dist-calc', 'edger-summary'):
        assert command in result.output


@pytest.mark.parametrize(
    'command', ['fetch-seqs', 'cross-count', 'dist-calc', 'edger-summary']
)
def test_subcommand_help(runner: CliRunner, command: str) -> None:
    """Each subcommand renders its own help without executing anything."""
    result = runner.invoke(main, [command, '--help'])

    assert result.exit_code == 0
    assert 'Usage:' in result.output


def test_fetch_seqs_command(
    runner: CliRunner,
    tmp_path: Path,
    transcripts: Path,
    cluster_map: Path,
    target_clusters: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """fetch-seqs runs end to end and writes its outputs."""
    monkeypatch.chdir(tmp_path)
    out_fasta = tmp_path / 'out.fa'

    result = runner.invoke(
        main,
        [
            '--loglevel',
            'ERROR',
            'fetch-seqs',
            '-i',
            str(transcripts),
            '-c',
            str(cluster_map),
            '-t',
            str(target_clusters),
            '-o',
            str(out_fasta),
        ],
    )

    assert result.exit_code == 0, result.output
    assert 'Cluster-46274.0_q0Qmv2lr6i' in read_fasta(out_fasta)
    # The not-found log is named after the target list.
    assert (tmp_path / 'NotFound_significantClusters.txt.log').is_file()


def test_cross_count_command(
    runner: CliRunner, tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """cross-count prints its summary and writes the count table."""
    out_file = tmp_path / 'counts.tab'

    result = runner.invoke(
        main,
        [
            '--loglevel',
            'ERROR',
            'cross-count',
            '-X',
            str(data_dir / 'transcripts_X.fa'),
            '-Y',
            str(data_dir / 'transcripts_Y.fa'),
            '-c',
            str(cluster_map),
            '-x',
            'SetX',
            '-y',
            'SetY',
            '-o',
            str(out_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert 'Clusters with 0 members from SetX: 1' in result.output
    assert out_file.is_file()


def test_dist_calc_command_reports_advice(
    runner: CliRunner, tmp_path: Path, data_dir: Path
) -> None:
    """dist-calc prints the minimum tolerable alignment score."""
    result = runner.invoke(
        main,
        [
            '--loglevel',
            'ERROR',
            'dist-calc',
            '-a',
            str(data_dir / 'pairSeqA.fa'),
            '-b',
            str(data_dir / 'pairSeqB.fa'),
            '-n',
            str(data_dir / 'pairs.txt'),
            '-o',
            str(tmp_path / 'stats.txt'),
            '--no-fig',
            '--proc',
            '1',
        ],
    )

    assert result.exit_code == 0, result.output
    assert 'minimum tolerable alignment score' in result.output
    assert '--score-min L,' in result.output
    assert 'your current setting' in result.output


def test_dist_calc_without_pairs_is_a_clean_error(
    runner: CliRunner, tmp_path: Path, data_dir: Path
) -> None:
    """A missing pair source is reported as a message, not a traceback."""
    result = runner.invoke(
        main,
        [
            'dist-calc',
            '-a',
            str(data_dir / 'pairSeqA.fa'),
            '-b',
            str(data_dir / 'pairSeqB.fa'),
            '-o',
            str(tmp_path / 'stats.txt'),
            '--no-fig',
        ],
    )

    assert result.exit_code == 1
    assert 'Provide a list of transcript pairs' in result.output
    assert 'Traceback' not in result.output


def test_edger_summary_label_mismatch_is_rejected(
    runner: CliRunner, tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """Mismatched --labels and --seqFiles counts fail before any work."""
    toptags = tmp_path / 'toptags.csv'
    toptags.write_text(
        ',logFC,logCPM,LR,PValue,FDR\nCluster-46274.0,3.5,8.1,40.2,1e-10,1e-8\n'
    )

    result = runner.invoke(
        main,
        [
            'edger-summary',
            '-i',
            str(toptags),
            '--clust',
            str(cluster_map),
            '--seqFiles',
            str(data_dir / 'transcripts_X.fa'),
            '--labels',
            'SetX',
            '--labels',
            'SetY',
            '-o',
            str(tmp_path / 'reports'),
        ],
    )

    assert result.exit_code != 0
    assert '2 labels for 1 sequence files' in result.output


def test_edger_summary_command(
    runner: CliRunner, tmp_path: Path, data_dir: Path, cluster_map: Path
) -> None:
    """edger-summary runs end to end and reports where it wrote things."""
    toptags = tmp_path / 'toptags.csv'
    toptags.write_text(
        ',logFC,logCPM,LR,PValue,FDR\n'
        'Cluster-46274.0,3.5,8.1,40.2,1e-10,1e-8\n'
        'Cluster-31069.0,-2.1,7.0,20.5,1e-5,0.001\n'
    )

    result = runner.invoke(
        main,
        [
            '--loglevel',
            'ERROR',
            'edger-summary',
            '-i',
            str(toptags),
            '--clust',
            str(cluster_map),
            '--seqFiles',
            str(data_dir / 'transcripts_X.fa'),
            '--seqFiles',
            str(data_dir / 'transcripts_Y.fa'),
            '--labels',
            'SetX',
            '--labels',
            'SetY',
            '-o',
            str(tmp_path / 'reports'),
        ],
    )

    assert result.exit_code == 0, result.output
    assert 'Reported 2 significant clusters' in result.output
    assert (tmp_path / 'reports' / 'toptags_significant_cluster_report.tab').is_file()
