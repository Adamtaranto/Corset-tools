"""
Command line interface for corset-tools.

Exposes the four helper tools as subcommands of a single ``corset-tools``
entry point. Each subcommand is a thin wrapper that validates arguments and
delegates to a ``run_*`` function, so the underlying logic stays importable and
testable without click.
"""

import logging
from typing import Optional

import click

from . import __version__
from .cross_count import run_cross_count
from .dist_calc.runner import run_dist_calc
from .dist_calc.scoring import format_advice
from .edger_summary import run_edger_summary
from .exceptions import CorsetToolsError
from .fetch_seqs import not_found_path, run_fetch_seqs
from .fileio import read_target_list
from .logs import init_logging

logger = logging.getLogger(__name__)

LOG_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']


class CorsetToolsGroup(click.Group):
    """Click group that reports package errors as clean CLI messages."""

    def invoke(self, ctx: click.Context) -> object:
        """
        Invoke a subcommand, translating package errors for the user.

        Parameters
        ----------
        ctx : click.Context
            The active click context.

        Returns
        -------
        object
            Whatever the invoked subcommand returns.

        Raises
        ------
        click.ClickException
            Wrapping any :class:`~corset_tools.exceptions.CorsetToolsError`, so
            the user sees a one-line message rather than a traceback.
        """
        try:
            return super().invoke(ctx)
        except CorsetToolsError as error:
            raise click.ClickException(str(error)) from error
        except ValueError as error:
            # The run_* functions signal bad argument combinations this way.
            raise click.ClickException(str(error)) from error


@click.group(
    cls=CorsetToolsGroup,
    context_settings={'help_option_names': ['-h', '--help']},
    help='Companion tools for working with transcript clusters from Corset.',
)
@click.version_option(__version__, '-V', '--version', prog_name='corset-tools')
@click.option(
    '--loglevel',
    type=click.Choice(LOG_LEVELS, case_sensitive=False),
    default='INFO',
    show_default=True,
    help='Console logging verbosity.',
)
@click.option(
    '--logfile',
    type=click.Path(dir_okay=False),
    default=None,
    help='Also write logs to this file.',
)
def main(loglevel: str, logfile: Optional[str]) -> None:
    """
    Initialise logging for whichever subcommand runs.

    Parameters
    ----------
    loglevel : str
        Console logging verbosity.
    logfile : str, optional
        Path to an additional log file.

    Returns
    -------
    None
        Logging is configured as a side effect.
    """
    init_logging(loglevel=loglevel, logfile=logfile)


@main.command('fetch-seqs')
@click.option(
    '-i',
    '--inFasta',
    'in_fasta',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Multi-fasta of transcript sequences.',
)
@click.option(
    '-c',
    '--clustMap',
    'cluster_map',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Corset cluster map: 'transcript_name Cluster_ID'.",
)
@click.option(
    '-t',
    '--targetClust',
    'target_clust',
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help='List of target cluster names, one per line. Default: all clusters.',
)
@click.option(
    '-o',
    '--outFasta',
    'out_fasta',
    default='filtered_seqs.fa',
    show_default=True,
    help='Output fasta path.',
)
@click.option(
    '-l',
    '--longest',
    is_flag=True,
    default=False,
    help='Report only the longest transcript of each cluster.',
)
@click.option(
    '--notFoundLog',
    'not_found_log',
    default=None,
    help='Log path for unresolved clusters and transcripts. '
    'Default: NotFound_<target list>.log',
)
@click.option(
    '--width',
    default=80,
    show_default=True,
    type=click.IntRange(min=1),
    help='Fasta line wrap width.',
)
def fetch_seqs_cmd(
    in_fasta: str,
    cluster_map: str,
    target_clust: Optional[str],
    out_fasta: str,
    longest: bool,
    not_found_log: Optional[str],
    width: int,
) -> None:
    """Extract cluster member transcripts, tagged with their cluster ID."""
    targets = read_target_list(target_clust) if target_clust else None
    if not_found_log is None:
        not_found_log = not_found_path(target_clust)

    result = run_fetch_seqs(
        in_fasta=in_fasta,
        cluster_map=cluster_map,
        targets=targets,
        out_fasta=out_fasta,
        longest=longest,
        not_found_log=not_found_log,
        width=width,
    )
    click.echo(f'Wrote {result.written} sequences to {out_fasta}')
    if result.missing_clusters:
        click.echo(
            f'{len(result.missing_clusters)} target clusters were not in the map'
        )
    if result.missing_transcripts:
        click.echo(
            f'{len(result.missing_transcripts)} cluster members were not in the fasta'
        )


@main.command('cross-count')
@click.option(
    '-X',
    '--transFastaX',
    'fasta_x',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Transcriptome X fasta.',
)
@click.option(
    '-Y',
    '--transFastaY',
    'fasta_y',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Transcriptome Y fasta.',
)
@click.option(
    '-c',
    '--clustMap',
    'cluster_map',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Corset cluster map: 'transcript_name Cluster_ID'.",
)
@click.option(
    '-x', '--transNameX', 'name_x', default=None, help='Label for transcriptome X.'
)
@click.option(
    '-y', '--transNameY', 'name_y', default=None, help='Label for transcriptome Y.'
)
@click.option(
    '-o',
    '--outFile',
    'out_file',
    default='CountClusterMembers.txt',
    show_default=True,
    help='Output count table.',
)
def cross_count_cmd(
    fasta_x: str,
    fasta_y: str,
    cluster_map: str,
    name_x: Optional[str],
    name_y: Optional[str],
    out_file: str,
) -> None:
    """Count cluster members contributed by each of two transcriptomes."""
    counts, summary = run_cross_count(
        fasta_x=fasta_x,
        fasta_y=fasta_y,
        cluster_map=cluster_map,
        name_x=name_x,
        name_y=name_y,
        out_file=out_file,
    )
    # Report whatever labels were actually used, including the fasta base names
    # that stand in when -x/-y are omitted.
    label_x = summary.label_x
    label_y = summary.label_y
    click.echo('')
    click.echo('Summary stats:')
    click.echo(f'Clusters examined: {summary.n_clusters}')
    click.echo(f'Clusters with 0 members from {label_x}: {summary.zero_x}')
    click.echo(f'Clusters with 0 members from {label_y}: {summary.zero_y}')
    click.echo(
        f'Clusters with 0 members from either {label_x} or {label_y}: '
        f'{summary.zero_both}'
    )
    click.echo(f'Wrote counts for {len(counts)} clusters to {out_file}')


@main.command('dist-calc')
@click.option(
    '-a',
    '--fastaA',
    'fasta_a',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Transcriptome A fasta.',
)
@click.option(
    '-b',
    '--fastaB',
    'fasta_b',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Transcriptome B fasta.',
)
@click.option(
    '-n',
    '--pairNames',
    'pair_names',
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help='Two-column table of matched transcript pairs.',
)
@click.option(
    '-x',
    '--blastAvB',
    'blast_a_vs_b',
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help='BLAST tabular output, A queried against B.',
)
@click.option(
    '-y',
    '--blastBvA',
    'blast_b_vs_a',
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help='BLAST tabular output, B queried against A.',
)
@click.option(
    '-r',
    '--readLength',
    'read_length',
    default=100,
    show_default=True,
    type=click.IntRange(min=1),
    help='Read length to model.',
)
@click.option(
    '-i',
    '--scoreMinIntercept',
    'score_min_intercept',
    default=-0.6,
    show_default=True,
    type=float,
    help='Intercept of the bowtie2 --score-min function.',
)
@click.option(
    '-s',
    '--scoreMinSlope',
    'score_min_slope',
    default=-0.6,
    show_default=True,
    type=float,
    help='Slope of the bowtie2 --score-min function.',
)
@click.option(
    '-g',
    '--gapOpen',
    'gap_open',
    default=-5.0,
    show_default=True,
    type=float,
    help='Penalty for opening a gap.',
)
@click.option(
    '-e',
    '--gapExtend',
    'gap_extend',
    default=-3.0,
    show_default=True,
    type=float,
    help='Penalty for extending a gap.',
)
@click.option(
    '-m',
    '--mismatch',
    default=-6.0,
    show_default=True,
    type=float,
    help='Penalty for a mismatched position.',
)
@click.option(
    '-p',
    '--percentile',
    default=95.0,
    show_default=True,
    type=click.FloatRange(min=0, max=100, min_open=True),
    help='Target percentage of transcript pairs that should cross-map.',
)
@click.option(
    '-o',
    '--outFile',
    'out_file',
    default='alignmentStats.txt',
    show_default=True,
    help='Output statistics table.',
)
@click.option(
    '-f',
    '--outFig',
    'out_fig',
    default='readPenaltyDist.pdf',
    show_default=True,
    help='Read penalty distribution figure. Use --no-fig to skip.',
)
@click.option(
    '--no-fig', 'no_fig', is_flag=True, default=False, help='Skip the figure.'
)
@click.option(
    '-l',
    '--minLen',
    'min_len',
    default=0,
    show_default=True,
    type=int,
    help='Minimum BLAST hit length to consider valid.',
)
@click.option(
    '-E',
    '--eVal',
    'e_value',
    default=0.001,
    show_default=True,
    type=float,
    help='Maximum BLAST e-value to consider valid.',
)
@click.option(
    '-w',
    '--recipFile',
    'write_pairs_file',
    is_flag=True,
    default=False,
    help='Write reciprocal best-hit pairs to file.',
)
@click.option(
    '-v',
    '--verbose',
    is_flag=True,
    default=False,
    help='Print formatted alignments to screen.',
)
@click.option(
    '--proc',
    'processes',
    default=4,
    show_default=True,
    type=click.IntRange(min=1),
    help='Processes to split the alignment job over.',
)
def dist_calc_cmd(
    fasta_a: str,
    fasta_b: str,
    pair_names: Optional[str],
    blast_a_vs_b: Optional[str],
    blast_b_vs_a: Optional[str],
    read_length: int,
    score_min_intercept: float,
    score_min_slope: float,
    gap_open: float,
    gap_extend: float,
    mismatch: float,
    percentile: float,
    out_file: str,
    out_fig: str,
    no_fig: bool,
    min_len: int,
    e_value: float,
    write_pairs_file: bool,
    verbose: bool,
    processes: int,
) -> None:
    """Estimate whether reads cross-map between two transcriptomes."""
    result = run_dist_calc(
        fasta_a=fasta_a,
        fasta_b=fasta_b,
        pair_names=pair_names,
        blast_a_vs_b=blast_a_vs_b,
        blast_b_vs_a=blast_b_vs_a,
        read_length=read_length,
        score_min_intercept=score_min_intercept,
        score_min_slope=score_min_slope,
        gap_open=gap_open,
        gap_extend=gap_extend,
        mismatch=mismatch,
        percentile=percentile,
        out_file=out_file,
        out_fig=None if no_fig else out_fig,
        processes=processes,
        min_len=min_len,
        e_value=e_value,
        write_pairs_file=write_pairs_file,
        verbose=verbose,
    )

    click.echo('')
    click.echo('Summary stats:')
    click.echo(f'Total pairwise comparisons made: {len(result.stats)}')
    click.echo(f'Median alignment length: {result.medians["len"]:g}')
    click.echo(f'Median gaps: {result.medians["gaps"]:g}')
    click.echo(f'Median mismatch: {result.medians["mismatch"]:g}')
    click.echo(f'Median alignment score: {result.medians["score"]:g}')
    click.echo(f'Median predicted read penalty: {result.medians["readScore"]:g}')
    click.echo('')
    click.echo(
        f'Current --score-min threshold of {result.min_score:g} is predicted to '
        f'allow cross-mapping for {result.cross_map_pass} of {len(result.stats)} pairs.'
    )
    click.echo('')
    click.echo(format_advice(result.advice))


@main.command('edger-summary')
@click.option(
    '-i',
    '--infile',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="edgeR topTags CSV: 'ClusterID,logFC,logCPM,LR,PValue,FDR'.",
)
@click.option(
    '--clust',
    'cluster_map',
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Corset cluster map: 'transcript_name Cluster_ID'.",
)
@click.option(
    '--seqFiles',
    'seq_files',
    required=True,
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Transcript fasta. Repeat once per transcriptome.',
)
@click.option(
    '--labels',
    required=True,
    multiple=True,
    help='Label for each --seqFiles entry, in the same order.',
)
@click.option(
    '--grepme',
    'grep_files',
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help='Annotation file to harvest lines from. Repeatable.',
)
@click.option(
    '-o',
    '--outDir',
    'out_dir',
    default='cluster_reports',
    show_default=True,
    help='Output directory.',
)
@click.option(
    '--maxfdr',
    'max_fdr',
    default=0.05,
    show_default=True,
    type=float,
    help='Only records at or below this FDR are considered.',
)
@click.option(
    '--minfc',
    'min_fc',
    default=1.0,
    show_default=True,
    type=float,
    help='Only records at or above this absolute log2 fold change are considered.',
)
@click.option(
    '-x',
    '--topx',
    default=None,
    type=int,
    help='Keep only the top x clusters, sorted by FDR.',
)
@click.option(
    '--writeClusters',
    'write_clusters',
    is_flag=True,
    default=False,
    help='Also write member transcripts to one fasta per cluster.',
)
def edger_summary_cmd(
    infile: str,
    cluster_map: str,
    seq_files: tuple[str, ...],
    labels: tuple[str, ...],
    grep_files: tuple[str, ...],
    out_dir: str,
    max_fdr: float,
    min_fc: float,
    topx: Optional[int],
    write_clusters: bool,
) -> None:
    """Summarise and annotate differentially expressed clusters from edgeR."""
    if len(labels) != len(seq_files):
        raise click.BadParameter(
            f'Got {len(labels)} labels for {len(seq_files)} sequence files.',
            param_hint='--labels',
        )

    result = run_edger_summary(
        infile=infile,
        cluster_map=cluster_map,
        seq_files=list(seq_files),
        labels=list(labels),
        grep_files=list(grep_files),
        out_dir=out_dir,
        max_fdr=max_fdr,
        min_fc=min_fc,
        topx=topx,
        write_clusters=write_clusters,
    )
    click.echo(f'Reported {len(result.table)} significant clusters')
    click.echo(f'Cluster report: {result.report_path}')
    for path in result.annotation_paths:
        click.echo(f'Annotation records: {path}')
    if result.fasta_paths:
        click.echo(f'Wrote {len(result.fasta_paths)} cluster fasta files')


if __name__ == '__main__':  # pragma: no cover
    main()
