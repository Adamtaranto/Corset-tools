"""
Summarise and annotate differentially expressed clusters reported by edgeR.

Takes an edgeR ``topTags`` table of differentially expressed Corset clusters,
filters it by FDR / log fold-change / rank, annotates each surviving cluster
with the number of members contributed by each input transcriptome, and
harvests matching lines from arbitrary annotation files. Optionally writes one
FASTA per significant cluster.
"""

from dataclasses import dataclass
import logging
import os
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from .fileio import (
    PathLike,
    label_transcripts,
    read_cluster_map,
    read_multi_fasta,
    write_fasta,
)

logger = logging.getLogger(__name__)

# edgeR topTags exports carry a nameless index column followed by the test
# statistics; rename positionally so differing edgeR versions still line up.
DE_COLUMNS = ['ClusterID', 'logFC', 'logCPM', 'LR', 'PValue', 'FDR']

# The original script wrote per-cluster FASTAs at 60 columns.
DEFAULT_LINE_WIDTH = 60


@dataclass
class EdgerSummaryResult:
    """
    Paths and counts produced by :func:`run_edger_summary`.

    Attributes
    ----------
    table : pandas.DataFrame
        The filtered DE table with membership count columns appended.
    report_path : pathlib.Path
        Path of the written cluster report TSV.
    annotation_paths : list of pathlib.Path
        One annotation record file per input annotation file.
    fasta_paths : list of pathlib.Path
        Per-cluster FASTA files, empty unless ``write_clusters`` was set.
    """

    table: pd.DataFrame
    report_path: Path
    annotation_paths: list[Path]
    fasta_paths: list[Path]


def filter_de_table(
    infile: PathLike,
    max_fdr: float = 0.05,
    min_fc: float = 1.0,
    topx: Optional[int] = None,
) -> pd.DataFrame:
    """
    Load and filter an edgeR topTags table.

    Parameters
    ----------
    infile : str or os.PathLike
        CSV exported from edgeR ``topTags``, with columns corresponding to
        ``ClusterID,logFC,logCPM,LR,PValue,FDR``.
    max_fdr : float, optional
        Keep only rows with FDR at or below this value, by default 0.05.
    min_fc : float, optional
        Keep only rows whose absolute log2 fold change is at least this value,
        by default 1.0 (a two-fold change).
    topx : int, optional
        If given, keep only the top ``x`` rows after sorting by FDR.

    Returns
    -------
    pandas.DataFrame
        The filtered table, sorted by ascending FDR with a fresh index.

    Raises
    ------
    ValueError
        If the input table does not have the expected six columns.
    """
    table = pd.read_csv(infile)
    if table.shape[1] != len(DE_COLUMNS):
        raise ValueError(
            f'Expected {len(DE_COLUMNS)} columns in {infile} '
            f'({",".join(DE_COLUMNS)}), found {table.shape[1]}.'
        )
    table.columns = DE_COLUMNS
    table['ClusterID'] = table['ClusterID'].astype(str)

    table = table[(table.logFC.abs() >= min_fc) & (table.FDR <= max_fdr)]
    table = table.sort_values(by='FDR', ascending=True).reset_index(drop=True)
    if topx:
        table = table.head(topx)

    logger.info('Retained %d differentially expressed clusters', len(table))
    return table


def summarise_membership(
    clusters: dict[str, list[str]],
    transcript_labels: dict[str, str],
    query_clusters: Sequence[str],
    labels: Sequence[str],
) -> dict[str, dict[str, int]]:
    """
    Count members of each query cluster by source transcriptome.

    Parameters
    ----------
    clusters : dict of str to list of str
        Cluster-to-members map.
    transcript_labels : dict of str to str
        Transcript-to-source-label map.
    query_clusters : sequence of str
        Clusters to summarise.
    labels : sequence of str
        Source labels, used so every cluster reports every label (zeros
        included).

    Returns
    -------
    dict of str to (dict of str to int)
        Mapping of cluster to per-label member counts.
    """
    summary: dict[str, dict[str, int]] = {}
    for cluster in query_clusters:
        counts = dict.fromkeys(labels, 0)
        members = clusters.get(cluster)
        if members is None:
            logger.warning('DE cluster %s not present in the cluster map', cluster)
        else:
            for transcript in members:
                label = transcript_labels.get(transcript)
                if label is None:
                    logger.warning(
                        'No transcriptome label for transcript %s in cluster %s',
                        transcript,
                        cluster,
                    )
                else:
                    counts[label] += 1
        summary[cluster] = counts
    return summary


def add_membership_columns(
    table: pd.DataFrame,
    membership: dict[str, dict[str, int]],
    labels: Sequence[str],
) -> pd.DataFrame:
    """
    Append per-transcriptome member counts and a total to the DE table.

    Parameters
    ----------
    table : pandas.DataFrame
        Filtered DE table containing a ``ClusterID`` column.
    membership : dict of str to (dict of str to int)
        Per-cluster counts from :func:`summarise_membership`.
    labels : sequence of str
        Source labels, in the order the columns should appear.

    Returns
    -------
    pandas.DataFrame
        A new table with one column per label plus a ``Totals`` column.
    """
    counts = pd.DataFrame.from_dict(membership, orient='index', dtype=int)
    # Guarantee column presence and order even if membership is empty.
    counts = counts.reindex(columns=list(labels), fill_value=0)
    counts['Totals'] = counts[list(labels)].sum(axis=1)
    counts.index.name = 'ClusterID'
    counts = counts.reset_index()
    return table.merge(counts, how='left', on='ClusterID')


def write_annotation_records(
    clusters: dict[str, list[str]],
    annotation_file: PathLike,
    out_file: PathLike,
    query_clusters: Sequence[str],
) -> Path:
    """
    Collect annotation lines mentioning members of the query clusters.

    Parameters
    ----------
    clusters : dict of str to list of str
        Cluster-to-members map.
    annotation_file : str or os.PathLike
        Free-form annotation file to scan for transcript names.
    out_file : str or os.PathLike
        Destination for the harvested records.
    query_clusters : sequence of str
        Clusters to report, in output order.

    Returns
    -------
    pathlib.Path
        Path of the written file.
    """
    # Only members of the query clusters are of interest. A transcript can only
    # belong to one cluster, so a flat name -> cluster lookup is enough.
    wanted: dict[str, str] = {}
    for cluster in query_clusters:
        for transcript in clusters.get(cluster, []):
            wanted[transcript] = cluster

    # Bucket by cluster in a single pass over the annotation file. The original
    # re-opened and re-read the whole file once per cluster.
    buckets: dict[str, list[str]] = {cluster: [] for cluster in query_clusters}
    names = list(wanted)
    with open(annotation_file, encoding='utf-8') as handle:
        for line in handle:
            for name in names:
                if name in line:
                    buckets[wanted[name]].append(line)

    out_path = Path(out_file)
    with open(out_path, 'w', encoding='utf-8') as handle:
        for cluster in query_clusters:
            handle.write(f'#{cluster}\n')
            handle.writelines(buckets[cluster])

    logger.info('Wrote annotation records to %s', out_path)
    return out_path


def write_cluster_fastas(
    out_dir: PathLike,
    sequence_sets: dict[str, dict[str, str]],
    clusters: dict[str, list[str]],
    query_clusters: Sequence[str],
    width: int = DEFAULT_LINE_WIDTH,
) -> list[Path]:
    r"""
    Write one FASTA of member sequences per query cluster.

    Record names are ``<transcript>\\t<set label>\\t<cluster>``, matching the
    original script's output.

    Parameters
    ----------
    out_dir : str or os.PathLike
        Directory to write the per-cluster FASTA files into.
    sequence_sets : dict of str to (dict of str to str)
        Mapping of set label to that set's sequences.
    clusters : dict of str to list of str
        Cluster-to-members map.
    query_clusters : sequence of str
        Clusters to write.
    width : int, optional
        FASTA line wrap width, by default 60.

    Returns
    -------
    list of pathlib.Path
        Paths of the FASTA files written.
    """
    written: list[Path] = []
    for cluster in query_clusters:
        out_path = Path(out_dir) / f'{cluster}.fa'
        if out_path.is_file():
            logger.info('Cluster fasta already exists, skipping: %s', out_path)
            continue
        members = clusters.get(cluster, [])
        with open(out_path, 'w', encoding='utf-8') as handle:
            for label, sequences in sequence_sets.items():
                for transcript in members:
                    if transcript in sequences:
                        write_fasta(
                            handle,
                            '\t'.join((transcript, label, cluster)),
                            sequences[transcript],
                            width,
                        )
        written.append(out_path)
    return written


def run_edger_summary(
    infile: PathLike,
    cluster_map: PathLike,
    seq_files: Sequence[PathLike],
    labels: Sequence[str],
    grep_files: Sequence[PathLike] = (),
    out_dir: PathLike = 'cluster_reports',
    max_fdr: float = 0.05,
    min_fc: float = 1.0,
    topx: Optional[int] = None,
    write_clusters: bool = False,
) -> EdgerSummaryResult:
    """
    Filter, annotate and report differentially expressed clusters.

    Parameters
    ----------
    infile : str or os.PathLike
        The edgeR topTags CSV.
    cluster_map : str or os.PathLike
        Corset transcript-to-cluster map.
    seq_files : sequence of str or os.PathLike
        Transcript FASTA files, one per source transcriptome.
    labels : sequence of str
        Labels for ``seq_files``, in the same order.
    grep_files : sequence of str or os.PathLike, optional
        Annotation files to harvest lines from.
    out_dir : str or os.PathLike, optional
        Output directory, by default ``cluster_reports``.
    max_fdr : float, optional
        FDR ceiling, by default 0.05.
    min_fc : float, optional
        Absolute log2 fold-change floor, by default 1.0.
    topx : int, optional
        Keep only the top ``x`` clusters by FDR.
    write_clusters : bool, optional
        If True, also write one FASTA of member sequences per cluster.

    Returns
    -------
    EdgerSummaryResult
        The annotated table and the paths of every file written.

    Raises
    ------
    ValueError
        If ``labels`` and ``seq_files`` differ in length.
    """
    if len(labels) != len(seq_files):
        raise ValueError(
            f'Got {len(labels)} labels for {len(seq_files)} sequence files; '
            'these must match.'
        )

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    clusters = read_cluster_map(cluster_map)
    sequence_sets = read_multi_fasta(list(seq_files), list(labels))
    transcript_labels = label_transcripts(sequence_sets)

    table = filter_de_table(infile, max_fdr=max_fdr, min_fc=min_fc, topx=topx)
    query_clusters = list(table['ClusterID'])

    membership = summarise_membership(
        clusters, transcript_labels, query_clusters, labels
    )
    table = add_membership_columns(table, membership, labels)

    infile_base = os.path.splitext(os.path.basename(str(infile)))[0]

    report_path = out_path / f'{infile_base}_significant_cluster_report.tab'
    if report_path.is_file():
        logger.info('Overwriting existing report file: %s', report_path)
    table.to_csv(report_path, sep='\t', index=False, encoding='utf-8')

    annotation_paths: list[Path] = []
    for annotation_file in grep_files:
        annot_base = os.path.splitext(os.path.basename(str(annotation_file)))[0]
        annot_out = (
            out_path
            / f'{annot_base}_{infile_base}_significant_cluster_annotation_records.txt'
        )
        annotation_paths.append(
            write_annotation_records(
                clusters, annotation_file, annot_out, query_clusters
            )
        )

    fasta_paths: list[Path] = []
    if write_clusters:
        fasta_paths = write_cluster_fastas(
            out_path, sequence_sets, clusters, query_clusters
        )

    return EdgerSummaryResult(
        table=table,
        report_path=report_path,
        annotation_paths=annotation_paths,
        fasta_paths=fasta_paths,
    )
