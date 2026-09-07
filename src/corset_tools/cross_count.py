"""
Count cluster members contributed by each of two transcriptomes.

When Corset is run over reads from two different transcriptome assemblies, a
cluster that draws all of its members from a single assembly is a warning sign:
the two assemblies' copies of the same gene failed to cross-map and were split
into separate clusters, which shows up downstream as spurious differential
expression. This module tallies, per cluster, how many members came from each
input transcriptome.
"""

from dataclasses import dataclass
import logging
import os
from typing import Optional

from .fileio import (
    PathLike,
    label_transcripts,
    read_cluster_map,
    read_fasta,
)

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ClusterCounts:
    """
    Per-cluster membership counts.

    Attributes
    ----------
    cluster : str
        Cluster identifier.
    count_x : int
        Members originating from transcriptome X.
    count_y : int
        Members originating from transcriptome Y.
    total : int
        Total members listed for the cluster, including any whose source
        transcriptome could not be determined.
    """

    cluster: str
    count_x: int
    count_y: int
    total: int

    @property
    def unassigned(self) -> int:
        """
        Count members whose source transcriptome was not found.

        Returns
        -------
        int
            ``total`` minus the members assigned to X or Y.
        """
        return self.total - self.count_x - self.count_y


@dataclass(frozen=True)
class CrossCountSummary:
    """
    Whole-run summary of single-source clusters.

    Attributes
    ----------
    label_x : str
        Resolved label used for transcriptome X.
    label_y : str
        Resolved label used for transcriptome Y.
    n_clusters : int
        Number of clusters examined.
    zero_x : int
        Clusters with no members from transcriptome X.
    zero_y : int
        Clusters with no members from transcriptome Y.
    zero_both : int
        Clusters with no members from either transcriptome.
    """

    label_x: str
    label_y: str
    n_clusters: int
    zero_x: int
    zero_y: int
    zero_both: int


def count_cluster_members(
    clusters: dict[str, list[str]],
    transcript_labels: dict[str, str],
    name_x: str,
    name_y: str,
) -> list[ClusterCounts]:
    """
    Tally per-cluster membership by source transcriptome.

    Parameters
    ----------
    clusters : dict of str to list of str
        Mapping of cluster identifier to member transcript names.
    transcript_labels : dict of str to str
        Mapping of transcript name to its source label.
    name_x : str
        Label used for transcriptome X.
    name_y : str
        Label used for transcriptome Y.

    Returns
    -------
    list of ClusterCounts
        One entry per cluster, in cluster-map order.
    """
    counts: list[ClusterCounts] = []
    for cluster, members in clusters.items():
        count_x = 0
        count_y = 0
        for transcript in members:
            label = transcript_labels.get(transcript)
            if label is None:
                logger.warning('No set name found for transcript: %s', transcript)
            elif label == name_x:
                count_x += 1
            elif label == name_y:
                count_y += 1
        counts.append(
            ClusterCounts(
                cluster=cluster,
                count_x=count_x,
                count_y=count_y,
                total=len(members),
            )
        )
    return counts


def summarise_counts(
    counts: list[ClusterCounts], label_x: str, label_y: str
) -> CrossCountSummary:
    """
    Summarise how many clusters lack members from one or both transcriptomes.

    Parameters
    ----------
    counts : list of ClusterCounts
        Per-cluster counts from :func:`count_cluster_members`.
    label_x : str
        Label used for transcriptome X.
    label_y : str
        Label used for transcriptome Y.

    Returns
    -------
    CrossCountSummary
        Aggregate counts of single-source clusters.
    """
    zero_x = sum(1 for c in counts if c.count_x == 0)
    zero_y = sum(1 for c in counts if c.count_y == 0)
    zero_both = sum(1 for c in counts if c.count_x == 0 and c.count_y == 0)
    return CrossCountSummary(
        label_x=label_x,
        label_y=label_y,
        n_clusters=len(counts),
        zero_x=zero_x,
        zero_y=zero_y,
        zero_both=zero_both,
    )


def _default_label(path: PathLike) -> str:
    """
    Derive a set label from a FASTA file name.

    Parameters
    ----------
    path : str or os.PathLike
        FASTA path.

    Returns
    -------
    str
        The file's base name with its extension removed.
    """
    return os.path.basename(os.path.splitext(str(path))[0])


def write_counts(
    path: PathLike,
    counts: list[ClusterCounts],
    name_x: str,
    name_y: str,
) -> None:
    """
    Write the per-cluster count table as a TSV.

    Parameters
    ----------
    path : str or os.PathLike
        Output file path.
    counts : list of ClusterCounts
        Per-cluster counts to write.
    name_x : str
        Column heading for transcriptome X.
    name_y : str
        Column heading for transcriptome Y.

    Returns
    -------
    None
        The table is written to ``path``.
    """
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(['clusterID', name_x, name_y, 'Total_Members']) + '\n')
        for row in counts:
            handle.write(
                '\t'.join(
                    [row.cluster, str(row.count_x), str(row.count_y), str(row.total)]
                )
                + '\n'
            )


def run_cross_count(
    fasta_x: PathLike,
    fasta_y: PathLike,
    cluster_map: PathLike,
    name_x: Optional[str] = None,
    name_y: Optional[str] = None,
    out_file: PathLike = 'CountClusterMembers.txt',
) -> tuple[list[ClusterCounts], CrossCountSummary]:
    """
    Count cluster members from two transcriptomes and write the count table.

    Parameters
    ----------
    fasta_x : str or os.PathLike
        Transcriptome X FASTA.
    fasta_y : str or os.PathLike
        Transcriptome Y FASTA.
    cluster_map : str or os.PathLike
        Corset transcript-to-cluster map.
    name_x : str, optional
        Label for transcriptome X. Defaults to the FASTA base name.
    name_y : str, optional
        Label for transcriptome Y. Defaults to the FASTA base name.
    out_file : str or os.PathLike, optional
        Output TSV path, by default ``CountClusterMembers.txt``.

    Returns
    -------
    tuple of (list of ClusterCounts, CrossCountSummary)
        The per-cluster counts and the aggregate summary.

    Raises
    ------
    ValueError
        If both transcriptomes resolve to the same label, which would make the
        two count columns indistinguishable.
    """
    label_x = name_x or _default_label(fasta_x)
    label_y = name_y or _default_label(fasta_y)
    if label_x == label_y:
        raise ValueError(
            f'Transcriptome labels must differ, both resolved to {label_x!r}. '
            'Set them explicitly with -x/-y.'
        )

    transcript_labels = label_transcripts(
        {label_x: read_fasta(fasta_x), label_y: read_fasta(fasta_y)}
    )
    clusters = read_cluster_map(cluster_map)

    counts = count_cluster_members(clusters, transcript_labels, label_x, label_y)
    summary = summarise_counts(counts, label_x, label_y)
    write_counts(out_file, counts, label_x, label_y)

    logger.info('Wrote counts for %d clusters to %s', len(counts), out_file)
    return counts, summary
