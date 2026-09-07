"""
Extract cluster member transcripts from a multi-FASTA.

Given a Corset transcript-to-cluster map and an optional list of target
clusters, write out the transcript sequences belonging to those clusters,
renamed as ``<clusterID>_<transcriptID>`` so that cluster membership survives
into downstream annotation.
"""

from dataclasses import dataclass, field
import logging
import os
from pathlib import Path
from typing import Optional

from .fileio import PathLike, read_cluster_map, read_fasta, write_fasta

logger = logging.getLogger(__name__)

# The original script wrapped output at 80 columns; keep that as the default so
# existing downstream files stay byte-comparable.
DEFAULT_LINE_WIDTH = 80


@dataclass
class FetchResult:
    """
    Summary of a :func:`run_fetch_seqs` call.

    Attributes
    ----------
    written : int
        Number of FASTA records written.
    missing_clusters : list of str
        Target clusters that were absent from the cluster map.
    missing_transcripts : list of (str, str)
        ``(cluster, transcript)`` pairs listed in the map but absent from the
        input FASTA.
    empty_clusters : list of str
        Clusters present in the map for which no member sequence was found.
    """

    written: int = 0
    missing_clusters: list[str] = field(default_factory=list)
    missing_transcripts: list[tuple[str, str]] = field(default_factory=list)
    empty_clusters: list[str] = field(default_factory=list)


def not_found_path(target_clusters: Optional[PathLike]) -> str:
    """
    Build the name of the not-found log for a given target list.

    Parameters
    ----------
    target_clusters : str or os.PathLike, optional
        Path to the target cluster list, or None when all clusters are used.

    Returns
    -------
    str
        File name for the not-found log.
    """
    if target_clusters is None:
        return 'NotFound_Clusters_Transcripts.log'
    return f'NotFound_{os.path.basename(str(target_clusters))}.log'


def run_fetch_seqs(
    in_fasta: PathLike,
    cluster_map: PathLike,
    targets: Optional[list[str]] = None,
    out_fasta: PathLike = 'filtered_seqs.fa',
    longest: bool = False,
    not_found_log: Optional[PathLike] = None,
    width: int = DEFAULT_LINE_WIDTH,
) -> FetchResult:
    """
    Write cluster-tagged transcript sequences for the requested clusters.

    Parameters
    ----------
    in_fasta : str or os.PathLike
        Multi-FASTA of transcript sequences.
    cluster_map : str or os.PathLike
        Corset transcript-to-cluster map.
    targets : list of str, optional
        Cluster names to report. If None, every cluster in the map is reported.
    out_fasta : str or os.PathLike, optional
        Output FASTA path, by default ``filtered_seqs.fa``.
    longest : bool, optional
        If True, emit only the longest transcript of each cluster.
    not_found_log : str or os.PathLike, optional
        Path for the log of clusters/transcripts that could not be resolved.
        Defaults to ``NotFound_Clusters_Transcripts.log``.
    width : int, optional
        FASTA line wrap width, by default 80.

    Returns
    -------
    FetchResult
        Counts and the lists of unresolved clusters and transcripts.
    """
    clusters = read_cluster_map(cluster_map)
    sequences = read_fasta(in_fasta)

    # No target list means "report everything", preserving map order.
    wanted = list(clusters) if targets is None else targets
    result = FetchResult()

    log_path = Path(not_found_log or 'NotFound_Clusters_Transcripts.log')

    with (
        open(out_fasta, 'w', encoding='utf-8') as fasta_handle,
        open(log_path, 'w', encoding='utf-8') as log_handle,
    ):
        for cluster in wanted:
            if cluster not in clusters:
                logger.warning('Target cluster not in map file: %s', cluster)
                result.missing_clusters.append(cluster)
                log_handle.write(f'{cluster}\n')
                continue

            # Collect the members that actually have a sequence, so the
            # "longest" branch never has to reach for an unbound variable.
            found: list[tuple[str, str]] = []
            for transcript in clusters[cluster]:
                if transcript in sequences:
                    found.append((transcript, sequences[transcript]))
                else:
                    logger.warning(
                        'Transcript not in reference fasta: %s: %s',
                        cluster,
                        transcript,
                    )
                    result.missing_transcripts.append((cluster, transcript))
                    log_handle.write(f'{cluster}\t{transcript}\n')

            if not found:
                # The original script raised UnboundLocalError here, or silently
                # re-emitted the previous cluster's sequence.
                logger.warning('No member sequences found for cluster: %s', cluster)
                result.empty_clusters.append(cluster)
                continue

            if longest:
                # max() with a >= tie-break equivalent: later members of equal
                # length won in the original, so scan in order and keep the last
                # maximum.
                best = found[0]
                for candidate in found[1:]:
                    if len(candidate[1]) >= len(best[1]):
                        best = candidate
                found = [best]

            for transcript, sequence in found:
                write_fasta(fasta_handle, f'{cluster}_{transcript}', sequence, width)
                result.written += 1

    logger.info(
        'Wrote %d sequences to %s (%d missing clusters, %d missing transcripts)',
        result.written,
        out_fasta,
        len(result.missing_clusters),
        len(result.missing_transcripts),
    )
    return result
