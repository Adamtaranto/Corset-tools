"""Distribution plot of predicted per-read alignment penalties."""

import logging
from typing import Optional, Sequence

import matplotlib

# Select a non-interactive backend before pyplot is imported so the command
# works over SSH and in CI, where no display is available.
matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from ..fileio import PathLike  # noqa: E402

logger = logging.getLogger(__name__)


def plot_penalty_distribution(
    read_penalties: Sequence[float],
    out_fig: PathLike,
    threshold: Optional[float] = None,
    bins: int = 20,
) -> None:
    """
    Plot the distribution of predicted per-read alignment penalties.

    Draws a histogram of the penalties with a fitted normal curve, and marks
    the recommended score-minimum threshold with a dashed vertical line.

    Parameters
    ----------
    read_penalties : sequence of float
        Expected per-read penalties, one per aligned transcript pair.
    out_fig : str or os.PathLike
        Path to write the figure to. The format is inferred from the suffix.
    threshold : float, optional
        Score-minimum threshold to mark on the plot.
    bins : int, optional
        Number of histogram bins, by default 20.

    Returns
    -------
    None
        The figure is written to ``out_fig``.

    Raises
    ------
    ValueError
        If ``read_penalties`` is empty.
    """
    if not len(read_penalties):
        raise ValueError('Cannot plot a penalty distribution with no alignments.')

    values = np.sort(np.asarray(read_penalties, dtype=float))

    figure, axes = plt.subplots()
    axes.hist(values, bins=bins, color='c', density=True)

    # A normal fit needs spread; a degenerate distribution would divide by zero.
    spread = float(np.std(values))
    if spread > 0:
        axes.plot(values, stats.norm.pdf(values, float(np.mean(values)), spread), 'b-')

    if threshold is not None:
        axes.axvline(threshold, color='b', linestyle='dashed', linewidth=2)

    axes.set_xlabel('Predicted per-read alignment penalty')
    axes.set_ylabel('Density')
    axes.set_title('Read penalty distribution')

    figure.savefig(out_fig, bbox_inches='tight')
    plt.close(figure)
    logger.info('Wrote penalty distribution figure to %s', out_fig)
