"""Tests for read-penalty scoring and score-min recommendations."""

import numpy as np
import pytest

from corset_tools.dist_calc.scoring import (
    format_advice,
    read_penalty,
    recommended_score_min,
    score_min,
)


def test_score_min_is_linear() -> None:
    """score_min evaluates bowtie2's L,intercept,slope function."""
    # The bowtie2 default of -L,-0.6,-0.6 at 100 bp gives -60.6.
    assert score_min(-0.6, -0.6, 100) == pytest.approx(-60.6)


def test_read_penalty_of_a_perfect_alignment_is_zero() -> None:
    """With no gaps or mismatches there is nothing to penalise."""
    assert read_penalty(0, 0, 100, -5, -3, -6, 100) == 0


def test_read_penalty_scales_with_read_length() -> None:
    """Densities are scaled to the read length, so a longer read costs more."""
    short = read_penalty(1, 2, 100, -5, -3, -6, 50)
    long = read_penalty(1, 2, 100, -5, -3, -6, 100)

    assert long == pytest.approx(2 * short)


def test_read_penalty_arithmetic() -> None:
    """The penalty combines mismatch and gap terms as bowtie2 would score them."""
    # 2 mismatches and 1 gap over 100 columns, modelled for a 100 bp read:
    # 2 * -6 (mismatch) + 1 * (-5 + -3) (gap open + extend).
    assert read_penalty(1, 2, 100, -5, -3, -6, 100) == pytest.approx(-20)


def test_read_penalty_rejects_zero_length() -> None:
    """A zero-length alignment has no density to scale."""
    with pytest.raises(ValueError, match='greater than zero'):
        read_penalty(0, 0, 0, -5, -3, -6, 100)


def test_recommended_score_min_covers_the_target_fraction() -> None:
    """The recommendation always achieves at least the requested coverage."""
    rng = np.random.default_rng(0)
    penalties = rng.normal(-30, 5, 1000)

    advice = recommended_score_min(
        penalties, percentile=95, read_length=100, slope=-0.6
    )

    assert advice.achieved_fraction >= 0.95
    # It should not be wildly over-permissive either.
    assert advice.achieved_fraction < 0.97
    assert advice.n_pairs == 1000


def test_recommended_score_min_intercept_round_trips() -> None:
    """Feeding the recommended intercept back into score_min recovers the threshold."""
    penalties = [-10.0, -20.0, -30.0, -40.0, -50.0]

    advice = recommended_score_min(
        penalties, percentile=95, read_length=100, slope=-0.6
    )

    assert score_min(advice.intercept, advice.slope, 100) == pytest.approx(
        advice.threshold
    )


def test_recommended_score_min_at_100_percent_uses_the_worst_pair() -> None:
    """Covering the whole population means tolerating the worst observed pair."""
    penalties = [-10.0, -20.0, -30.0, -40.0, -50.0]

    advice = recommended_score_min(
        penalties, percentile=100, read_length=100, slope=-0.6
    )

    assert advice.threshold == pytest.approx(-50.0)
    assert advice.achieved_fraction == 1.0


def test_recommended_score_min_compares_current_setting() -> None:
    """A supplied current intercept is evaluated against the same distribution."""
    penalties = [-10.0, -20.0, -30.0, -40.0, -50.0]

    permissive = recommended_score_min(
        penalties, percentile=95, read_length=100, slope=-0.6, current_intercept=-0.6
    )
    # -0.6 + (-0.6 * 100) = -60.6, more permissive than the worst pair.
    assert permissive.current_threshold == pytest.approx(-60.6)
    assert permissive.current_fraction == 1.0
    assert permissive.current_is_sufficient

    strict = recommended_score_min(
        penalties, percentile=95, read_length=100, slope=-0.1, current_intercept=-5.0
    )
    # -5.0 + (-0.1 * 100) = -15.0, which only the best pair clears.
    assert strict.current_threshold == pytest.approx(-15.0)
    assert strict.current_fraction == pytest.approx(0.2)
    assert not strict.current_is_sufficient


def test_score_min_flag_is_pasteable() -> None:
    """The advice renders as a bowtie2 flag the user can copy."""
    advice = recommended_score_min(
        [-10.0, -20.0], percentile=100, read_length=100, slope=-0.6
    )

    assert advice.score_min_flag.startswith('--score-min L,')
    assert advice.score_min_flag.endswith(',-0.6')


def test_current_is_sufficient_without_a_current_setting() -> None:
    """With no current setting supplied there is nothing to call sufficient."""
    advice = recommended_score_min(
        [-10.0, -20.0], percentile=95, read_length=100, slope=-0.6
    )

    assert advice.current_threshold is None
    assert not advice.current_is_sufficient


def test_recommended_score_min_rejects_empty_input() -> None:
    """A recommendation needs at least one aligned pair."""
    with pytest.raises(ValueError, match='zero alignments'):
        recommended_score_min([], percentile=95, read_length=100, slope=-0.6)


@pytest.mark.parametrize('percentile', [0, -5, 101])
def test_recommended_score_min_rejects_bad_percentile(percentile: float) -> None:
    """The target must be a percentage in (0, 100]."""
    with pytest.raises(ValueError, match='Percentile'):
        recommended_score_min(
            [-10.0], percentile=percentile, read_length=100, slope=-0.6
        )


def test_format_advice_mentions_the_current_setting() -> None:
    """The rendered report includes both the recommendation and the comparison."""
    advice = recommended_score_min(
        [-10.0, -20.0, -30.0],
        percentile=95,
        read_length=100,
        slope=-0.6,
        current_intercept=-0.6,
    )
    text = format_advice(advice)

    assert 'minimum tolerable alignment score' in text
    assert '--score-min L,' in text
    assert 'your current setting' in text
    assert 'sufficient' in text
