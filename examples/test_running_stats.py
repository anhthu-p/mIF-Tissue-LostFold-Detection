"""Behavioural tests for :mod:`running_stats` (Repo Rules 4 and 5).

These capture the *fundamental* contract of :class:`RunningStats` -- the
properties that must never regress -- rather than line-by-line internals.

Layout conventions demonstrated here:

* One test class per unit of behaviour, so a failing node id such as
  ``TestMerge::test_merge_matches_single_pass`` localises the regression.
* A terse one-line docstring on every test stating the property checked.
* Shared setup lifted into a ``fixture``; input variation expressed with
  ``parametrize`` instead of loops.
* Float comparisons go through ``pytest.approx``.

Run it from the repo root::

    uv run pytest examples/
"""

from __future__ import annotations

import math

import pytest

from running_stats import RunningStats

# Hand-computed reference: n = 8, mean = 5.0, sum of squared deviations
# from the mean = 32, so the Bessel-corrected sample variance is 32 / 7.
KNOWN_SAMPLE = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]


@pytest.fixture
def known() -> RunningStats:
    """Accumulator populated from :data:`KNOWN_SAMPLE`."""
    return RunningStats.from_iterable(KNOWN_SAMPLE)


class TestEmpty:
    """State of a freshly constructed accumulator."""

    def test_count_is_zero(self) -> None:
        """No observations recorded yet."""
        assert RunningStats().count == 0

    def test_mean_defaults_to_zero(self) -> None:
        """Mean is defined as 0.0 in the absence of data."""
        assert RunningStats().mean == 0.0

    def test_variance_is_zero_below_two_samples(self) -> None:
        """Sample variance is undefined for n < 2 and reported as 0.0."""
        assert RunningStats.from_iterable([42.0]).variance == 0.0


class TestUpdate:
    """Incremental incorporation of observations."""

    def test_single_value_sets_mean(self) -> None:
        """One observation makes the mean equal to that value."""
        stats = RunningStats()
        stats.update(3.5)
        assert stats.count == 1
        assert stats.mean == pytest.approx(3.5)

    def test_mean_matches_batch_arithmetic_mean(self, known: RunningStats) -> None:
        """Running mean agrees with sum / n over the same data."""
        assert known.mean == pytest.approx(sum(KNOWN_SAMPLE) / len(KNOWN_SAMPLE))

    def test_update_many_equals_repeated_update(self) -> None:
        """`update_many` is a plain fold of `update` over the iterable."""
        folded = RunningStats.from_iterable(KNOWN_SAMPLE)
        manual = RunningStats()
        for x in KNOWN_SAMPLE:
            manual.update(x)
        assert manual == folded


class TestVariance:
    """Second-moment estimates."""

    def test_sample_variance_is_bessel_corrected(self, known: RunningStats) -> None:
        """Divisor is n - 1, matching the hand-computed 32 / 7."""
        assert known.variance == pytest.approx(32 / 7)

    def test_std_is_sqrt_of_variance(self, known: RunningStats) -> None:
        """std is the non-negative square root of variance."""
        assert known.std == pytest.approx(math.sqrt(32 / 7))

    @pytest.mark.parametrize("constant", [-3.0, 0.0, 10.0])
    def test_constant_stream_has_zero_variance(self, constant: float) -> None:
        """A stream of identical values has no dispersion."""
        stats = RunningStats.from_iterable([constant] * 25)
        assert stats.variance == pytest.approx(0.0, abs=1e-12)


class TestMerge:
    """Composition of independently computed partials."""

    @pytest.mark.parametrize("split", [0, 1, 3, 7, 8])
    def test_merge_matches_single_pass(self, split: int) -> None:
        """Merging two partials equals one pass over the concatenation."""
        left = RunningStats.from_iterable(KNOWN_SAMPLE[:split])
        right = RunningStats.from_iterable(KNOWN_SAMPLE[split:])
        merged = left.merge(right)
        whole = RunningStats.from_iterable(KNOWN_SAMPLE)
        assert merged.count == whole.count
        assert merged.mean == pytest.approx(whole.mean)
        assert merged.variance == pytest.approx(whole.variance)

    def test_merge_does_not_mutate_operands(self, known: RunningStats) -> None:
        """Both operands are left unchanged by merge."""
        snapshot = (known.count, known.mean, known.m2)
        known.merge(RunningStats.from_iterable([1.0, 2.0]))
        assert (known.count, known.mean, known.m2) == snapshot
