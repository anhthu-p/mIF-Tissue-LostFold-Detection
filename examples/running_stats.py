"""Online first- and second-moment estimation.

Provides :class:`RunningStats`, a single-responsibility accumulator that
tracks the count, mean, and sample variance of a stream of scalars in
constant memory using Welford's algorithm.  It exists here as the unit
under test for the pytest examples in ``test_running_stats.py``.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass


@dataclass
class RunningStats:
    """Constant-memory accumulator for the mean and variance of a stream.

    Attributes
    ----------
    count:
        Number of observations incorporated so far.
    mean:
        Arithmetic mean of the observations; ``0.0`` when ``count == 0``.
    m2:
        Welford's aggregate of squared deviations from the running mean;
        an implementation detail exposed only to support :meth:`merge`.
    """

    count: int = 0
    mean: float = 0.0
    m2: float = 0.0

    def update(self, x: float) -> None:
        """Incorporate one observation via Welford's recurrence.

        The mean is advanced before ``m2`` consumes the post-update
        residual, which keeps the second moment stable over long streams.
        """
        self.count += 1
        delta = x - self.mean
        self.mean += delta / self.count
        self.m2 += delta * (x - self.mean)

    def update_many(self, xs: Iterable[float]) -> None:
        """Incorporate every observation in ``xs`` in iteration order."""
        for x in xs:
            self.update(x)

    @property
    def variance(self) -> float:
        """Unbiased (Bessel-corrected) sample variance.

        Returns ``0.0`` for fewer than two observations, where the sample
        variance is undefined.
        """
        if self.count < 2:
            return 0.0
        return self.m2 / (self.count - 1)

    @property
    def std(self) -> float:
        """Sample standard deviation, the non-negative root of :attr:`variance`."""
        return self.variance ** 0.5

    @classmethod
    def from_iterable(cls, xs: Iterable[float]) -> RunningStats:
        """Construct an accumulator already populated from ``xs``."""
        stats = cls()
        stats.update_many(xs)
        return stats

    def merge(self, other: RunningStats) -> RunningStats:
        """Return a new accumulator equivalent to concatenating both streams.

        Implements Chan et al.'s pairwise variance combination so that
        partials computed independently can be composed without revisiting
        the underlying data.  Neither operand is mutated.
        """
        if other.count == 0:
            return RunningStats(self.count, self.mean, self.m2)
        if self.count == 0:
            return RunningStats(other.count, other.mean, other.m2)
        total = self.count + other.count
        delta = other.mean - self.mean
        mean = self.mean + delta * other.count / total
        m2 = self.m2 + other.m2 + delta ** 2 * self.count * other.count / total
        return RunningStats(total, mean, m2)
