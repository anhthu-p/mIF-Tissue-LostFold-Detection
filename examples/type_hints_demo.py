"""Reference examples for Repo Rule 2 (type hints on function arguments).

Every callable below annotates all of its parameters and its return
type.  ``from __future__ import annotations`` defers annotation
evaluation, which lets us use the ``X | None`` union syntax and forward
references even on older interpreters.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from typing import TypeVar

T = TypeVar("T")


def scale(value: float, factor: float = 2.0) -> float:
    """Return ``value * factor``; shows a plain scalar arg with a default."""
    return value * factor


def normalise(xs: Sequence[float]) -> list[float]:
    """Divide each element by the max.

    Accept the widest reasonable input type (``Sequence``) and return the
    most concrete useful type (``list``): permissive in, precise out.
    """
    top = max(xs)
    return [x / top for x in xs]


def group_by_parity(xs: Iterable[int]) -> dict[str, list[int]]:
    """Partition integers into ``{"even": [...], "odd": [...]}``."""
    buckets: dict[str, list[int]] = {"even": [], "odd": []}
    for x in xs:
        buckets["even" if x % 2 == 0 else "odd"].append(x)
    return buckets


def first_match(xs: Iterable[T], predicate: Callable[[T], bool]) -> T | None:
    """Return the first element satisfying ``predicate``, or ``None``.

    ``T`` ties the element type of ``xs`` to the predicate's argument and
    to the return type; ``Callable[[T], bool]`` annotates the higher-order
    argument.
    """
    for x in xs:
        if predicate(x):
            return x
    return None


def summarise(label: str, *values: float, precision: int = 2) -> str:
    """Format a labelled mean; annotates ``*args`` and a keyword-only arg.

    ``*values`` is typed with the element type, so each variadic argument
    is checked as a ``float``.
    """
    mean = sum(values) / len(values)
    return f"{label}: {mean:.{precision}f}"


@dataclass
class Window:
    """Half-open interval ``[start, end)`` on the real line.

    Dataclass fields are annotated exactly like function arguments; the
    annotations drive ``__init__``'s generated signature.
    """

    start: float
    end: float

    def contains(self, x: float) -> bool:
        """Test membership of ``x`` in the half-open interval."""
        return self.start <= x < self.end


if __name__ == "__main__":
    print(scale(21))
    print(normalise([1.0, 2.0, 4.0]))
    print(group_by_parity(range(6)))
    print(first_match(range(10), lambda n: n > 6))
    print(summarise("latency_ms", 12.0, 15.0, 9.0))
    print(Window(0.0, 1.0).contains(0.5))
