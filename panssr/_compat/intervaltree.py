"""Minimal intervaltree-compatible fallback used when dependency is unavailable."""

from dataclasses import dataclass
from typing import Any, List, Set


@dataclass(frozen=True)
class Interval:
    begin: int
    end: int
    data: Any = None


class IntervalTree:
    def __init__(self):
        self._intervals: List[Interval] = []

    def __setitem__(self, key, value):
        if not isinstance(key, slice):
            raise TypeError("IntervalTree indices must be slices")
        if key.start is None or key.stop is None:
            raise ValueError("Interval boundaries must be defined")
        self._intervals.append(Interval(int(key.start), int(key.stop), value))

    def overlap(self, begin: int, end: int) -> Set[Interval]:
        begin = int(begin)
        end = int(end)
        return {iv for iv in self._intervals if iv.begin < end and begin < iv.end}
