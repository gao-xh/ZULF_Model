"""Built-in lightweight timers.

Usage:

    from zulf_model.timing import TIMER
    with TIMER.section("render"):
        ...
    TIMER.report()        # {"render": {"count": 10, "total_s": ..., "mean_ms": ...}}

`Timer` instances are cheap; pipelines keep their own and expose `report()`.
Timing measures wall-clock time with `time.perf_counter` and is not CPU time.
"""
from __future__ import annotations

import contextlib
import functools
import threading
import time
from collections import defaultdict
from typing import Callable, Dict, Optional


class Timer:
    def __init__(self, enabled: bool = True):
        self.enabled = enabled
        self._total: Dict[str, float] = defaultdict(float)
        self._count: Dict[str, int] = defaultdict(int)
        self._max: Dict[str, float] = defaultdict(float)
        self._lock = threading.Lock()

    @contextlib.contextmanager
    def section(self, name: str):
        if not self.enabled:
            yield
            return
        start = time.perf_counter()
        try:
            yield
        finally:
            self.add(name, time.perf_counter() - start)

    def add(self, name: str, seconds: float) -> None:
        with self._lock:
            self._total[name] += seconds
            self._count[name] += 1
            self._max[name] = max(self._max[name], seconds)

    def timed(self, name: Optional[str] = None) -> Callable:
        def decorator(func):
            label = name or f"{func.__module__}.{func.__qualname__}"

            @functools.wraps(func)
            def wrapper(*args, **kwargs):
                with self.section(label):
                    return func(*args, **kwargs)
            return wrapper
        return decorator

    def report(self) -> Dict[str, dict]:
        with self._lock:
            return {name: {"count": self._count[name], "total_s": self._total[name],
                           "mean_ms": 1000 * self._total[name] / max(1, self._count[name]),
                           "max_ms": 1000 * self._max[name]}
                    for name in sorted(self._total, key=lambda k: -self._total[k])}

    def summary(self) -> str:
        lines = [f"{'section':40s} {'count':>8s} {'total_s':>10s} {'mean_ms':>10s} {'max_ms':>10s}"]
        for name, row in self.report().items():
            lines.append(f"{name:40s} {row['count']:8d} {row['total_s']:10.3f} {row['mean_ms']:10.2f} {row['max_ms']:10.2f}")
        return "\n".join(lines)

    def reset(self) -> None:
        with self._lock:
            self._total.clear()
            self._count.clear()
            self._max.clear()

    def merge(self, other: "Timer") -> None:
        for name, row in other.report().items():
            with self._lock:
                self._total[name] += row["total_s"]
                self._count[name] += row["count"]
                self._max[name] = max(self._max[name], row["max_ms"] / 1000)


TIMER = Timer()
