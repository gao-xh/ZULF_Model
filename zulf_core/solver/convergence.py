"""Stopping a least-squares fit at a requested coupling precision (D55).

scipy's least_squares only has relative tolerances on the whole vector. `CouplingPrecision` is a callback
(least_squares(..., callback=...)) that stops the fit when, for `patience` consecutive accepted steps, every
coupling changed by less than `precision_hz` and the cost fell by less than `cost_rtol` of itself: the couplings
are settled to the requested precision and nothing else (decay rates, field, phase) is still improving the fit.
precision_hz = 0 never stops (the fit runs to its tolerances or max_nfev, the behaviour before D55).
"""
from __future__ import annotations

import math
from typing import Callable, Optional

import numpy as np

DEFAULT_PRECISION_HZ = 0.01


def decimals_for(precision_hz: float) -> int:
    """Digits after the decimal point to report couplings fitted to `precision_hz` (one more than it resolves)."""
    if not precision_hz or precision_hz <= 0:
        return 3
    return max(1, int(math.ceil(-math.log10(precision_hz))) + 1)


class CouplingPrecision:
    def __init__(self, couplings_of: Callable[[np.ndarray], np.ndarray], precision_hz: float = DEFAULT_PRECISION_HZ,
                 patience: int = 2, cost_rtol: float = 1e-6):
        self.couplings_of, self.precision_hz = couplings_of, float(precision_hz)
        self.patience, self.cost_rtol = int(patience), float(cost_rtol)
        self.previous: Optional[tuple] = None
        self.settled = 0
        self.iterations = 0
        self.stopped = False

    def __call__(self, intermediate_result):
        self.iterations += 1
        if self.precision_hz <= 0:
            return
        values = np.asarray(self.couplings_of(np.asarray(intermediate_result.x)), float)
        cost = float(intermediate_result.cost)
        if self.previous is not None:
            prev_values, prev_cost = self.previous
            step = float(np.max(np.abs(values - prev_values))) if values.size else 0.0
            gain = (prev_cost - cost) / max(abs(cost), 1e-300)
            self.settled = self.settled + 1 if (step < self.precision_hz and gain < self.cost_rtol) else 0
        self.previous = (values, cost)
        if self.settled >= self.patience:
            self.stopped = True
            raise StopIteration


def index_couplings(indices) -> Callable[[np.ndarray], np.ndarray]:
    """couplings_of for a vector whose couplings sit at fixed positions (in Hz)."""
    idx = np.asarray(list(indices), int)
    return lambda x: np.asarray(x)[idx]
