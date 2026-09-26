"""Spectral evaluation grids and frequency bands."""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np

from .acquisition import Acquisition


@dataclass(frozen=True)
class SpectrumGrid:
    frequencies_hz: np.ndarray
    zero_fill: int = 1

    def __post_init__(self):
        f = np.asarray(self.frequencies_hz, float)
        if f.ndim != 1 or not len(f) or np.any(f <= 0) or np.any(np.diff(f) <= 0):
            raise ValueError("Grid frequencies must be positive and strictly increasing.")
        object.__setattr__(self, "frequencies_hz", f)

    def __len__(self) -> int:
        return len(self.frequencies_hz)

    @classmethod
    def for_acquisition(cls, acquisition: Acquisition, f_min_hz: float, f_max_hz: float,
                        zero_fill: int = 1) -> "SpectrumGrid":
        """Bins k fs / (n z) inside [f_min, f_max], excluding zero frequency."""
        step = acquisition.sampling_rate_hz / (acquisition.n * zero_fill)
        k0 = max(1, int(math.ceil(f_min_hz / step - 1e-9)))
        k1 = int(math.floor(min(f_max_hz, acquisition.nyquist_hz) / step + 1e-9))
        if k1 < k0:
            raise ValueError("Grid range contains no bins.")
        return cls(np.arange(k0, k1 + 1) * step, zero_fill)

    @classmethod
    def uniform(cls, f_min_hz: float, f_max_hz: float, spacing_hz: float) -> "SpectrumGrid":
        """Multiples of `spacing_hz` inside [f_min, f_max], excluding zero."""
        k0 = max(1, int(math.ceil(f_min_hz / spacing_hz - 1e-9)))
        k1 = int(math.floor(f_max_hz / spacing_hz + 1e-9))
        if k1 < k0:
            raise ValueError("Grid range contains no points.")
        return cls(np.arange(k0, k1 + 1) * spacing_hz)

    @classmethod
    def from_spec(cls, spec, acquisition: Acquisition) -> "SpectrumGrid":
        """Model grid: explicit spacing if configured, else the full record's bins times zero fill.

        `spec` is any object with f_min_hz, f_max_hz, zero_fill and spacing_hz (for example
        zulf_model's GridSpec). Using the full (uncropped) record length keeps the grid fixed when
        crops vary.
        """
        spacing = spec.spacing_hz or acquisition.sampling_rate_hz / (acquisition.points * spec.zero_fill)
        return cls.uniform(spec.f_min_hz, min(spec.f_max_hz, acquisition.nyquist_hz), spacing)

    def band_membership(self, ranges: Sequence[Tuple[float, float]]) -> np.ndarray:
        """Band index per grid point, -1 outside every band. Ranges must be disjoint."""
        membership = np.full(len(self), -1, dtype=int)
        previous = -math.inf
        for i, (lo, hi) in enumerate(ranges):
            if not lo < hi or lo <= previous:
                raise ValueError("Ranges must be increasing and disjoint.")
            previous = hi
            membership[(self.frequencies_hz >= lo) & (self.frequencies_hz <= hi)] = i
        return membership

    def subset(self, mask: np.ndarray) -> "SpectrumGrid":
        return SpectrumGrid(self.frequencies_hz[np.asarray(mask, bool)], self.zero_fill)

    def ranges_indices(self, ranges: Sequence[Tuple[float, float]]) -> List[np.ndarray]:
        membership = self.band_membership(ranges)
        return [np.flatnonzero(membership == i) for i in range(len(ranges))]
