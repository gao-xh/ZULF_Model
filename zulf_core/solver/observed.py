"""Observed complex spectra for refinement and held-out validation."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

from ..render.acquisition import Acquisition, evaluate_spectrum, process_record
from ..render.grid import SpectrumGrid


@dataclass
class ObservedSpectrum:
    frequencies_hz: np.ndarray
    values: np.ndarray
    acquisition: Optional[Acquisition]      # None for the continuous (infinite-record) route
    band_index: np.ndarray = None            # band per point; points with -1 are excluded
    label: str = ""
    metadata: dict = field(default_factory=dict)
    fid: Optional[np.ndarray] = None          # raw averaged FID, enables re-processing (continuation)

    def __post_init__(self):
        self.frequencies_hz = np.asarray(self.frequencies_hz, float)
        self.values = np.asarray(self.values, complex)
        if self.frequencies_hz.shape != self.values.shape:
            raise ValueError("Frequencies and values must match.")
        if self.band_index is None:
            self.band_index = np.zeros(len(self.values), int)
        self.band_index = np.asarray(self.band_index, int)
        if not np.isfinite(self.values).all():
            raise ValueError("Observed values must be finite.")

    @property
    def selected(self) -> np.ndarray:
        return self.band_index >= 0

    @property
    def n_bands(self) -> int:
        return int(self.band_index.max()) + 1 if np.any(self.selected) else 0

    @classmethod
    def from_fid(cls, fid: np.ndarray, acquisition: Acquisition, ranges: Sequence[Tuple[float, float]],
                 zero_fill: int = 1, label: str = "") -> "ObservedSpectrum":
        """Process an averaged FID with the acquisition recipe; keep the given bands."""
        lo = min(r[0] for r in ranges)
        hi = max(r[1] for r in ranges)
        grid = SpectrumGrid.for_acquisition(acquisition, lo, hi, zero_fill)
        membership = grid.band_membership(ranges)
        keep = membership >= 0
        values = evaluate_spectrum(process_record(fid, acquisition), acquisition, grid.frequencies_hz[keep])
        return cls(grid.frequencies_hz[keep], values, acquisition, membership[keep], label,
                   {"ranges": [list(r) for r in ranges], "zero_fill": zero_fill}, np.asarray(fid, float))

    @property
    def reprocessable(self) -> bool:
        return self.fid is not None and self.acquisition is not None

    def with_acquisition(self, acquisition: Acquisition) -> "ObservedSpectrum":
        """Re-process the stored FID with another recipe on the same frequencies."""
        if self.fid is None:
            raise ValueError("No stored FID; cannot re-process.")
        values = evaluate_spectrum(process_record(self.fid, acquisition), acquisition, self.frequencies_hz)
        return ObservedSpectrum(self.frequencies_hz, values, acquisition, self.band_index, self.label,
                                dict(self.metadata), self.fid)

    def restricted(self, ranges: Sequence[Tuple[float, float]]) -> "ObservedSpectrum":
        grid = SpectrumGrid(self.frequencies_hz)
        membership = grid.band_membership(ranges)
        keep = membership >= 0
        return ObservedSpectrum(self.frequencies_hz[keep], self.values[keep], self.acquisition, membership[keep],
                                self.label, dict(self.metadata, ranges=[list(r) for r in ranges]), self.fid)
