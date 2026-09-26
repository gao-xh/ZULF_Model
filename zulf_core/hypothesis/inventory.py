"""Band inventory: the observed lines of a spectrum, grouped into bands.

Lines are local maxima of abs(values) whose narrow excess (abs minus a
running median) exceeds `threshold` noise sigma (the same robust estimate as
the solver's signal weighting). Lines within `instrument_tol_hz` of a listed
instrument frequency are set aside, and lines closer than `merge_gap_hz`
form one band. `support(f)` gives the evidence for a line at any frequency,
for scoring predicted partner lines that did not pass the threshold.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from ..solver.forward import narrow_excess, signal_regions
from ..solver.observed import ObservedSpectrum


@dataclass
class Line:
    frequency_hz: float
    height: float          # abs(value)
    snr: float             # narrow excess / noise sigma
    band: int = -1


@dataclass
class Band:
    index: int
    lo_hz: float
    hi_hz: float
    lines: List[Line]

    @property
    def peak(self) -> Line:
        return max(self.lines, key=lambda l: l.snr)

    @property
    def center_hz(self) -> float:
        w = np.array([l.height for l in self.lines])
        return float(np.sum(w * [l.frequency_hz for l in self.lines]) / w.sum())

    def contains(self, f: float, margin_hz: float = 0.0) -> bool:
        return self.lo_hz - margin_hz <= f <= self.hi_hz + margin_hz

    def to_dict(self) -> dict:
        return {"index": self.index, "lo_hz": self.lo_hz, "hi_hz": self.hi_hz, "peak_hz": self.peak.frequency_hz,
                "peak_snr": self.peak.snr, "lines": [(l.frequency_hz, l.snr) for l in self.lines]}


@dataclass
class Inventory:
    bands: List[Band]
    noise_sigma: float
    frequencies_hz: np.ndarray = field(repr=False)
    excess: np.ndarray = field(repr=False)
    ranges: List[Tuple[float, float]] = field(default_factory=list)
    instrument_lines: List[Line] = field(default_factory=list)
    instrument_hz: Tuple[float, ...] = ()
    instrument_tol_hz: float = 0.6

    @property
    def lines(self) -> List[Line]:
        return [l for b in self.bands for l in b.lines]

    def observable(self, f: float) -> bool:
        """Inside a fitted range and away from instrument lines."""
        inside = any(lo <= f <= hi for lo, hi in self.ranges) if self.ranges else True
        return inside and all(abs(f - x) > self.instrument_tol_hz for x in self.instrument_hz)

    def support(self, f: float, tol_hz: float = 1.0) -> float:
        """Largest narrow excess within tol of f, in noise sigma (0 outside the data)."""
        m = np.abs(self.frequencies_hz - f) <= tol_hz
        return float(self.excess[m].max() / self.noise_sigma) if m.any() else 0.0

    def band_of(self, f: float, margin_hz: float = 1.0) -> Optional[Band]:
        for b in self.bands:
            if b.contains(f, margin_hz):
                return b
        return None

    def summary(self) -> dict:
        return {"noise_sigma": self.noise_sigma, "bands": [b.to_dict() for b in self.bands],
                "instrument_lines": [(l.frequency_hz, l.snr) for l in self.instrument_lines]}


def band_inventory(observed: ObservedSpectrum, instrument_hz: Sequence[float] = (), instrument_tol_hz: float = 0.6,
                   threshold: float = 4.0, baseline_hz: float = 8.0, merge_gap_hz: float = 4.0,
                   min_sigma_fraction: float = 1e-3) -> Inventory:
    """`min_sigma_fraction`: floor for the noise sigma relative to the largest value, so near noise-free
    (simulated) spectra do not turn every wiggle of a multiplet into a line."""
    sel = observed.selected
    f = observed.frequencies_hz[sel]
    y = observed.values[sel]
    band_index = observed.band_index[sel]
    _, sigma = signal_regions(f, y, band_index, threshold, baseline_hz)
    sigma = max(sigma, min_sigma_fraction * float(np.abs(y).max()))
    excess = narrow_excess(f, y, band_index, baseline_hz)
    a = np.abs(y)
    found, instrument = [], []
    for i in range(1, len(f) - 1):
        if band_index[i - 1] != band_index[i] or band_index[i + 1] != band_index[i]:
            continue
        if a[i] >= a[i - 1] and a[i] > a[i + 1] and excess[i] > threshold * sigma:
            line = Line(float(f[i]), float(a[i]), float(excess[i] / sigma))
            (instrument if any(abs(f[i] - x) <= instrument_tol_hz for x in instrument_hz) else found).append(line)
    found.sort(key=lambda l: l.frequency_hz)
    bands: List[Band] = []
    for line in found:
        if bands and line.frequency_hz - bands[-1].hi_hz <= merge_gap_hz:
            bands[-1].lines.append(line)
            bands[-1].hi_hz = line.frequency_hz
        else:
            bands.append(Band(len(bands), line.frequency_hz, line.frequency_hz, [line]))
        line.band = bands[-1].index
    ranges = [(float(lo), float(hi)) for lo, hi in (observed.metadata or {}).get("ranges", [])]
    if not ranges:
        ranges = [(float(f[band_index == b].min()), float(f[band_index == b].max())) for b in np.unique(band_index)]
    return Inventory(bands, float(sigma), f, excess, ranges, instrument, tuple(instrument_hz), instrument_tol_hz)
