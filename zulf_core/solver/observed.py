"""Observed spectra for refinement and held-out validation.

Two entry points:

* `ObservedSpectrum.from_fid`: an averaged FID processed with an explicit
  acquisition recipe (preferred; enables linewidth continuation and global
  search, which re-process the FID);
* `ObservedSpectrum.from_spectrum`: a spectrum that was already processed
  elsewhere, complex or real (phase-corrected absorption). Supplying the
  processing `record` (sampling rate, record length, crop, SG, mean removal)
  lets the model render exactly the same finite-record lineshapes; without it
  ideal Lorentzian (infinite-record) lines are used. If the spectrum was
  phase-corrected, `phasing` gives the zero-order phase and delay that were
  applied, in the convention of `render.phasing`, and the model is corrected
  identically. Real spectra are compared on the real part only.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np

from ..render.acquisition import Acquisition, evaluate_spectrum, process_record
from ..render.grid import SpectrumGrid
from ..render.phasing import correction_phasor, reference_delay_s


@dataclass
class ObservedSpectrum:
    frequencies_hz: np.ndarray
    values: np.ndarray
    acquisition: Optional[Acquisition]      # None for the continuous (infinite-record) route
    band_index: np.ndarray = None            # band per point; points with -1 are excluded
    label: str = ""
    metadata: dict = field(default_factory=dict)
    fid: Optional[np.ndarray] = None          # raw averaged FID, enables re-processing (continuation)
    real_only: bool = False                  # compare real parts only (phase-corrected absorption spectra)
    phasing: Optional[dict] = None           # {"phase0_rad", "delay_s"} already applied to `values`

    def __post_init__(self):
        self.frequencies_hz = np.asarray(self.frequencies_hz, float)
        self.values = np.asarray(self.values, complex)
        if self.real_only:
            self.values = self.values.real.astype(complex)
        if self.phasing is not None:
            self.phasing = {"phase0_rad": float(self.phasing.get("phase0_rad", 0.0)),
                            "delay_s": float(self.phasing.get("delay_s", 0.0))}
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

    @classmethod
    def from_spectrum(cls, frequencies_hz: Sequence[float], values: Sequence, ranges: Optional[Sequence[Tuple[float, float]]] = None,
                      record=None, phasing: Optional[dict] = None, real_only: Optional[bool] = None,
                      label: str = "") -> "ObservedSpectrum":
        """A spectrum processed elsewhere.

        `values` complex, or real for a phase-corrected absorption spectrum
        (`real_only` defaults to True for real input). `record` is an
        `Acquisition` or its dict (sampling_rate_hz, points, start_sample,
        sg_window, sg_order, remove_mean, ...) describing how the spectrum was
        produced; None selects ideal Lorentzian lines. `ranges` selects bands
        (default: everything above 0 Hz).
        """
        f = np.asarray(frequencies_hz, float)
        v = np.asarray(values)
        if f.shape != v.shape or f.ndim != 1:
            raise ValueError("frequencies_hz and values must be matching 1D arrays.")
        if real_only is None:
            real_only = not np.iscomplexobj(v)
        order = np.argsort(f)
        f, v = f[order], v[order]
        positive = f > 0
        f, v = f[positive], v[positive]
        if len(f) < 3 or np.any(np.diff(f) <= 0):
            raise ValueError("Need at least three distinct positive frequencies.")
        acquisition = Acquisition.from_dict(record) if isinstance(record, dict) else record
        if ranges is None:
            ranges = [(float(f[0]), float(f[-1]))]
        membership = SpectrumGrid(f).band_membership(ranges)
        keep = membership >= 0
        if not np.any(keep):
            raise ValueError("No spectrum points inside the requested ranges.")
        return cls(f[keep], np.asarray(v[keep], complex), acquisition, membership[keep], label,
                   {"ranges": [list(r) for r in ranges], "source": "processed_spectrum"}, None, bool(real_only),
                   phasing)

    def model_correction(self, frequencies_hz: np.ndarray) -> Optional[np.ndarray]:
        """Phasor applied to model spectra so they match an already phase-corrected observation."""
        if self.phasing is None:
            return None
        reference = reference_delay_s(self.acquisition) if self.acquisition is not None else 0.0
        return correction_phasor(frequencies_hz, self.phasing["phase0_rad"], self.phasing["delay_s"] + reference)

    @property
    def reprocessable(self) -> bool:
        return self.fid is not None and self.acquisition is not None

    def with_acquisition(self, acquisition: Acquisition) -> "ObservedSpectrum":
        """Re-process the stored FID with another recipe on the same frequencies."""
        if self.fid is None:
            raise ValueError("No stored FID; cannot re-process.")
        values = evaluate_spectrum(process_record(self.fid, acquisition), acquisition, self.frequencies_hz)
        return ObservedSpectrum(self.frequencies_hz, values, acquisition, self.band_index, self.label,
                                dict(self.metadata), self.fid, self.real_only, self.phasing)

    def restricted(self, ranges: Sequence[Tuple[float, float]]) -> "ObservedSpectrum":
        grid = SpectrumGrid(self.frequencies_hz)
        membership = grid.band_membership(ranges)
        keep = membership >= 0
        return ObservedSpectrum(self.frequencies_hz[keep], self.values[keep], self.acquisition, membership[keep],
                                self.label, dict(self.metadata, ranges=[list(r) for r in ranges]), self.fid,
                                self.real_only, self.phasing)
