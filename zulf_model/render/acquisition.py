"""Acquisition description and the single experimental preprocessing operator.

`process_record` is the only implementation of baseline subtraction, crop and
mean removal. `evaluate_spectrum` is the only implementation of the finite
record Fourier sum. Simulated and experimental data both pass through them
(the analytic renderer is tested against them).

Every processing step is optional. The default `Acquisition` is pure: no crop,
no SG baseline subtraction, no mean removal, time origin zero. Experimental
recipes switch steps on explicitly in configuration.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.signal import savgol_coeffs, savgol_filter


@dataclass(frozen=True)
class Acquisition:
    sampling_rate_hz: float
    points: int
    start_sample: int = 0
    stop_sample: Optional[int] = None
    sg_window: int = 0
    sg_order: int = 2
    remove_mean: bool = False
    time_origin_s: float = 0.0

    def __post_init__(self):
        if not (math.isfinite(self.sampling_rate_hz) and self.sampling_rate_hz > 0):
            raise ValueError("sampling_rate_hz must be positive.")
        if type(self.points) is not int or self.points < 32:
            raise ValueError("points must be an integer >= 32.")
        stop = self.points if self.stop_sample is None else self.stop_sample
        if type(self.start_sample) is not int or type(stop) is not int or not 0 <= self.start_sample < stop <= self.points:
            raise ValueError("Require 0 <= start_sample < stop_sample <= points.")
        if stop - self.start_sample < 32:
            raise ValueError("Retain at least 32 samples.")
        object.__setattr__(self, "stop_sample", stop)
        if type(self.sg_window) is not int or type(self.sg_order) is not int or self.sg_window < 0 or self.sg_order < 0:
            raise ValueError("SG window and order must be nonnegative integers.")
        if self.sg_window and (self.sg_window % 2 == 0 or self.sg_window <= self.sg_order or self.sg_window > self.points):
            raise ValueError("SG window must be odd, larger than the order and no longer than the record.")
        if not math.isfinite(self.time_origin_s):
            raise ValueError("time_origin_s must be finite.")

    @classmethod
    def pure(cls, sampling_rate_hz: float, points: int) -> "Acquisition":
        """Sampling only: every processing step off."""
        return cls(sampling_rate_hz, points)

    @property
    def is_pure(self) -> bool:
        return (self.start_sample == 0 and self.stop_sample == self.points and not self.sg_window
                and not self.remove_mean and self.time_origin_s == 0.0)

    def without_processing(self) -> "Acquisition":
        return Acquisition.pure(self.sampling_rate_hz, self.points)

    def with_processing(self, **changes) -> "Acquisition":
        data = self.to_dict()
        data.update(changes)
        return Acquisition.from_dict(data)

    @property
    def n(self) -> int:
        """Retained sample count."""
        return self.stop_sample - self.start_sample

    @property
    def native_spacing_hz(self) -> float:
        return self.sampling_rate_hz / self.n

    @property
    def nyquist_hz(self) -> float:
        return self.sampling_rate_hz / 2

    def times(self) -> np.ndarray:
        """Acquisition times of the full record in seconds."""
        return self.time_origin_s + np.arange(self.points) / self.sampling_rate_hz

    def sg_coefficients(self) -> np.ndarray:
        """Dot-product SG coefficients for samples m-h .. m+h (empty when disabled)."""
        if not self.sg_window:
            return np.zeros(0)
        return savgol_coeffs(self.sg_window, self.sg_order, use="dot")

    @classmethod
    def from_times(cls, sampling_rate_hz: float, points: int, start_s: float = 0.0,
                   stop_s: Optional[float] = None, **kwargs) -> "Acquisition":
        start = int(math.ceil(start_s * sampling_rate_hz - 1e-9))
        stop = points if stop_s is None else min(points, int(math.ceil(stop_s * sampling_rate_hz - 1e-9)))
        return cls(sampling_rate_hz, points, start, stop, **kwargs)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Acquisition":
        data = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}
        return cls(**data)

    @classmethod
    def load(cls, path) -> "Acquisition":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def process_record(fid: np.ndarray, acquisition: Acquisition) -> np.ndarray:
    """SG baseline subtraction (full record, mirror edges), crop, optional mean removal.

    Accepts a real array of shape (points,) or (..., points). Returns the
    retained processed samples.
    """
    x = np.asarray(fid, dtype=float)
    if x.shape[-1] != acquisition.points:
        raise ValueError(f"Expected {acquisition.points} samples, got {x.shape[-1]}.")
    if acquisition.sg_window:
        x = x - savgol_filter(x, acquisition.sg_window, acquisition.sg_order, mode="mirror", axis=-1)
    y = x[..., acquisition.start_sample:acquisition.stop_sample].copy()
    if acquisition.remove_mean:
        y -= y.mean(axis=-1, keepdims=True)
    return y


def _zero_fill_factor(frequencies: np.ndarray, acquisition: Acquisition, max_factor: int = 64) -> Optional[int]:
    """Return z if all frequencies lie on the grid k fs / (n z), else None."""
    f = np.asarray(frequencies, float)
    if not len(f):
        return 1
    for z in range(1, max_factor + 1):
        k = f * acquisition.n * z / acquisition.sampling_rate_hz
        if np.allclose(k, np.round(k), atol=1e-6, rtol=0) and np.all(np.round(k) >= 0):
            return z
    return None


def evaluate_spectrum(processed: np.ndarray, acquisition: Acquisition, frequencies_hz: np.ndarray,
                      chunk: int = 4096) -> np.ndarray:
    """X(f) = (1/n) sum_m y[m] exp(-2 pi i f m / fs), m counted from the crop start.

    Uses a zero-filled real FFT when all frequencies lie on such a grid, and an
    exact direct sum otherwise. Supports leading batch dimensions.
    """
    y = np.asarray(processed, float)
    n = acquisition.n
    if y.shape[-1] != n:
        raise ValueError("Processed record length does not match the acquisition.")
    f = np.asarray(frequencies_hz, float)
    z = _zero_fill_factor(f, acquisition)
    if z is not None:
        spectrum = np.fft.rfft(y, n=n * z, axis=-1) / n
        k = np.round(f * n * z / acquisition.sampling_rate_hz).astype(int)
        return spectrum[..., k]
    m = np.arange(n)
    out = np.empty(y.shape[:-1] + (len(f),), dtype=complex)
    for first in range(0, len(f), chunk):
        kernel = np.exp(-2j * np.pi * np.outer(m, f[first:first + chunk]) / acquisition.sampling_rate_hz)
        out[..., first:first + chunk] = (y @ kernel) / n
    return out


def spectrum_from_fid(fid: np.ndarray, acquisition: Acquisition, frequencies_hz: np.ndarray) -> np.ndarray:
    return evaluate_spectrum(process_record(fid, acquisition), acquisition, frequencies_hz)
