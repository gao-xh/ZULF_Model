"""Zero- and first-order phase correction of complex spectra.

Convention: a line with complex amplitude a (signal Re(a exp(2 pi i f t)))
appears in the processed spectrum multiplied by exp(i (phase0 + 2 pi f delay)),
where `delay` is the physical delay between the spin evolution origin and the
acquisition time origin plus the Fourier reference shift of the crop
(`reference_delay_s`). Correction multiplies by the inverse phasor, after which
the real part is the absorption spectrum with the physical line signs.

For experimental data the operator supplies `phase0_rad` and the instrument
`delay_s` (the first-order phase expressed as a time); the crop contribution is
added from the acquisition, so the same values apply to any crop.
"""
from __future__ import annotations

import numpy as np

from .acquisition import Acquisition


def reference_delay_s(acquisition: Acquisition) -> float:
    """First-order phase (as a time) introduced by the Fourier reference of the processed record."""
    if acquisition.phase_reference == "crop_start":
        return acquisition.time_origin_s + acquisition.start_sample / acquisition.sampling_rate_hz
    return 0.0


def correction_phasor(frequencies_hz: np.ndarray, phase0_rad: float, delay_s: float) -> np.ndarray:
    """exp(-i (phase0 + 2 pi f delay)); multiply a spectrum by this to remove that phase."""
    f = np.asarray(frequencies_hz, float)
    return np.exp(-1j * (float(phase0_rad) + 2 * np.pi * f * float(delay_s)))


def phase_correct(values: np.ndarray, frequencies_hz: np.ndarray, phase0_rad: float, delay_s: float,
                  acquisition: Acquisition = None) -> np.ndarray:
    """Apply a zero/first-order correction; includes the crop reference when an acquisition is given."""
    total = float(delay_s) + (reference_delay_s(acquisition) if acquisition is not None else 0.0)
    return np.asarray(values, complex) * correction_phasor(frequencies_hz, phase0_rad, total)


def phase1_to_delay_s(phase1_deg: float, spectral_width_hz: float) -> float:
    """Convert a first-order phase given in degrees across a spectral width into a delay in seconds."""
    return float(phase1_deg) / 360.0 / float(spectral_width_hz)


def remove_smooth_background(values: np.ndarray, frequencies_hz: np.ndarray, width_hz: float = 3.0) -> np.ndarray:
    """Subtract a running-median background (real and imaginary parts separately).

    `values` is a processed spectrum in its native (crop) frame, where the
    broadband tail of an abruptly cropped record is smooth; lines much narrower
    than `width_hz` survive. Apply before any first-order phase correction.
    Contiguous frequency segments are filtered separately.
    """
    from scipy.ndimage import median_filter
    f = np.asarray(frequencies_hz, float)
    v = np.asarray(values, complex)
    out = v.copy()
    if len(f) < 3:
        return out
    step = np.median(np.diff(f))
    breaks = np.flatnonzero(np.diff(f) > 1.5 * step) + 1
    for segment in np.split(np.arange(len(f)), breaks):
        size = max(3, int(round(width_hz / step)) | 1)
        if len(segment) <= size:
            continue
        frame = v[segment]
        out[segment] = frame - (median_filter(frame.real, size, mode="nearest")
                                + 1j * median_filter(frame.imag, size, mode="nearest"))
    return out


def estimate_phase(values: np.ndarray, frequencies_hz: np.ndarray, acquisition: Acquisition = None,
                   background_width_hz: float = 3.0, peaks: int = 16, min_separation_hz: float = 0.5,
                   window_hz: float = 0.15, delay_range_s=(-0.01, 0.01), delay_step_s: float = 5e-6) -> dict:
    """Zero/first-order phase for spectra whose lines may be positive or negative.

    Works in the crop frame (the frame of the processed spectrum), where a line
    at f_k has the constant phase phase0 + 2 pi f_k (reference + delay), modulo
    pi because a line may be negative. For each of the strongest peaks the
    complex values are summed over a symmetric window (the dispersive part of a
    Lorentzian cancels) and the centre is refined by parabolic interpolation of
    the magnitude. The doubled-angle weighted circular mean gives the best
    phase0 for each delay on a grid; the delay with the smallest weighted misfit
    wins, and the overall sign makes the weighted peak sum positive. Returns
    phase0_rad and delay_s in the convention of `phase_correct` (the crop
    reference is not included in delay_s).
    """
    from scipy.signal import find_peaks
    f = np.asarray(frequencies_hz, float)
    v = remove_smooth_background(values, f, background_width_hz)
    ref = reference_delay_s(acquisition) if acquisition is not None else 0.0
    magnitude = np.abs(v)
    step = float(np.median(np.diff(f))) if len(f) > 1 else 1.0
    index, props = find_peaks(magnitude, distance=max(1, int(round(min_separation_hz / step))), prominence=0.0)
    if not len(index):
        raise ValueError("No peaks found for phase estimation.")
    index = index[np.argsort(props["prominences"])[::-1][:peaks]]
    half = max(1, int(round(window_hz / step)))
    fk, theta, w = [], [], []
    for i in index:
        lo, hi = max(0, i - half), min(len(f), i + half + 1)
        if i - lo != hi - 1 - i:
            continue
        a, b, c = magnitude[i - 1:i + 2] if 0 < i < len(f) - 1 else (0.0, magnitude[i], 0.0)
        shift = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
        fk.append(f[i] + float(np.clip(shift, -0.5, 0.5)) * step)
        total = v[lo:hi].sum()
        theta.append(np.angle(total))
        w.append(abs(total) ** 2)
    fk, theta, w = np.array(fk), np.array(theta), np.array(w)
    best = None
    for delay in np.arange(delay_range_s[0], delay_range_s[1] + delay_step_s / 2, delay_step_s):
        alpha = theta - 2 * np.pi * fk * (ref + delay)
        resultant = np.sum(w * np.exp(2j * alpha))
        cost = float(np.sum(w) - np.abs(resultant))
        if best is None or cost < best[0]:
            best = (cost, float(delay), 0.5 * float(np.angle(resultant)))
    cost, delay, phase0 = best
    signs = np.cos(theta - phase0 - 2 * np.pi * fk * (ref + delay))
    if np.sum(np.sqrt(w) * np.sign(signs)) < 0:
        phase0 += np.pi
    phase0 = float((phase0 + np.pi) % (2 * np.pi) - np.pi)
    return {"phase0_rad": phase0, "delay_s": delay, "peaks_hz": fk.tolist(),
            "misfit": cost / float(np.sum(w)), "reference_delay_s": ref,
            "note": "Automatic estimate; overall sign follows the majority of peak weight."}
