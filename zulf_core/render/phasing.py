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

from typing import Optional, Sequence, Tuple

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


def _fit_line_window(f: np.ndarray, v: np.ndarray, centers: Sequence[float], rate_bounds=(0.3, 30.0),
                     shift_hz: float = 0.3) -> Optional[dict]:
    """Complex Lorentzians c_j / (R_j + 2 pi i (f - f_j)) plus a complex constant and slope, fitted to one
    window by variable projection (linear c_j and background, nonlinear f_j and R_j)."""
    from scipy.optimize import least_squares
    centers = np.asarray(centers, float)
    n = len(centers)
    fc = float(np.mean(f))

    def design(x):
        fj, rj = x[:n], np.exp(x[n:])
        cols = [1.0 / (rj[j] + 2j * np.pi * (f - fj[j])) for j in range(n)]
        cols += [np.ones(len(f), complex), (f - fc).astype(complex)]
        return np.column_stack(cols)

    def solve(x):
        a = design(x)
        coef, *_ = np.linalg.lstsq(a, v, rcond=None)
        return a, coef

    def residual(x):
        a, coef = solve(x)
        r = a @ coef - v
        return np.r_[r.real, r.imag]

    x0 = np.r_[centers, np.full(n, np.log(3.0))]
    lo = np.r_[centers - shift_hz, np.full(n, np.log(rate_bounds[0]))]
    hi = np.r_[centers + shift_hz, np.full(n, np.log(rate_bounds[1]))]
    try:
        sol = least_squares(residual, x0, bounds=(lo, hi), max_nfev=200)
    except (ValueError, np.linalg.LinAlgError):
        return None
    a, coef = solve(sol.x)
    fit = a @ coef
    quality = float(np.linalg.norm(fit - v) / max(np.linalg.norm(v - coef[n] - coef[n + 1] * (f - fc)), 1e-300))
    return {"f": sol.x[:n], "rate": np.exp(sol.x[n:]), "c": coef[:n], "quality": quality}


def estimate_phase_lines(values: np.ndarray, frequencies_hz: np.ndarray, acquisition: Acquisition = None,
                         exclude_hz: Sequence[float] = (), exclude_width_hz: float = 1.0, threshold: float = 8.0,
                         window_hz: float = 1.2, min_separation_hz: float = 0.35, max_lines: int = 40,
                         max_quality: float = 0.35, delay_range_s: Tuple[float, float] = (-0.01, 0.01),
                         delay_step_s: float = 2e-6, outlier_deg: float = 25.0) -> dict:
    """Model-free zero/first-order phase of one spectrum from its resolved lines. Coarse; for well-resolved
    spectra only. Validated limits (D43): on a noise-free rendering of the e66a4b08 fit (N-ethylmethylamine, dense
    multiplets with lines of both signs closer than one line width) single-line phases were off by 10-40 deg and
    the delay by 0.5 ms; on the real spectrum the delay was 4.7 ms off. For J-spectra use the phase of a complex
    fit (zulf_hypothesis.report.fit_phasing), which models all lines jointly.

    1. Lines: peaks of the background-free magnitude above `threshold` robust noise levels (instrument lines at
       `exclude_hz` removed), strongest first.
    2. Each line is fitted in a +-`window_hz` window with complex Lorentzians (its neighbours inside the window
       included) and a complex linear background; its phase theta_k = arg(c_k). Fits whose relative residual
       exceeds `max_quality` (overlap, distortion) are dropped.
    3. theta_k = phase0 + 2 pi f_k (delay + reference) modulo pi (lines may be negative). The delay maximising
       the doubled-angle resultant of the weighted lines (weights |c_k|^2) is taken on a grid over
       `delay_range_s`, refined, then lines whose residual phase exceeds `outlier_deg` are removed and the fit
       repeated (robust).
    Each spectrum is phased on its own. The overall sign (phase0 vs phase0 + pi) is a convention: the strongest
    line is made positive. Returns phase0_rad and delay_s in the convention of `phase_correct` (crop reference
    excluded from delay_s), the lines used with their residual phases, the coherence R (1 = all lines agree),
    and `delay_interval_s` (delays with R within 1 % of the maximum)."""
    from scipy.signal import find_peaks
    f = np.asarray(frequencies_hz, float)
    v = remove_smooth_background(np.asarray(values, complex), f, 3.0)
    mag = np.abs(v)
    noise = 1.4826 * float(np.median(np.abs(v.real - np.median(v.real)))) or 1e-300
    step = float(np.median(np.diff(f))) if len(f) > 1 else 1.0
    idx, props = find_peaks(mag, height=threshold * noise, distance=max(1, int(round(min_separation_hz / step))))
    keep = [i for i in idx if all(abs(f[i] - x) > exclude_width_hz for x in exclude_hz)]
    keep = sorted(keep, key=lambda i: -mag[i])[:max_lines]
    if len(keep) < 3:
        raise ValueError("Fewer than three lines for phase estimation.")
    peaks = np.array(sorted(f[keep]))
    lines = []
    for i in keep:
        f0 = f[i]
        m = np.abs(f - f0) <= window_hz
        if m.sum() < 6:
            continue
        near = [p for p in peaks if abs(p - f0) <= window_hz * 0.8]
        fit = _fit_line_window(f[m], v[m], near)
        if fit is None:
            continue
        j = int(np.argmin(np.abs(fit["f"] - f0)))
        rate = float(fit["rate"][j])
        if fit["quality"] > max_quality or not (0.31 < rate < 29.0):     # rejected: poor fit or rate at a bound
            continue
        # c / R is the line height: a broad fit (background or unresolved overlap) has a large c but no height
        lines.append((float(fit["f"][j]), complex(fit["c"][j]) / rate, rate, fit["quality"]))
    if len(lines) >= 3:
        # single transitions share roughly one width; much broader or narrower "lines" are composites
        med = float(np.median([l[2] for l in lines]))
        lines = [l for l in lines if med / 2 <= l[2] <= med * 2]
    if len(lines) < 3:
        raise ValueError("Fewer than three resolved lines for phase estimation.")
    ref = reference_delay_s(acquisition) if acquisition is not None else 0.0
    fk = np.array([l[0] for l in lines])
    theta = np.angle(np.array([l[1] for l in lines]))
    w = np.abs(np.array([l[1] for l in lines])) ** 2
    use = np.ones(len(fk), bool)
    delays = np.arange(delay_range_s[0], delay_range_s[1] + delay_step_s / 2, delay_step_s)

    def coherence(d, sel):
        alpha = theta[sel][None, :] - 2 * np.pi * fk[sel][None, :] * (np.atleast_1d(d)[:, None] + ref)
        return np.abs((w[sel][None, :] * np.exp(2j * alpha)).sum(axis=1)) / w[sel].sum()

    for _ in range(3):
        R = coherence(delays, use)
        k = int(np.argmax(R))
        # parabolic refinement on the grid
        if 0 < k < len(delays) - 1:
            y0, y1, y2 = R[k - 1:k + 2]
            den = y0 - 2 * y1 + y2
            delay = float(delays[k] + (0.5 * (y0 - y2) / den * delay_step_s if den != 0 else 0.0))
        else:
            delay = float(delays[k])
        alpha = theta - 2 * np.pi * fk * (delay + ref)
        phase0 = float(np.angle((w[use] * np.exp(2j * alpha[use])).sum()) / 2)
        resid = (alpha - phase0 + np.pi / 2) % np.pi - np.pi / 2
        new_use = np.abs(np.degrees(resid)) <= outlier_deg
        if new_use.sum() < 3 or np.array_equal(new_use, use):
            break
        use = new_use
    R = coherence(delays, use)
    near = delays[R >= 0.99 * R.max()]
    corrected = np.array([l[1] for l in lines]) * np.exp(-1j * (phase0 + 2 * np.pi * fk * (delay + ref)))
    strongest = int(np.argmax(w * use))
    if corrected[strongest].real < 0:
        phase0 += np.pi
    phase0 = float((phase0 + np.pi) % (2 * np.pi) - np.pi)
    return {"phase0_rad": phase0, "delay_s": delay, "coherence": float(R.max()), "reference_delay_s": ref,
            "delay_interval_s": (float(near.min()), float(near.max())),
            "lines": [{"f_hz": l[0], "amplitude": abs(l[1]), "rate_per_s": l[2], "quality": l[3],
                       "residual_deg": float(np.degrees(r)), "used": bool(u)}
                      for l, r, u in zip(lines, resid, use)],
            "method": "per-line complex Lorentzian phases, doubled-angle delay fit, outliers removed",
            "note": "Overall sign: strongest line positive (convention); phase0 + pi is equally consistent."}
