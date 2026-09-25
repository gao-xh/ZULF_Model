"""Per-dataset FID diagnostics and candidate processing recipes.

Processing parameters differ between experiments, so they are proposed from
each dataset's own FID rather than fixed in code. Diagnostics report:

* the first-sample anomaly (reported, never silently removed);
* the end of the ADC plateau (saturation) at the start of the record;
* the end of fast ringing (short-window high-pass envelope reaching the
  late-time noise level);
* how well a few real exponentials describe the slow baseline after that;
* a list of candidate recipes (crop start, SG window or explicit exponential
  nuisance terms) to be compared on held-out data. None is chosen here.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import List, Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


@dataclass
class FIDDiagnostics:
    sampling_rate_hz: float
    points: int
    first_point_anomaly: bool
    first_point_value: float
    plateau_end_s: float
    plateau_level: Optional[float]
    ringing_end_s: float
    late_noise_rms: float
    baseline_fits: List[dict] = field(default_factory=list)
    candidate_recipes: List[dict] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _odd(n: int) -> int:
    n = max(3, int(n))
    return n if n % 2 else n + 1


def plateau_end(x: np.ndarray, relative_tolerance: float = 1e-3, min_length: int = 3, search: int = 2000) -> tuple:
    """Index after the initial flat run of near-extreme samples (skipping sample 0)."""
    head = x[1:search]
    if len(head) < min_length:
        return 0, None
    level = float(np.median(head[:min_length]))
    tail = x[-min(len(x) // 4, 4000):]
    noise = float(np.std(np.diff(tail))) / np.sqrt(2) if len(tail) > 10 else 0.0
    tolerance = max(relative_tolerance * abs(level), 3 * noise, 1e-12)
    run = 0
    for value in head:
        if abs(value - level) <= tolerance:
            run += 1
        else:
            break
    if run >= min_length and abs(level) >= 0.5 * np.abs(x[:search]).max():
        return run + 1, float(level)
    return 0, None


def ringing_end(x: np.ndarray, fs: float, start: int, block_s: float = 0.002, window_s: float = 0.02,
                slope_window_s: float = 0.01, rate_threshold_per_s: float = 30.0, sustain: int = 3) -> tuple:
    """End of the fast initial burst, judged by the local envelope decay rate.

    The short-window high-pass envelope is computed in blocks; its local log
    slope is fitted over `slope_window_s`. Instrument ringing decays at
    hundreds per second, molecular zero-field signals at a few per second, so
    the burst ends where the local decay rate stays below `rate_threshold_per_s`.
    """
    window = _odd(window_s * fs)
    if window >= len(x) - start:
        return start / fs, float("nan")
    high = x - savgol_filter(x, window, 2, mode="mirror")
    block = max(4, int(block_s * fs))
    n_blocks = (len(x) - start) // block
    rms = np.sqrt(np.mean(high[start:start + n_blocks * block].reshape(n_blocks, block) ** 2, axis=1))
    late = float(np.median(rms[int(0.75 * n_blocks):])) if n_blocks > 8 else float(np.median(rms))
    log_env = np.log(np.maximum(rms, 1e-12))
    span = max(3, int(round(slope_window_s / block_s)))
    times = np.arange(span) * block_s
    quiet = 0
    for k in range(n_blocks - span):
        slope = np.polyfit(times, log_env[k:k + span], 1)[0]
        near_noise = rms[k] <= 3 * late
        if -slope <= rate_threshold_per_s or near_noise:
            quiet += 1
            if quiet >= sustain:
                return (start + (k - sustain + 1) * block) / fs, late
        else:
            quiet = 0
    return (start + n_blocks * block) / fs, late


def fit_exponentials(t: np.ndarray, y: np.ndarray, count: int) -> dict:
    """Least-squares sum of `count` real exponentials plus a constant (diagnostic only)."""
    rates0 = np.geomspace(0.3, 30.0, count)
    p0 = []
    for r in rates0:
        p0 += [float(y[0]) / count, float(r)]
    p0.append(float(y[-1]))

    def model(tt, *p):
        out = np.full_like(tt, p[-1])
        for k in range(count):
            out = out + p[2 * k] * np.exp(-p[2 * k + 1] * tt)
        return out

    lower = [-np.inf, 1e-3] * count + [-np.inf]
    upper = [np.inf, 1e3] * count + [np.inf]
    try:
        p, _ = curve_fit(model, t, y, p0=p0, bounds=(lower, upper), maxfev=20000)
        residual = y - model(t, *p)
        return {"count": count, "rates_per_s": [float(p[2 * k + 1]) for k in range(count)],
                "amplitudes": [float(p[2 * k]) for k in range(count)], "constant": float(p[-1]),
                "residual_rms": float(residual.std()), "converged": True}
    except (RuntimeError, ValueError) as exc:
        return {"count": count, "converged": False, "error": str(exc)}


def diagnose_fid(fid: np.ndarray, sampling_rate_hz: float, max_exponentials: int = 3) -> FIDDiagnostics:
    x = np.asarray(fid, float)
    fs = float(sampling_rate_hz)
    notes = []
    neighbours = x[1:6]
    spread = max(float(np.median(np.abs(neighbours - np.median(neighbours)))), 1e-3 * abs(float(np.median(neighbours))), 1.0)
    first_anomaly = bool(len(x) > 6 and abs(x[0] - np.median(neighbours)) > 20 * spread)
    if first_anomaly:
        notes.append("First decoded sample differs strongly from its neighbours; it is excluded by any crop > 0.")
    p_end, level = plateau_end(x)
    if level is not None:
        notes.append("Flat near-extreme plateau at the start: likely amplifier/ADC saturation; exclude it.")
    r_end, late = ringing_end(x, fs, max(p_end, 1))
    start = int(math.ceil(max(r_end, p_end / fs) * fs))
    t = np.arange(start, len(x)) / fs
    fits = [fit_exponentials(t, x[start:], n) for n in range(1, max_exponentials + 1)]
    sg_rms = {}
    for seconds in (0.05, 0.2):
        w = _odd(seconds * fs)
        if w < len(x):
            sg_rms[w] = float((x - savgol_filter(x, w, 2, mode="mirror"))[start:].std())
    candidates = []
    for multiplier in (1.0, 2.0, 4.0):
        crop = int(math.ceil(max(r_end, p_end / fs) * multiplier * fs)) + 1
        candidates.append({"label": f"crop x{multiplier:g}, SG off, exponential nuisance",
                           "acquisition": {"start_sample": crop, "sg_window": 0, "remove_mean": False},
                           "nuisance_exponentials": max_exponentials})
        for w, rms in sg_rms.items():
            candidates.append({"label": f"crop x{multiplier:g}, SG {w}",
                               "acquisition": {"start_sample": crop, "sg_window": w, "sg_order": 2,
                                               "remove_mean": True},
                               "diagnostic_residual_rms": rms})
    notes.append("Candidates are hypotheses to compare on held-out acquisitions; none is selected automatically.")
    return FIDDiagnostics(fs, len(x), first_anomaly, float(x[0]), p_end / fs, level, r_end, late, fits,
                          candidates, notes)
