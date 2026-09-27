"""Per-dataset phase: delay from the switching edge, global search over (phase0, delay), local fine-tune.

Convention: `zulf_core.render.phasing.phase_correct(values, f, phase0, delay, acquisition)` multiplies by
exp(-i (phase0 + 2 pi f (delay + crop reference))); the crop reference comes from each dataset's own plan, so a
dataset with another crop start or SG window gets its own, consistent correction.

Criteria are registered (`register_phase_criterion`); each maps (phase0, delay) to a cost from the processed
spectrum alone. Validation so far (ANALYSIS_LOG, D44): on the confirmed J-spectra the model-free criteria fix
the delay only near the switching edge and leave the zero-order phase uncertain by 10-80 deg (lines of both signs
overlap); the edge delay is measured per dataset from the raw FID. The criteria are the part to optimise next.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

from zulf_core.render.acquisition import Acquisition
from zulf_core.render.phasing import phase_correct, reference_delay_s, remove_smooth_background


def signal_regions_hz(values: np.ndarray, frequencies_hz: np.ndarray, threshold: float = 8.0, margin_hz: float = 3.0,
                      exclude_hz: Sequence[float] = (), exclude_width_hz: float = 1.0,
                      background_width_hz: float = 3.0) -> List[Tuple[float, float]]:
    """Intervals holding narrow signal: background-free magnitude above `threshold` robust noise levels, widened
    by `margin_hz`, instrument lines at `exclude_hz` removed. Data only."""
    f = np.asarray(frequencies_hz, float)
    v = remove_smooth_background(values, f, background_width_hz)
    noise = 1.4826 * float(np.median(np.abs(v.real - np.median(v.real)))) or 1e-300
    hit = np.abs(v) > threshold * noise
    for line in exclude_hz:
        hit &= np.abs(f - line) > exclude_width_hz
    if not hit.any():
        return []
    step = float(np.median(np.diff(f))) if len(f) > 1 else 1.0
    grow = max(1, int(round(margin_hz / step)))
    mask = np.convolve(hit.astype(float), np.ones(2 * grow + 1), mode="same") > 0
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    return [(float(f[a]), float(f[b - 1])) for a, b in zip(edges[0::2], edges[1::2])]


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



class PhaseCriterion:
    """Cost of a trial phase (lower is better). `prepare` sees the spectrum once; `grid` evaluates one delay for
    many zero-order phases (vectorised); `cost` one point (fine-tune)."""
    name = "criterion"

    def prepare(self, values, f, acquisition, exclude_hz=()):
        raise NotImplementedError

    def grid(self, delay_s: float, phases_rad: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def cost(self, phase0_rad: float, delay_s: float) -> float:
        return float(self.grid(delay_s, np.array([phase0_rad]))[0])


class DerivativeEntropy(PhaseCriterion):
    """Entropy of |d Re / df| inside the signal regions (ACME, Chen et al. 2002, without the negativity penalty:
    J-spectrum lines have either sign). Invariant under phase0 -> phase0 + pi."""
    name = "entropy"

    def prepare(self, values, f, acquisition, exclude_hz=()):
        regions = signal_regions_hz(values, f, exclude_hz=exclude_hz)
        segs = [np.flatnonzero((f >= lo) & (f <= hi)) for lo, hi in regions]
        segs = [s for s in segs if len(s) >= 4]
        if not segs:
            raise ValueError("No signal regions for the entropy criterion.")
        idx = np.concatenate(segs)
        self.f, self.v = f[idx], np.asarray(values, complex)[idx]
        starts = np.cumsum([0] + [len(s) for s in segs])
        self.keep = np.ones(len(idx) - 1, bool)
        self.keep[starts[1:-1] - 1] = False
        self.ref = reference_delay_s(acquisition) if acquisition is not None else 0.0
        self.regions = regions
        return self

    def grid(self, delay_s, phases_rad):
        base = self.v * np.exp(-2j * np.pi * self.f * (delay_s + self.ref))
        real = (base[None, :] * np.exp(-1j * np.asarray(phases_rad))[:, None]).real
        der = np.abs(np.diff(real, axis=1))[:, self.keep]
        p = der / np.maximum(der.sum(axis=1, keepdims=True), 1e-300)
        return -(p * np.log(p + 1e-300)).sum(axis=1)


class LineCoherence(PhaseCriterion):
    """Doubled-angle coherence of per-line phases (complex Lorentzian fits of resolved lines; sign-agnostic):
    cost = -sum w cos 2(theta_k - phase0 - 2 pi f_k (delay + reference)) / sum w."""
    name = "lines"

    def prepare(self, values, f, acquisition, exclude_hz=()):
        lines = extract_lines(values, f, exclude_hz=exclude_hz)
        if len(lines) < 3:
            raise ValueError("Fewer than three resolved lines for the line criterion.")
        self.fk = np.array([l[0] for l in lines])
        c = np.array([l[1] for l in lines])
        self.theta, self.w = np.angle(c), np.abs(c) ** 2
        self.ref = reference_delay_s(acquisition) if acquisition is not None else 0.0
        return self

    def grid(self, delay_s, phases_rad):
        alpha = self.theta - 2 * np.pi * self.fk * (delay_s + self.ref)
        return -(self.w[None, :] * np.cos(2 * (alpha[None, :] - np.asarray(phases_rad)[:, None]))).sum(axis=1) / \
            self.w.sum()


PHASE_CRITERIA: Dict[str, Callable[[], PhaseCriterion]] = {}


def register_phase_criterion(name: str, factory: Callable[[], PhaseCriterion]) -> None:
    PHASE_CRITERIA[name] = factory


register_phase_criterion("entropy", DerivativeEntropy)
register_phase_criterion("lines", LineCoherence)


def extract_lines(values, frequencies_hz, exclude_hz=(), exclude_width_hz=1.0, threshold=8.0, window_hz=1.2,
                  min_separation_hz=0.35, max_lines=40, max_quality=0.35) -> List[tuple]:
    """Resolved lines as (f_hz, complex height, rate, fit quality): peaks of the background-free magnitude, each
    fitted with complex Lorentzians (neighbours in the window included); poor fits, rates at bounds and widths
    far from the median width (composites) are dropped."""
    from scipy.signal import find_peaks
    f = np.asarray(frequencies_hz, float)
    v = remove_smooth_background(np.asarray(values, complex), f, 3.0)
    mag = np.abs(v)
    noise = 1.4826 * float(np.median(np.abs(v.real - np.median(v.real)))) or 1e-300
    step = float(np.median(np.diff(f))) if len(f) > 1 else 1.0
    idx, _ = find_peaks(mag, height=threshold * noise, distance=max(1, int(round(min_separation_hz / step))))
    keep = [i for i in idx if all(abs(f[i] - x) > exclude_width_hz for x in exclude_hz)]
    keep = sorted(keep, key=lambda i: -mag[i])[:max_lines]
    peaks = np.array(sorted(f[keep])) if keep else np.zeros(0)
    lines = []
    for i in keep:
        m = np.abs(f - f[i]) <= window_hz
        if m.sum() < 6:
            continue
        fit = _fit_line_window(f[m], v[m], [p for p in peaks if abs(p - f[i]) <= window_hz * 0.8])
        if fit is None:
            continue
        j = int(np.argmin(np.abs(fit["f"] - f[i])))
        rate = float(fit["rate"][j])
        if fit["quality"] > max_quality or not (0.31 < rate < 29.0):
            continue
        lines.append((float(fit["f"][j]), complex(fit["c"][j]) / rate, rate, fit["quality"]))
    if len(lines) >= 3:
        med = float(np.median([l[2] for l in lines]))
        lines = [l for l in lines if med / 2 <= l[2] <= med * 2]
    return lines


@dataclass
class PhaseResult:
    phase0_rad: float
    delay_s: float
    criterion: str
    cost: float
    delay_mode: str
    edge_delay_s: Optional[float]
    candidates: List[dict] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def delay_minus_edge_s(self) -> Optional[float]:
        return None if self.edge_delay_s is None else self.delay_s - self.edge_delay_s

    def to_dict(self) -> dict:
        out = asdict(self)
        out["delay_minus_edge_s"] = self.delay_minus_edge_s
        return out


def _local_minima(profile: np.ndarray, k: int) -> List[int]:
    idx = [i for i in range(len(profile)) if (i == 0 or profile[i] <= profile[i - 1]) and
           (i == len(profile) - 1 or profile[i] <= profile[i + 1])]
    return sorted(idx, key=lambda i: profile[i])[:k]


def phase_dataset(values: np.ndarray, frequencies_hz: np.ndarray, acquisition: Acquisition,
                  edge_delay_s: Optional[float] = None, criterion: str = "entropy", delay_mode: str = "edge_prior",
                  delay_span_s: float = 0.0005, delay_step_s: float = 1e-5, phase_step_deg: float = 1.0,
                  top_k: int = 4, exclude_hz: Sequence[float] = (), fine_tune: bool = True) -> PhaseResult:
    """Phase of one processed spectrum: global search, then fine-tune.

    delay_mode: "edge_prior" (default: delay within +-`delay_span_s` of the switching edge measured in the raw FID,
    zero-order phase over its whole range), "global" (delay over +-`delay_span_s` around zero; needs a wide span),
    "edge_fixed" (delay = edge delay, zero-order phase only). Without the edge, delays that differ by about
    1 / (2 x line spacing) are nearly equivalent for any model-free criterion when the bands sit near harmonics
    (methyl lines at J and 2J: the entropy optimum moved 4 ms on a synthetic spectrum); the edge removes that
    ambiguity, which is why it is the default. The global stage evaluates the criterion on a
    (delay, phase0) grid and keeps the `top_k` best local minima over delay; each is fine-tuned by Nelder-Mead in
    (phase0, delay) (delay held in "edge_fixed"); the best is returned with all candidates. The overall sign
    (phase0 vs phase0 + pi) is a convention: the strongest point of the corrected real part is made positive."""
    f = np.asarray(frequencies_hz, float)
    crit = PHASE_CRITERIA[criterion]().prepare(values, f, acquisition, exclude_hz)
    phases = np.radians(np.arange(0.0, 180.0, phase_step_deg))
    notes = []
    if delay_mode == "edge_prior" and edge_delay_s is None:
        notes.append("no switching edge: delay searched around zero")
    if delay_mode == "edge_fixed":
        if edge_delay_s is None:
            raise ValueError("delay_mode 'edge_fixed' needs the edge delay.")
        delays = np.array([edge_delay_s])
    else:
        centre = edge_delay_s if (delay_mode == "edge_prior" and edge_delay_s is not None) else 0.0
        delays = np.arange(centre - delay_span_s, centre + delay_span_s + delay_step_s / 2, delay_step_s)
    table = np.array([crit.grid(d, phases) for d in delays])           # (delays, phases)
    best_phase = phases[np.argmin(table, axis=1)]
    profile = table.min(axis=1)
    starts = _local_minima(profile, top_k) if len(delays) > 1 else [0]
    candidates = []
    for i in starts:
        x0 = np.array([best_phase[i], delays[i]])
        if fine_tune:
            if len(delays) > 1:
                lo, hi = delays[0], delays[-1]
                fun = lambda x: crit.cost(x[0], float(np.clip(x[1], lo, hi)))
                res = minimize(fun, x0, method="Nelder-Mead",
                               options={"xatol": 1e-6, "fatol": 1e-9, "initial_simplex":
                                        [x0, x0 + [np.radians(3.0), 0], x0 + [0, 5 * delay_step_s]]})
                p, d, c = float(res.x[0]), float(np.clip(res.x[1], lo, hi)), float(res.fun)
            else:
                res = minimize(lambda x: crit.cost(x[0], float(delays[0])), x0[:1], method="Nelder-Mead",
                               options={"xatol": 1e-6, "fatol": 1e-9})
                p, d, c = float(res.x[0]), float(delays[0]), float(res.fun)
        else:
            p, d, c = float(x0[0]), float(x0[1]), float(profile[i])
        candidates.append({"phase0_rad": p, "delay_s": d, "cost": c, "grid_delay_s": float(delays[i])})
    candidates.sort(key=lambda c: c["cost"])
    best = candidates[0]
    phase0 = best["phase0_rad"]
    corrected = phase_correct(np.asarray(values, complex), f, phase0, best["delay_s"], acquisition).real
    if corrected[np.argmax(np.abs(corrected))] < 0:
        phase0 += np.pi
    phase0 = float((phase0 + np.pi) % (2 * np.pi) - np.pi)
    if edge_delay_s is not None and abs(best["delay_s"] - edge_delay_s) > 3e-4:
        notes.append(f"delay {best['delay_s'] * 1e3:.3f} ms is {abs(best['delay_s'] - edge_delay_s) * 1e3:.2f} ms from "
                     f"the switching edge ({edge_delay_s * 1e3:.3f} ms)")
    return PhaseResult(phase0, best["delay_s"], criterion, best["cost"], delay_mode, edge_delay_s, candidates, notes)


def calibrated_phase(edge_delay_s: Optional[float], calibration: dict) -> PhaseResult:
    """Instrument calibration plus this dataset's own delay: delay = switching edge + `delay_offset_s`, zero-order
    phase = `phase0_deg` (absolute sign kept: it comes from fits, not from a sign convention). The calibration is
    per instrument/sequence (configs/confirmed_samples.json "phase_calibration"); the edge is measured per dataset.
    On the confirmed samples this is within 14 deg of each sample's own complex-fit phase at 125-250 Hz (leave-one-out,
    ANALYSIS_LOG), closer than any model-free criterion so far."""
    if edge_delay_s is None:
        raise ValueError("calibrated phase needs the switching edge of this dataset.")
    delay = float(edge_delay_s + calibration.get("delay_offset_s", 0.0))
    phase0 = float(np.radians(calibration["phase0_deg"]))
    phase0 = float((phase0 + np.pi) % (2 * np.pi) - np.pi)
    return PhaseResult(phase0, delay, "calibration", 0.0, "edge_calibrated", float(edge_delay_s),
                       notes=[f"calibration: {calibration.get('source', 'unspecified')}"])
