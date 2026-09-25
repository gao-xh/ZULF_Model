"""Phase-insensitive global search that proposes starting points for refinement.

Local complex refinement (variable projection in `forward.py`) converges only
inside the basin of the correct line pattern: with lines about 0.3 Hz wide, a
start that is a few hertz off has no gradient pointing to the right answer.
This module searches the whole bounded parameter box with a cheaper objective
that depends on line positions and relative intensities but not on phases:

* data: magnitude spectrum of the processed record after a smooth start ramp
  and end taper (removes the broadband 1/f tail of a record that starts
  abruptly), evaluated on the observed bands;
* model: per component, the magnitude of the complex sum over exact
  transitions of amplitude times the lineshape of the same taper with the
  component decay rate (so antiphase neighbours inside a multiplet cancel as
  in the data); only the phase between components is ignored;
* linear part: nonnegative component gains plus a linear baseline per band,
  solved by NNLS; bands are weighted equally.

The objective ignores relative phases between components, the Savitzky-Golay
high-pass and negative-frequency images, so it is used only to generate starts. Every start is refined with the
exact complex forward model, and ranking never uses this objective.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import differential_evolution, dual_annealing, minimize, nnls

from ..physics.protocol import SUDDEN_DROP, Protocol
from ..physics.transitions import TransitionCache
from ..render.acquisition import process_record
from ..timing import Timer
from .observed import ObservedSpectrum
from .parameterization import Parameterization


@dataclass(frozen=True)
class SearchSettings:
    rise_s: float = 0.3            # cosine start ramp applied to the processed record (search only)
    end_taper_s: float = 1.0       # cosine end taper
    zero_fill: int = 2
    decimate: int = 1              # keep every n-th grid point of the search grid
    popsize: int = 15
    maxiter: int = 100
    seed: int = 0
    solutions: int = 4             # distinct starts returned
    distinct_hz: float = 0.5       # minimum coupling difference between returned starts
    max_seconds: float = 600.0
    polish: bool = True
    polish_maxiter: int = 60
    line_margin_hz: float = 15.0   # transitions outside the bands by more than this are ignored
    rate_grid: int = 32
    kernel_zero_fill: int = 16     # lineshape table density relative to the record length
    include_rates: bool = True
    method: str = "differential_evolution"   # or "dual_annealing" (generalized simulated annealing)
    annealing_maxiter: int = 300

    @classmethod
    def from_dict(cls, data: Optional[dict]) -> "SearchSettings":
        data = dict(data or {})
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


def _cosine_taper(n: int, fs: float, rise_s: float, end_s: float) -> np.ndarray:
    t = np.arange(n) / fs
    w = np.ones(n)
    if rise_s > 0:
        a = t < rise_s
        w[a] = 0.5 * (1 - np.cos(np.pi * t[a] / rise_s))
    if end_s > 0 and n > 1:
        start = t[-1] - end_s
        b = t > start
        w[b] *= 0.5 * (1 + np.cos(np.pi * (t[b] - start) / end_s))
    return w


class PatternObjective:
    """Magnitude line-pattern mismatch for a parameterization and an observed FID."""

    def __init__(self, parameterization: Parameterization, observed: ObservedSpectrum,
                 protocol: Protocol = SUDDEN_DROP, settings: SearchSettings = SearchSettings()):
        if not observed.reprocessable:
            raise ValueError("Global search needs the observed FID and acquisition (ObservedSpectrum.from_fid).")
        self.p = parameterization
        self.protocol = protocol
        self.settings = settings
        self.cache = TransitionCache(512)
        acq = observed.acquisition
        fs = acq.sampling_rate_hz
        y = process_record(observed.fid, acq)
        n = len(y)
        if settings.rise_s + settings.end_taper_s >= n / fs:
            raise ValueError("Search tapers are longer than the processed record.")
        self.window = _cosine_taper(n, fs, settings.rise_s, settings.end_taper_s)
        if acq.apodization_rate_per_s:
            self.window = self.window * np.exp(-acq.apodization_rate_per_s * np.arange(n) / fs)
        m = max(1, int(settings.zero_fill)) * n
        freqs = np.fft.rfftfreq(m, 1 / fs)
        spectrum = np.abs(np.fft.rfft(y * self.window, m)) / fs
        ranges = []
        for b in range(observed.n_bands):
            fb = observed.frequencies_hz[observed.band_index == b]
            ranges.append((float(fb.min()), float(fb.max())))
        self.ranges = ranges
        band = np.full(len(freqs), -1)
        for b, (lo, hi) in enumerate(ranges):
            band[(freqs >= lo) & (freqs <= hi)] = b
        keep = np.flatnonzero(band >= 0)[:: max(1, int(settings.decimate))]
        self.f = freqs[keep]
        self.band = band[keep]
        data = spectrum[keep]
        self.row_weight = np.ones(len(data))
        for b in range(len(ranges)):
            sel = self.band == b
            self.row_weight[sel] = 1.0 / max(float(np.linalg.norm(data[sel])), 1e-30)
        self.d = data * self.row_weight
        self.d_norm = float(np.linalg.norm(self.d)) or 1.0
        cols = []
        for b, (lo, hi) in enumerate(ranges):
            sel = self.band == b
            x = np.where(sel, (self.f - 0.5 * (lo + hi)) / max(hi - lo, 1e-9), 0.0)
            for c in (sel.astype(float), x):
                cols += [c * self.row_weight, -c * self.row_weight]
        self.baseline = np.array(cols).T
        # Complex lineshapes K(df) = sum_m w_m exp(-r t_m) exp(-2 pi i df t_m) / fs for a log grid of rates,
        # tabulated densely with the window-centre linear phase removed so that interpolation is smooth.
        lo_rate, hi_rate = parameterization.policy.rate_bounds_per_s
        self.log_rates = np.linspace(math.log(lo_rate), math.log(hi_rate), max(2, settings.rate_grid))
        t = np.arange(n) / fs
        self.t_center = 0.5 * n / fs
        self.t_start = acq.time_origin_s + acq.start_sample / fs
        mk = max(1, int(settings.kernel_zero_fill)) * n
        offsets = np.fft.fftshift(np.fft.fftfreq(mk, 1 / fs))
        near = np.abs(offsets) <= settings.line_margin_hz
        self.kernel_offsets = offsets[near]
        derotate = np.exp(2j * np.pi * self.kernel_offsets * self.t_center)
        self.kernels = np.array([np.fft.fftshift(np.fft.fft(self.window * np.exp(-math.exp(lr) * t), mk))[near] / fs
                                 * derotate for lr in self.log_rates])
        self.f_lo = min(r[0] for r in ranges) - settings.line_margin_hz
        self.f_hi = max(r[1] for r in ranges) + settings.line_margin_hz
        self.evaluations = 0

    def _kernel(self, rate: float, offsets: np.ndarray) -> np.ndarray:
        lr = float(np.clip(math.log(max(rate, 1e-12)), self.log_rates[0], self.log_rates[-1]))
        k = min(int(np.searchsorted(self.log_rates, lr)), len(self.log_rates) - 1)
        k0 = max(k - 1, 0)
        span = self.log_rates[k] - self.log_rates[k0]
        w = (lr - self.log_rates[k0]) / span if span > 0 else 0.0
        table = self.kernels[k0] if w == 0.0 else (1 - w) * self.kernels[k0] + w * self.kernels[k]
        smooth = (np.interp(offsets, self.kernel_offsets, table.real, left=0.0, right=0.0)
                  + 1j * np.interp(offsets, self.kernel_offsets, table.imag, left=0.0, right=0.0))
        return smooth * np.exp(-2j * np.pi * offsets * self.t_center)

    def columns(self, values: Dict[str, float]) -> List[np.ndarray]:
        edges = np.asarray(self.p.policy.family_edges_hz, float)
        cols = []
        for c, system in enumerate(self.p.systems(values)):
            tl = self.cache.get(system, self.protocol)
            fr = np.asarray(tl.frequencies_hz)
            amp = np.asarray(tl.amplitudes, complex)
            sel = (fr >= self.f_lo) & (fr <= self.f_hi)
            fr = fr[sel]
            amp = amp[sel] * np.exp(2j * np.pi * fr * self.t_start)   # amplitude at the crop start
            rates = self.p.rates(values, c)
            family = np.searchsorted(edges, fr) if len(edges) else np.zeros(len(fr), int)
            col = np.zeros(len(self.f), complex)
            for fam in np.unique(family):
                m = family == fam
                offsets = self.f[:, None] - fr[None, m]
                col += self._kernel(float(rates[fam]), offsets) @ amp[m]
            cols.append(np.abs(col) * self.row_weight)
        return cols

    def fit(self, values: Dict[str, float]):
        self.evaluations += 1
        a = np.column_stack(self.columns(values) + [self.baseline])
        coef, residual = nnls(a, self.d, maxiter=50 * a.shape[1])
        return residual / self.d_norm, coef, a

    def __call__(self, values: Dict[str, float]) -> float:
        return self.fit(values)[0]


@dataclass
class SearchResult:
    names: List[str]                 # searched parameter names (subset of free names)
    points: List[np.ndarray]         # full free vectors, best first
    costs: List[float]
    evaluations: int
    elapsed_s: float
    flags: List[str] = field(default_factory=list)

    def summary(self, parameterization: Parameterization) -> dict:
        return {"searched": self.names, "costs": self.costs, "evaluations": self.evaluations,
                "elapsed_s": self.elapsed_s, "flags": self.flags,
                "starts": [{n: float(v) for n, v in parameterization.values(x).items() if n in self.names}
                           for x in self.points],
                "note": "Phase-insensitive magnitude objective; used only to choose starts."}


class _Stop(Exception):
    pass


def global_search(parameterization: Parameterization, observed: ObservedSpectrum,
                  settings: SearchSettings = SearchSettings(), protocol: Protocol = SUDDEN_DROP,
                  timer: Optional[Timer] = None) -> SearchResult:
    """Differential evolution on the pattern objective; returns distinct starts, best first."""
    timer = timer or Timer(enabled=False)
    objective = PatternObjective(parameterization, observed, protocol, settings)
    free = parameterization.free_names
    kinds = ("coupling", "log_rate") if settings.include_rates else ("coupling",)
    searched = [n for n in free if parameterization.parameters[n].kind in kinds]
    index = [free.index(n) for n in searched]
    base = parameterization.vector()
    lower, upper = parameterization.bounds()
    start = time.perf_counter()
    flags: List[str] = []
    if not searched:
        return SearchResult([], [base], [objective(parameterization.values(base))], objective.evaluations,
                            time.perf_counter() - start, ["nothing_to_search"])

    def full(sub: np.ndarray) -> np.ndarray:
        x = base.copy()
        x[index] = sub
        return x

    def cost(sub: np.ndarray) -> float:
        with timer.section("search.evaluate"):
            return objective(parameterization.values(full(sub)))

    def callback(*_args, **_kwargs):
        return time.perf_counter() - start > settings.max_seconds

    bounds = list(zip(lower[index], upper[index]))
    if settings.method == "differential_evolution":
        de = differential_evolution(cost, bounds, popsize=settings.popsize, maxiter=settings.maxiter,
                                    seed=settings.seed, tol=1e-8, polish=False, init="sobol", callback=callback,
                                    x0=base[index])
        population = getattr(de, "population", np.atleast_2d(de.x))
        energies = getattr(de, "population_energies", np.array([de.fun]))
    elif settings.method == "dual_annealing":
        # Independent annealing runs (different seeds) until enough distinct minima or the budget is used.
        found_x, found_f = [], []
        for run in range(max(1, settings.solutions) * 2):
            if time.perf_counter() - start > settings.max_seconds:
                break
            res = dual_annealing(cost, bounds, maxiter=settings.annealing_maxiter, seed=settings.seed + run,
                                 x0=base[index] if run == 0 else None, no_local_search=True,
                                 callback=lambda x, f, context: time.perf_counter() - start > settings.max_seconds)
            found_x.append(np.asarray(res.x, float))
            found_f.append(float(res.fun))
        population, energies = np.array(found_x), np.array(found_f)
    else:
        raise ValueError("SearchSettings.method must be 'differential_evolution' or 'dual_annealing'.")
    if time.perf_counter() - start > settings.max_seconds:
        flags.append("search_budget_exhausted")
    order = np.argsort(energies)
    couplings = np.array([parameterization.parameters[n].kind == "coupling" for n in searched])
    chosen: List[np.ndarray] = []
    chosen_cost: List[float] = []
    for k in order:
        candidate = np.asarray(population[k], float)
        if any(np.max(np.abs(candidate[couplings] - c[couplings]), initial=0.0) < settings.distinct_hz for c in chosen):
            continue
        chosen.append(candidate)
        chosen_cost.append(float(energies[k]))
        if len(chosen) >= settings.solutions:
            break
    if settings.polish:
        polished = []
        for sub, value in zip(chosen, chosen_cost):
            if time.perf_counter() - start > settings.max_seconds:
                polished.append((sub, value))
                continue
            res = minimize(cost, sub, method="L-BFGS-B", bounds=bounds,
                           options={"maxiter": settings.polish_maxiter})
            polished.append((res.x, float(res.fun)) if res.fun < value else (sub, value))
        polished.sort(key=lambda item: item[1])
        chosen = [p[0] for p in polished]
        chosen_cost = [p[1] for p in polished]
    return SearchResult(searched, [full(s) for s in chosen], chosen_cost, objective.evaluations,
                        time.perf_counter() - start, flags)
