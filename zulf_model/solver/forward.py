"""Mixture forward model with variable projection of linear gains.

For nonlinear parameters (couplings, rates, widths, phase delay) the model
renders each component through the observed acquisition operator and solves
the linear response in closed form:

* `gain_model="complex"`: one free complex gain per component (amplitude and
  phase), two real columns each;
* `gain_model="shared_phase"`: one global phase (nonlinear parameter searched
  on a grid inside the projection) and nonnegative real amplitudes per
  component, as expected when all components share one acquisition.

Optional per-band complex polynomial backgrounds (`background_order` >= 0)
absorb smooth residual baselines, such as the low-frequency tail of an
incompletely removed FID baseline. Off by default (-1) because they can bias
decay estimates; any use is recorded in the result.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import brentq, nnls

from ..physics.protocol import SUDDEN_DROP, Protocol
from ..physics.transitions import TransitionCache
from ..render.renderer import ContinuousRenderer, Renderer
from ..timing import Timer
from .observed import ObservedSpectrum
from .parameterization import Parameterization


@dataclass
class Prediction:
    model: np.ndarray                 # complex prediction on observed points
    component_spectra: List[np.ndarray]
    gains: np.ndarray                 # complex gain per component
    background: np.ndarray            # complex constant per band (zeros when disabled)
    residual: np.ndarray              # weighted real residual vector
    score: float


class MixtureForward:
    def __init__(self, parameterization: Parameterization, observed: ObservedSpectrum,
                 protocol: Protocol = SUDDEN_DROP, gain_model: str = "shared_phase", background: int = -1,
                 band_weighting: str = "equal", timer: Optional[Timer] = None, phase_grid: int = 72):
        if gain_model not in ("complex", "shared_phase"):
            raise ValueError("gain_model must be 'complex' or 'shared_phase'.")
        if band_weighting not in ("equal", "none"):
            raise ValueError("band_weighting must be 'equal' or 'none'.")
        self.p = parameterization
        self.obs = observed
        self.protocol = protocol
        self.gain_model = gain_model
        self.background_order = int(background) if not isinstance(background, bool) else (0 if background else -1)
        self.background = self.background_order >= 0
        self.timer = timer or Timer(enabled=False)
        self.phase_grid = phase_grid
        self.cache = TransitionCache(256)
        sel = observed.selected
        self.f = observed.frequencies_hz[sel]
        self.y = observed.values[sel]
        self.band = observed.band_index[sel]
        if observed.acquisition is None:
            self.renderer = ContinuousRenderer()
        else:
            self.renderer = Renderer(observed.acquisition, timer=self.timer)
        scale = np.ones(len(self.y))
        if band_weighting == "equal":
            for b in np.unique(self.band):
                m = self.band == b
                scale[m] = max(float(np.sqrt(np.mean(np.abs(self.y[m]) ** 2))), 1e-30) * math.sqrt(m.sum())
        self.weight = 1.0 / scale
        self.norm = float(np.linalg.norm(self.y * self.weight)) or 1.0

    # -- columns ----------------------------------------------------------------------
    def component_columns(self, values: Dict[str, float]) -> List[np.ndarray]:
        cols = []
        delay = self.p.phase_delay(values)
        for c, system in enumerate(self.p.systems(values)):
            with self.timer.section("solver.transitions"):
                tl = self.cache.get(system, self.protocol)
            edges = self.p.policy.family_edges_hz
            if edges:
                tl = tl.split_families(edges)
            with self.timer.section("solver.render"):
                cols.append(self.renderer.render_pair(tl, self.p.rates(values, c), self.f,
                                                      delay, self.p.sigma(values, c)))
        return cols

    def nuisance_columns(self, values: Dict[str, float]) -> np.ndarray:
        """One complex column per real nuisance amplitude, through the observed operator."""
        terms = self.p.policy.nuisance
        if not terms:
            return np.zeros((len(self.f), 0), complex)
        acq = self.obs.acquisition
        if acq is None:
            raise ValueError("Nuisance terms need a finite-record acquisition.")
        fs = acq.sampling_rate_hz
        t0 = acq.time_origin_s
        cols = []
        with self.timer.section("solver.nuisance"):
            for i, term in enumerate(terms):
                kind = term["kind"]
                if kind == "exponential":
                    rate = math.exp(values[f"n{i}.log_rate"])
                    logz = np.array([-rate / fs])
                    c = np.array([math.exp(-rate * t0)])
                    cols.append(self.renderer.render_modes(logz, c, self.f)[:, 0])
                elif kind == "damped_sinusoid":
                    rate = math.exp(values[f"n{i}.log_rate"])
                    lam = -rate + 2j * math.pi * values[f"n{i}.frequency"]
                    logz = np.array([lam, np.conj(lam)]) / fs
                    base = np.exp(np.array([lam, np.conj(lam)]) * t0) / 2
                    cos_c = base
                    sin_c = base * np.array([-1j, 1j])
                    cols.extend(self.renderer.render_modes(logz, np.column_stack([cos_c, sin_c]), self.f).T)
                elif kind == "template":
                    from ..render.acquisition import evaluate_spectrum, process_record
                    template = np.asarray(term["template"], float)
                    shift = values.get(f"n{i}.shift_s", 0.0)
                    if shift:
                        times = np.arange(len(template)) / fs
                        template = np.interp(times - shift, times, template, left=template[0], right=template[-1])
                    cols.append(evaluate_spectrum(process_record(template, acq), acq, self.f))
        return np.column_stack(cols) if cols else np.zeros((len(self.f), 0), complex)

    def _background_columns(self) -> np.ndarray:
        """Complex polynomial (Legendre-scaled coordinate) per band: columns 1 and i per power."""
        bands = np.unique(self.band)
        order = self.background_order
        out = np.zeros((len(self.f), 2 * len(bands) * (order + 1)), complex)
        col = 0
        for b in bands:
            m = self.band == b
            f = self.f[m]
            u = (f - (f.max() + f.min()) / 2) / max((f.max() - f.min()) / 2, 1e-12)
            for power in range(order + 1):
                basis = np.polynomial.legendre.legval(u, [0] * power + [1])
                out[m, col] = basis
                out[m, col + 1] = 1j * basis
                col += 2
        return out

    @staticmethod
    def _background_real(background: np.ndarray, n_columns: int) -> np.ndarray:
        """Real coefficient vector from the complex packing used in results.

        Background coefficients are stored as complex numbers packing pairs of
        real column coefficients (re, im); a trailing odd column stores its
        real coefficient in the real part.
        """
        real = np.ravel(np.column_stack([background.real, background.imag]))
        return real[:n_columns]

    @staticmethod
    def _real_system(design: np.ndarray, target: np.ndarray, weight: np.ndarray):
        a = np.vstack([(design * weight[:, None]).real, (design * weight[:, None]).imag])
        b = np.r_[(target * weight).real, (target * weight).imag]
        return a, b

    def predict(self, x: Optional[np.ndarray] = None, values: Optional[Dict[str, float]] = None,
                fixed_gains: Optional[np.ndarray] = None, fixed_background: Optional[np.ndarray] = None) -> Prediction:
        values = values if values is not None else self.p.values(x)
        cols = self.component_columns(values)
        n_comp = len(cols)
        bg_cols = self._background_columns() if self.background else np.zeros((len(self.f), 0), complex)
        nuisance = self.nuisance_columns(values)
        if nuisance.shape[1]:
            bg_cols = np.column_stack([bg_cols, nuisance])
        if fixed_gains is not None:
            gains = np.asarray(fixed_gains, complex)
            background = np.zeros((bg_cols.shape[1] + 1) // 2, complex) if fixed_background is None else fixed_background
        elif self.gain_model == "complex":
            design = np.column_stack(cols + [bg_cols]) if n_comp else bg_cols
            a, b = self._real_system(design, self.y, self.weight)
            norms = np.maximum(np.linalg.norm(a, axis=0), 1e-30)
            coef = np.linalg.lstsq(a / norms, b, rcond=1e-10)[0] / norms
            gains = coef[0:2 * n_comp:2] + 1j * coef[1:2 * n_comp:2]
            rest = coef[2 * n_comp:]
            if len(rest) % 2:
                rest = np.r_[rest, 0.0]
            background = rest[0::2] + 1j * rest[1::2]
        else:
            gains, background = self._shared_phase(cols, bg_cols)
        model = sum(col @ np.array([g.real, g.imag]) for col, g in zip(cols, gains)) if n_comp else 0
        component_spectra = [col @ np.array([g.real, g.imag]) for col, g in zip(cols, gains)]
        if len(background):
            model = model + bg_cols @ self._background_real(background, bg_cols.shape[1])
        diff = (model - self.y) * self.weight / self.norm
        residual = np.r_[diff.real, diff.imag]
        return Prediction(np.asarray(model, complex), component_spectra, gains, np.asarray(background, complex),
                          residual, float(residual @ residual))

    def _shared_phase(self, cols, bg_cols):
        """Common phase with nonnegative component amplitudes; background unconstrained.

        A real FID's spectrum is not complex-linear in the gain (mirror term), so
        the response to gain exp(i phi) is cos(phi) col0 + sin(phi) col1. The
        background is projected out once; for each phase only a K x K
        nonnegative least-squares problem remains (Gram-matrix form), so the
        grid-then-golden phase search is cheap.
        """
        c0 = np.column_stack([col[:, 0] for col in cols])
        c1 = np.column_stack([col[:, 1] for col in cols])
        a0, y = self._real_system(c0, self.y, self.weight)
        a1, _ = self._real_system(c1, self.y, self.weight)
        if bg_cols.shape[1]:
            ab, _ = self._real_system(bg_cols, self.y, self.weight)
            q, _ = np.linalg.qr(ab)
            project = lambda m: m - q @ (q.T @ m)
            a0, a1, y_res = project(a0), project(a1), project(y[:, None])[:, 0]
        else:
            y_res = y
        g00, g01, g11 = a0.T @ a0, a0.T @ a1, a1.T @ a1
        r0, r1 = a0.T @ y_res, a1.T @ y_res
        yy = float(y_res @ y_res)
        k = len(cols)

        def solve(phi):
            c, s_ = np.cos(phi), np.sin(phi)
            gram = c * c * g00 + c * s_ * (g01 + g01.T) + s_ * s_ * g11
            rhs = c * r0 + s_ * r1
            scale = np.sqrt(np.maximum(np.diag(gram), 1e-300))
            gram_n = gram / np.outer(scale, scale)
            rhs_n = rhs / scale
            try:
                chol = np.linalg.cholesky(gram_n + 1e-12 * np.eye(k))
            except np.linalg.LinAlgError:
                chol = np.linalg.cholesky(gram_n + 1e-8 * np.eye(k))
            target = np.linalg.solve(chol, rhs_n)
            amp_n, _ = nnls(chol.T, target, maxiter=50 * k)
            amp = amp_n / scale
            value = yy - 2 * amp @ rhs + amp @ gram @ amp
            return amp, value

        def derivative(phi):
            """dV/dphi for the unconstrained amplitudes (no cancellation of O(1) terms at the optimum)."""
            c, s_ = np.cos(phi), np.sin(phi)
            gram = c * c * g00 + c * s_ * (g01 + g01.T) + s_ * s_ * g11
            d_gram = 2 * c * s_ * (g11 - g00) + (c * c - s_ * s_) * (g01 + g01.T)
            rhs = c * r0 + s_ * r1
            d_rhs = -s_ * r0 + c * r1
            amp = np.linalg.solve(gram + 1e-14 * np.trace(gram) * np.eye(k), rhs)
            return -2 * amp @ d_rhs + amp @ d_gram @ amp

        best = None
        for phi in np.linspace(-np.pi, np.pi, self.phase_grid, endpoint=False):
            amp, value = solve(phi)
            if best is None or value < best[2]:
                best = (phi, amp, value)
        step = 2 * np.pi / self.phase_grid
        lo, hi = best[0] - step, best[0] + step
        phi = best[0]
        try:
            d_lo, d_hi = derivative(lo), derivative(hi)
            if d_lo < 0 < d_hi:
                phi = brentq(derivative, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps, maxiter=200)
        except (np.linalg.LinAlgError, ValueError):
            pass
        amp, value = solve(phi)
        if np.any(amp <= 0):
            # Nonnegativity active: fall back to a golden-section search on the NNLS value.
            for _ in range(60):
                m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
                if solve(m1)[1] < solve(m2)[1]:
                    hi = m2
                else:
                    lo = m1
            phi = (lo + hi) / 2
            amp, _ = solve(phi)
        gains = amp * np.exp(1j * phi)
        background = np.zeros(0, complex)
        if bg_cols.shape[1]:
            fitted = sum(np.cos(phi) * a * col[:, 0] + np.sin(phi) * a * col[:, 1] for a, col in zip(amp, cols))
            ab, target = self._real_system(bg_cols, self.y - fitted, self.weight)
            coef = np.linalg.lstsq(ab, target, rcond=1e-12)[0]
            if len(coef) % 2:
                coef = np.r_[coef, 0.0]
            background = coef[0::2] + 1j * coef[1::2]
        return gains, background
