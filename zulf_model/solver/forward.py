"""Mixture forward model with variable projection of linear gains.

For nonlinear parameters (couplings, rates, widths, phase delay) the model
renders each component through the observed acquisition operator and solves
the linear response in closed form:

* `gain_model="complex"`: one free complex gain per component (amplitude and
  phase), two real columns each;
* `gain_model="shared_phase"`: one global phase (nonlinear parameter searched
  on a grid inside the projection) and nonnegative real amplitudes per
  component, as expected when all components share one acquisition.

Optional per-band complex constants absorb residual baselines (off by default,
because they can bias decay estimates).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import nnls

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
                 protocol: Protocol = SUDDEN_DROP, gain_model: str = "shared_phase", background: bool = False,
                 band_weighting: str = "equal", timer: Optional[Timer] = None, phase_grid: int = 72):
        if gain_model not in ("complex", "shared_phase"):
            raise ValueError("gain_model must be 'complex' or 'shared_phase'.")
        if band_weighting not in ("equal", "none"):
            raise ValueError("band_weighting must be 'equal' or 'none'.")
        self.p = parameterization
        self.obs = observed
        self.protocol = protocol
        self.gain_model = gain_model
        self.background = background
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

    def _background_columns(self) -> np.ndarray:
        bands = np.unique(self.band)
        out = np.zeros((len(self.f), 2 * len(bands)), complex)
        for k, b in enumerate(bands):
            m = self.band == b
            out[m, 2 * k] = 1.0
            out[m, 2 * k + 1] = 1j
        return out

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
        if fixed_gains is not None:
            gains = np.asarray(fixed_gains, complex)
            background = np.zeros(bg_cols.shape[1] // 2, complex) if fixed_background is None else fixed_background
        elif self.gain_model == "complex":
            design = np.column_stack(cols + [bg_cols]) if n_comp else bg_cols
            a, b = self._real_system(design, self.y, self.weight)
            norms = np.maximum(np.linalg.norm(a, axis=0), 1e-30)
            coef = np.linalg.lstsq(a / norms, b, rcond=1e-10)[0] / norms
            gains = coef[0:2 * n_comp:2] + 1j * coef[1:2 * n_comp:2]
            rest = coef[2 * n_comp:]
            background = rest[0::2] + 1j * rest[1::2]
        else:
            gains, background = self._shared_phase(cols, bg_cols)
        model = sum(col @ np.array([g.real, g.imag]) for col, g in zip(cols, gains)) if n_comp else 0
        component_spectra = [col @ np.array([g.real, g.imag]) for col, g in zip(cols, gains)]
        if len(background):
            model = model + bg_cols @ np.ravel(np.column_stack([background.real, background.imag]))
        diff = (model - self.y) * self.weight / self.norm
        residual = np.r_[diff.real, diff.imag]
        return Prediction(np.asarray(model, complex), component_spectra, gains, np.asarray(background, complex),
                          residual, float(residual @ residual))

    def _shared_phase(self, cols, bg_cols):
        """Grid-then-refine search of a common phase with nonnegative amplitudes.

        A real FID's spectrum is not complex-linear in the gain (mirror term), so
        the response to gain exp(i phi) is cos(phi) col0 + sin(phi) col1.
        """
        best = None

        def solve(phi):
            design = np.column_stack([np.cos(phi) * col[:, 0] + np.sin(phi) * col[:, 1] for col in cols])
            a, b = self._real_system(design, self.y, self.weight)
            if bg_cols.shape[1]:
                abg, _ = self._real_system(bg_cols, self.y, self.weight)
                # Background is signed: split into positive and negative parts for NNLS.
                a = np.hstack([a, abg, -abg])
            norms = np.maximum(np.linalg.norm(a, axis=0), 1e-30)
            coef, rnorm = nnls(a / norms, b, maxiter=50 * a.shape[1])
            return coef / norms, rnorm

        grid = np.linspace(-np.pi, np.pi, self.phase_grid, endpoint=False)
        for phi in grid:
            coef, rnorm = solve(phi)
            if best is None or rnorm < best[2]:
                best = (phi, coef, rnorm)
        step = 2 * np.pi / self.phase_grid
        lo, hi = best[0] - step, best[0] + step
        for _ in range(40):  # golden-section refinement
            m1, m2 = lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo)
            if solve(m1)[1] < solve(m2)[1]:
                hi = m2
            else:
                lo = m1
        phi = (lo + hi) / 2
        coef, _ = solve(phi)
        n = len(cols)
        gains = coef[:n] * np.exp(1j * phi)
        background = np.zeros(0, complex)
        if bg_cols.shape[1]:
            k = bg_cols.shape[1]
            signed = coef[n:n + k] - coef[n + k:n + 2 * k]
            background = signed[0::2] + 1j * signed[1::2]
        return gains, background
