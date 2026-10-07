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


def narrow_excess(frequencies_hz: np.ndarray, values: np.ndarray, band: np.ndarray, baseline_hz: float = 8.0) -> np.ndarray:
    """abs(values) minus its running median (window `baseline_hz`, per band): narrow features only."""
    from scipy.ndimage import median_filter
    f = np.asarray(frequencies_hz, float)
    a = np.abs(np.asarray(values))
    excess = np.zeros_like(a)
    for b in np.unique(band):
        m = np.flatnonzero(band == b)
        if len(m) < 3:
            continue
        spacing = float(np.median(np.diff(f[m])))
        width = max(3, int(round(baseline_hz / max(spacing, 1e-12))) | 1)
        width = min(width, len(m) if len(m) % 2 else len(m) - 1)
        excess[m] = a[m] - median_filter(a[m], size=width, mode="nearest")
    return excess


def signal_regions(frequencies_hz: np.ndarray, values: np.ndarray, band: np.ndarray, threshold: float = 4.0,
                   baseline_hz: float = 8.0):
    """Data-only peak cores and noise level for signal-focused weighting.

    Narrow features are separated from broad background by subtracting a
    running median of abs(values); the noise sigma is a robust (MAD) scale of
    that excess; points whose excess exceeds `threshold` sigma are cores. The
    cores depend on the data only.
    """
    excess = narrow_excess(frequencies_hz, values, band, baseline_hz)
    sigma = max(1.4826 * float(np.median(np.abs(excess - np.median(excess)))), 1e-30)
    return excess > threshold * sigma, sigma


def signal_weights(frequencies_hz: np.ndarray, band: np.ndarray, cores: np.ndarray, outside: float = 0.2,
                   taper_hz: float = 2.0) -> np.ndarray:
    """Relative weights in [outside, 1]: 1 at peak cores, Gaussian fall-off (width `taper_hz`) with the distance
    to the nearest core in the same band, `outside` far away. No hard edges."""
    f = np.asarray(frequencies_hz, float)
    d = np.full(len(f), np.inf)
    for b in np.unique(band):
        m = band == b
        cf = f[m & cores]
        if len(cf):
            idx = np.searchsorted(cf, f[m])
            left = np.abs(f[m] - cf[np.clip(idx - 1, 0, len(cf) - 1)])
            right = np.abs(cf[np.clip(idx, 0, len(cf) - 1)] - f[m])
            d[m] = np.minimum(left, right)
    return outside + (1.0 - outside) * np.exp(-0.5 * (d / max(taper_hz, 1e-12)) ** 2)


def height_factor(frequencies_hz: np.ndarray, values: np.ndarray, band: np.ndarray, cores: np.ndarray,
                  power: float, floor: float = 0.1) -> np.ndarray:
    """Relative weight h^-power for signal weighting, h = height of the nearest peak core.

    Each contiguous run of core points is one peak; its height is the largest abs(values) in it over the band
    maximum (at least `floor`). Every point takes the height of the nearest core in its band, so a small
    peak and its flanks count (h_small)^-power times more than the tallest peak (power 1: equal relative
    errors), while empty gaps follow their neighbouring peak instead of the noise level. Data only."""
    f = np.asarray(frequencies_hz, float)
    a = np.abs(np.asarray(values))
    out = np.ones(len(f))
    if not power:
        return out
    for b in np.unique(band):
        m = np.flatnonzero(band == b)
        c = np.flatnonzero(cores[m])
        if not len(c):
            continue
        top = max(float(a[m].max()), 1e-30)
        run = np.concatenate([[0], np.cumsum(np.diff(c) > 1)])          # run index of every core point
        height = np.array([a[m][c[run == r]].max() for r in range(run[-1] + 1)]) / top
        idx = np.clip(np.searchsorted(f[m][c], f[m]), 1, len(c)) - 1
        nxt = np.clip(idx + 1, 0, len(c) - 1)
        near = np.where(np.abs(f[m] - f[m][c][idx]) <= np.abs(f[m][c][nxt] - f[m]), idx, nxt)
        out[m] = np.maximum(height[run[near]], floor) ** (-power)
    return out


def _mix(pair: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """pair[..., 0] w0 + pair[..., 1] w1 (complex times real; avoids the slow mixed-type matmul path)."""
    return pair[..., 0] * weights[0] + pair[..., 1] * weights[1]


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
                 band_weighting: str = "equal", timer: Optional[Timer] = None, phase_grid: int = 72,
                 signal_threshold: float = 4.0, signal_taper_hz: float = 2.0, signal_outside_weight: float = 0.2,
                 signal_baseline_hz: float = 8.0, signal_extra_hz: Sequence[float] = (),
                 signal_height_power: float = 0.0, signal_height_floor: float = 0.1,
                 amplitude_ratios: Optional[Sequence[float]] = None,
                 amplitude_map: Optional[Sequence[Sequence[float]]] = None):
        if gain_model not in ("complex", "shared_phase"):
            raise ValueError("gain_model must be 'complex' or 'shared_phase'.")
        if band_weighting not in ("equal", "none", "signal"):
            raise ValueError("band_weighting must be 'equal', 'none' or 'signal'.")
        self.p = parameterization
        self.obs = observed
        # Linear amplitude map M (components x free amplitudes): component gains are M @ g. The columns are merged
        # as sum_c M[c, j] col_c before the linear solve, so each column of M costs one free (complex) amplitude.
        # `amplitude_ratios` r is the one-column case M = r[:, None] (fixed isotopologue abundances); a map with one
        # column per parent lets minor isotopologues follow their parent at a fixed abundance factor.
        self.amp_map = None
        n_layout = len(parameterization.layouts)
        if amplitude_ratios is not None and amplitude_map is not None:
            raise ValueError("Give amplitude_ratios or amplitude_map, not both.")
        if amplitude_ratios is not None:
            ratios = np.asarray(amplitude_ratios, float)
            if ratios.shape != (n_layout,) or not np.all(np.isfinite(ratios)) or np.any(ratios < 0) \
                    or not ratios.any():
                raise ValueError("amplitude_ratios needs one nonnegative value per component, not all zero.")
            self.amp_map = ratios[:, None]
        elif amplitude_map is not None:
            m = np.asarray(amplitude_map, float)
            if m.ndim != 2 or m.shape[0] != n_layout or not np.all(np.isfinite(m)) or np.any(m < 0) \
                    or np.any(~m.any(axis=0)):
                raise ValueError("amplitude_map needs shape (components, free amplitudes), nonnegative, "
                                 "no empty column.")
            self.amp_map = m
        self.protocol = protocol
        self.gain_model = gain_model
        self.background_order = int(background) if not isinstance(background, bool) else (0 if background else -1)
        self.background = self.background_order >= 0
        self.timer = timer or Timer(enabled=False)
        self.phase_grid = phase_grid
        self.cache = TransitionCache(256)
        self.exchange_cache = None
        # Optional (lo, hi) in Hz: only transitions inside are rendered and differentiated (local refinement of a
        # window: lines farther away contribute only their tails there; None = every transition)
        self.line_band_hz = None
        # Optional set of free parameter names: jacobian() computes only these columns (the others are zero), e.g.
        # for a local fit that frees a few couplings and rates
        self.jacobian_only = None
        sel = observed.selected
        self.f = observed.frequencies_hz[sel]
        self.y = observed.values[sel]
        self.band = observed.band_index[sel]
        if observed.acquisition is None:
            self.renderer = ContinuousRenderer()
        else:
            # Models are rendered over the shortest record that fixes the retained samples (exact; SG is local).
            self.renderer = Renderer(observed.acquisition.local_record(), timer=self.timer)
        # Processed-spectrum observations: the model is phase-corrected exactly as the data were, and a
        # real (absorption) observation is compared on the real part only.
        self.real_only = bool(getattr(observed, "real_only", False))
        self.correction = observed.model_correction(self.f) if hasattr(observed, "model_correction") else None
        scale = np.ones(len(self.y))
        self.signal_mask = None
        self.signal_cores = None
        self.data_signal_mask = None     # data cores only: the same region for every model
        if band_weighting == "equal":
            for b in np.unique(self.band):
                m = self.band == b
                scale[m] = max(float(np.sqrt(np.mean(np.abs(self.y[m]) ** 2))), 1e-30) * math.sqrt(m.sum())
        elif band_weighting == "signal":
            cores, sigma = signal_regions(self.f, self.y, self.band, signal_threshold, signal_baseline_hz)
            self.data_signal_mask = (signal_weights(self.f, self.band, cores, signal_outside_weight, signal_taper_hz)
                                     >= 0.6) if cores.any() else None
            spacing = float(np.median(np.diff(self.f))) if len(self.f) > 1 else 1.0
            for fx in signal_extra_hz:           # model-predicted line positions (second pass)
                cores |= np.abs(self.f - fx) <= 0.5 * spacing + 1e-9
            self.signal_cores, self.noise_sigma = cores, sigma
            self.signal_settings = dict(threshold=signal_threshold, baseline_hz=signal_baseline_hz,
                                        taper_hz=signal_taper_hz, outside=signal_outside_weight,
                                        height_power=signal_height_power, height_floor=signal_height_floor)
            scale = self._signal_scale()
        self.weight = 1.0 / scale
        self.norm = float(np.linalg.norm(self.y * self.weight)) or 1.0

    def _signal_scale(self) -> np.ndarray:
        """Per-point scale sigma / w of signal weighting for the current cores (see signal_weights)."""
        st, cores, sigma = self.signal_settings, self.signal_cores, self.noise_sigma
        if not cores.any():
            self.signal_mask = None
            return np.full(len(self.y), sigma)
        factor = signal_weights(self.f, self.band, cores, st["outside"], st["taper_hz"])
        self.signal_mask = factor >= 0.6
        if st["height_power"]:
            factor = factor * height_factor(self.f, self.y, self.band, cores, st["height_power"], st["height_floor"])
        return sigma / factor

    def add_signal_cores(self, points: np.ndarray) -> int:
        """Join grid points (bool mask, e.g. model_line_points) to the peak cores of signal weighting and
        recompute the weights and the residual norm: the model-line pass of a fitter that keeps one forward.
        Equivalent to building the forward with signal_extra_hz at these points. Returns the number of new
        core points (0 = nothing changed)."""
        if self.signal_cores is None:
            raise ValueError("add_signal_cores needs band_weighting='signal'.")
        new = np.asarray(points, bool) & ~self.signal_cores
        if new.any():
            self.signal_cores = self.signal_cores | new
            self.weight = 1.0 / self._signal_scale()
            self.norm = float(np.linalg.norm(self.y * self.weight)) or 1.0
        return int(new.sum())

    # -- columns ----------------------------------------------------------------------
    def protocol_for(self, values: Dict[str, float]) -> Protocol:
        """The protocol of this evaluation: the fixed one, with the fitted field when the parameterization has one."""
        field_ut = self.p.field_ut(values) if hasattr(self.p, "field_ut") else None
        protocol = self.protocol if field_ut is None else self.protocol.with_field(field_ut)
        gammas = self.p.gamma_overrides(values) if hasattr(self.p, "gamma_overrides") else {}
        return protocol.with_gamma(gammas) if gammas else protocol

    def transitions(self, values: Dict[str, float], c: int, system) -> "TransitionList":
        """Transition list of component c: static, or with chemical exchange of its exchanging groups
        (physics.exchange; the lines then carry their own decay rates)."""
        ex = self.p.exchange_rates(values, c) if hasattr(self.p, "exchange_rates") else {}
        if not ex:
            tl = self.cache.get(system, self.protocol_for(values))
            return tl.within(*self.line_band_hz) if self.line_band_hz is not None else tl
        if self.exchange_cache is None:
            from ..physics.exchange import ExchangeCache
            self.exchange_cache = ExchangeCache(64)
        spins = {spin: k for g, k in ex.items() for spin in system.groups[g]}
        return self.exchange_cache.get(system, spins, self.protocol_for(values))

    def component_columns(self, values: Dict[str, float]) -> List[np.ndarray]:
        cols = []
        delay = self.p.phase_delay(values)
        for c, system in enumerate(self.p.systems(values)):
            with self.timer.section("solver.transitions"):
                tl = self.transitions(values, c, system)
            edges = self.p.policy.family_edges_hz
            if edges:
                tl = tl.split_families(edges)
            with self.timer.section("solver.render"):
                pair = self.renderer.render_pair(tl, self.p.rates(values, c), self.f, delay, self.p.sigma(values, c))
            if self.correction is not None:
                pair = pair * self.correction[:, None]
            cols.append(pair)
        return cols

    def model_line_points(self, values: Dict[str, float], gains: np.ndarray, threshold: float) -> np.ndarray:
        """Grid points where the model's line envelope (`model_envelope`) exceeds `threshold` noise sigma."""
        if self.signal_cores is None:
            return np.zeros(len(self.f), bool)
        return self.model_envelope(values, gains) > threshold * self.noise_sigma

    def model_envelope(self, values: Dict[str, float], gains: np.ndarray, half_width_hz: float = 20.0) -> np.ndarray:
        """Incoherent line envelope of the model on the fitted grid, read from the transition lists (no peak picking).

            E(f) = sum_k abs(g_c a_k) P_R(f - f_k)

        where P_R is the rendered magnitude profile of a unit line with the same
        decay rate and width in this acquisition (finite record, apodization),
        cut at `half_width_hz`. Line phases are ignored, so E marks every place
        the model puts signal: isolated lines of either sign, many broad
        overlapping lines, and holes where lines cancel in the coherent sum. At
        an isolated line E equals the rendered peak magnitude.
        """
        from ..render.renderer import _per_transition
        from ..physics.transitions import TransitionList
        envelope = np.zeros(len(self.f))
        if not len(self.f):
            return envelope
        spacing = float(np.median(np.diff(self.f))) if len(self.f) > 1 else 1.0
        step = spacing / 4
        offsets = np.arange(-half_width_hz, half_width_hz + step / 2, step)
        delay = self.p.phase_delay(values)
        for c, system in enumerate(self.p.systems(values)):
            tl = self.transitions(values, c, system)
            edges = self.p.policy.family_edges_hz
            if edges:
                tl = tl.split_families(edges)
            if not len(tl) or c >= len(gains):
                continue
            rates = _per_transition(self.p.rates(values, c), tl, "decay rate") + tl.rates_of_lines()
            sigma = self.p.sigma(values, c)
            weight = np.abs(gains[c] * tl.amplitudes)
            for r in np.unique(rates):
                members = np.flatnonzero((rates == r) & (weight > 0))
                if not len(members):
                    continue
                f0 = float(np.median(tl.frequencies_hz[members]))
                pair = self.renderer.render_pair(TransitionList(np.array([f0]), np.array([1.0 + 0j])), float(r),
                                                 f0 + offsets, delay, sigma)
                profile = np.abs(pair[:, 0])
                lo = np.searchsorted(self.f, tl.frequencies_hz[members] - half_width_hz)
                hi = np.searchsorted(self.f, tl.frequencies_hz[members] + half_width_hz, side="right")
                for k, i0, i1 in zip(members, lo, hi):
                    if i1 > i0:
                        envelope[i0:i1] += weight[k] * np.interp(self.f[i0:i1] - tl.frequencies_hz[k], offsets,
                                                                 profile, left=0.0, right=0.0)
        return envelope

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
        out = np.column_stack(cols) if cols else np.zeros((len(self.f), 0), complex)
        if self.correction is not None and out.shape[1]:
            out = out * self.correction[:, None]
        return out

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
        if self.correction is not None:
            # Baseline and ringing leftovers are smooth in the record's own frame, not in the phase-corrected
            # frame, where the first-order correction turns them into a ripple of period 1/(crop + delay).
            out = out * self.correction[:, None]
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

    def mismatch(self, model: np.ndarray, mask: Optional[np.ndarray] = None) -> np.ndarray:
        """Model minus data as compared by the objective (real part only for real observations)."""
        d = np.asarray(model) - self.y
        if mask is not None:
            d = d[mask]
        return d.real if self.real_only else d

    def _real_system(self, design: np.ndarray, target: np.ndarray, weight: np.ndarray):
        wd = design * weight[:, None]
        wt = target * weight
        if self.real_only:
            return wd.real, wt.real
        return np.vstack([wd.real, wd.imag]), np.r_[wt.real, wt.imag]

    def _orthonormal_basis(self, a: np.ndarray) -> np.ndarray:
        """Orthonormal basis of the column space.

        Real-only systems contain exactly null columns (the imaginary background terms), which QR would
        turn into arbitrary directions; there the basis comes from an SVD that drops them. Complex
        systems keep QR.
        """
        if not a.shape[1]:
            return np.zeros((a.shape[0], 0))
        if not self.real_only:
            return np.linalg.qr(a)[0]
        u, sv, _ = np.linalg.svd(a, full_matrices=False)
        keep = sv > max(sv.max(), 1e-300) * 1e-12
        return u[:, keep]

    def predict(self, x: Optional[np.ndarray] = None, values: Optional[Dict[str, float]] = None,
                fixed_gains: Optional[np.ndarray] = None, fixed_background: Optional[np.ndarray] = None) -> Prediction:
        values = values if values is not None else self.p.values(x)
        raw_cols = self.component_columns(values)
        cols = raw_cols
        if self.amp_map is not None and fixed_gains is None:
            cols = [sum(w * col for w, col in zip(self.amp_map[:, j], raw_cols) if w)
                    for j in range(self.amp_map.shape[1])]
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
        solved_gains = np.asarray(gains, complex)
        if self.amp_map is not None and fixed_gains is None:
            gains = self.amp_map @ solved_gains               # component gains from the free amplitudes
        component_spectra = [col @ np.array([g.real, g.imag]) for col, g in zip(raw_cols, gains)]
        if len(background):
            model = model + bg_cols @ self._background_real(background, bg_cols.shape[1])
        diff = (model - self.y) * self.weight / self.norm
        residual = diff.real.copy() if self.real_only else np.r_[diff.real, diff.imag]
        if fixed_gains is None:
            self._last = {"values": dict(values), "cols": cols, "bg_cols": bg_cols, "gains": solved_gains,
                          "background": np.asarray(background, complex), "residual": residual}
        return Prediction(np.asarray(model, complex), component_spectra, gains, np.asarray(background, complex),
                          residual, float(residual @ residual))

    def background_only_residual(self, values: Dict[str, float]) -> float:
        """Norm of what remains after fitting only background and nuisance columns.

        Used as the denominator of the signal-relative residual, so a large
        smooth background cannot make a poor molecular fit look good.
        """
        bg_cols = self._background_columns() if self.background else np.zeros((len(self.f), 0), complex)
        nuisance = self.nuisance_columns(values)
        if nuisance.shape[1]:
            bg_cols = np.column_stack([bg_cols, nuisance])
        if not bg_cols.shape[1]:
            return float(np.linalg.norm(self.y))
        a, b = self._real_system(bg_cols, self.y, np.ones(len(self.y)))
        coef = np.linalg.lstsq(a, b, rcond=1e-12)[0]
        return float(np.linalg.norm(a @ coef - b))

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
            q = self._orthonormal_basis(ab)
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

        def derivative(phi, idx=None):
            """dV/dphi for unconstrained amplitudes of the components `idx` (all by default); the others are 0.

            Analytic, so there is no cancellation of O(1) terms at the optimum.
            """
            sel = np.arange(k) if idx is None else np.asarray(idx)
            c, s_ = np.cos(phi), np.sin(phi)
            sub = np.ix_(sel, sel)
            gram = c * c * g00[sub] + c * s_ * (g01 + g01.T)[sub] + s_ * s_ * g11[sub]
            d_gram = 2 * c * s_ * (g11 - g00)[sub] + (c * c - s_ * s_) * (g01 + g01.T)[sub]
            rhs = c * r0[sel] + s_ * r1[sel]
            d_rhs = -s_ * r0[sel] + c * r1[sel]
            amp = np.linalg.solve(gram + 1e-14 * np.trace(gram) * np.eye(len(sel)), rhs)
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
            # Golden section resolves phi only to about sqrt(eps); polish on the active set with the analytic
            # derivative so that the reduced objective is smooth to working precision.
            active = np.flatnonzero(amp > 0)
            if len(active):
                try:
                    # Golden section on a flat minimum leaves an error of order sqrt(eps) rad; bracket wider.
                    for width in (1e-6, 1e-4, 1e-2):
                        a_lo, a_hi = phi - width, phi + width
                        if derivative(a_lo, active) < 0 < derivative(a_hi, active):
                            break
                    if derivative(a_lo, active) < 0 < derivative(a_hi, active):
                        polished = brentq(lambda v: derivative(v, active), a_lo, a_hi, xtol=1e-15,
                                          rtol=4 * np.finfo(float).eps, maxiter=200)
                        amp_p, _ = solve(polished)
                        if np.array_equal(amp_p > 0, amp > 0):
                            phi, amp = polished, amp_p
                except (np.linalg.LinAlgError, ValueError):
                    pass
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

    # -- Jacobian --------------------------------------------------------------------
    def _realify(self, z: np.ndarray) -> np.ndarray:
        """Weighted complex rows (F, ...) to the real residual space."""
        w = z * self.weight.reshape((-1,) + (1,) * (z.ndim - 1))
        return w.real if self.real_only else np.concatenate([w.real, w.imag], axis=0)

    def _component_derivatives(self, values: Dict[str, float], c: int, names: Sequence[str]) -> Dict[str, np.ndarray]:
        """Analytic pair-column derivatives (F, 2) of component c for couplings, log rates and the phase delay."""
        from ..physics.derivatives import transition_derivatives
        system = self.p.systems(values)[c]
        coupling_names = [n for n in names if self.p.parameters[n].kind == "coupling"]
        pairs = [self.p.parameters[n].detail for n in coupling_names]
        ex = self.p.exchange_rates(values, c) if hasattr(self.p, "exchange_rates") else {}
        if ex:
            return self._exchange_component_derivatives(values, c, names, system, coupling_names, pairs, ex)
        with self.timer.section("solver.transition_derivatives"):
            d = transition_derivatives(system, pairs, self.protocol_for(values))
        if self.line_band_hz is not None:
            from ..physics.derivatives import TransitionDerivatives
            keep = (d.frequencies_hz >= self.line_band_hz[0]) & (d.frequencies_hz <= self.line_band_hz[1])
            d = TransitionDerivatives(d.frequencies_hz[keep], d.amplitudes[keep], d.d_amplitudes[:, keep],
                                      d.t_weights[:, keep], d.pairs)
        edges = self.p.policy.family_edges_hz
        families = np.searchsorted(np.asarray(edges, float), d.frequencies_hz, side="right") if edges else \
            np.zeros(len(d), int)
        rates = self.p.rates(values, c)[families]
        sigma = np.full(len(d), float(self.p.sigma(values, c)))
        delay = self.p.phase_delay(values)
        a = d.amplitudes
        rows_b, rows_c, keys = [], [], []
        for n, k in zip(coupling_names, range(len(pairs))):
            rows_b.append(d.d_amplitudes[k] + 2j * np.pi * delay * d.t_weights[k])
            rows_c.append(2j * np.pi * d.t_weights[k])
            keys.append(n)
        for n in names:
            prm = self.p.parameters[n]
            if prm.kind == "log_rate":
                mask = families == prm.detail[0]
                rows_b.append(np.zeros(len(d), complex))
                rows_c.append(-rates * a * mask)
                keys.append(n)
            elif prm.kind == "phase_delay":
                rows_b.append(2j * np.pi * d.frequencies_hz * a)
                rows_c.append(np.zeros(len(d), complex))
                keys.append(n)
        if not keys:
            return {}
        with self.timer.section("solver.render_derivatives"):
            cols = self.renderer.render_pair_directions(d.frequencies_hz, rates, sigma, np.array(rows_b),
                                                        np.array(rows_c), self.f, delay)
        if self.correction is not None:
            cols = cols * self.correction[:, None, None]
        return {k: cols[:, :, i] for i, k in enumerate(keys)}

    def _exchange_component_derivatives(self, values, c, names, system, coupling_names, pairs, ex):
        """Analytic pair-column derivatives of a component with chemical exchange (physics.exchange)."""
        from ..physics.exchange import exchange_transition_derivatives
        spins = {spin: k for g, k in ex.items() for spin in system.groups[g]}
        ex_names = [n for n in names if self.p.parameters[n].kind == "log_exchange"]
        spin_pairs = [[(p, q) for p in system.groups[a] for q in system.groups[b]] for a, b in pairs]
        ex_spins = [list(system.groups[self.p.parameters[n].detail[0]]) for n in ex_names]
        with self.timer.section("solver.exchange_derivatives"):
            d = exchange_transition_derivatives(system, spins, spin_pairs, ex_spins, self.protocol_for(values))
        tl = d.transitions
        edges = self.p.policy.family_edges_hz
        families = np.searchsorted(np.asarray(edges, float), tl.frequencies_hz, side="right") if edges else \
            np.zeros(len(tl), int)
        comp_rates = self.p.rates(values, c)[families]
        rates = comp_rates + tl.rates_of_lines()
        sigma = np.full(len(tl), float(self.p.sigma(values, c)))
        delay = self.p.phase_delay(values)
        a = tl.amplitudes
        rows_b, rows_c, keys = [], [], []
        for k, n in enumerate(coupling_names + ex_names):
            rows_b.append(d.d_amplitudes[k] + 2j * np.pi * delay * d.phase_terms[k])
            rows_c.append(d.t_terms[k])
            keys.append(n)
        for n in names:
            prm = self.p.parameters[n]
            if prm.kind == "log_rate":
                mask = families == prm.detail[0]
                rows_b.append(np.zeros(len(tl), complex))
                rows_c.append(-comp_rates * a * mask)
                keys.append(n)
            elif prm.kind == "phase_delay":
                rows_b.append(2j * np.pi * tl.frequencies_hz * a)
                rows_c.append(np.zeros(len(tl), complex))
                keys.append(n)
        if not keys or not len(tl):
            return {}
        with self.timer.section("solver.render_derivatives"):
            cols = self.renderer.render_pair_directions(tl.frequencies_hz, rates, sigma, np.array(rows_b),
                                                        np.array(rows_c), self.f, delay)
        if self.correction is not None:
            cols = cols * self.correction[:, None, None]
        return {k: cols[:, :, i] for i, k in enumerate(keys)}

    def _numeric_column_derivatives(self, values: Dict[str, float], names: Sequence[str]):
        """Central differences of the component pair columns and nuisance columns (no closed form)."""
        h = 1e-6 * max(1.0, abs(values[names[0]]))
        out = []
        for sign in (1.0, -1.0):
            v = dict(values)
            for n in names:
                v[n] = values[n] + sign * h
            kinds = {self.p.parameters[n].kind for n in names}
            cols = self.component_columns(v) if kinds - {"nuisance"} else None
            nuis = self.nuisance_columns(v) if "nuisance" in kinds else None
            out.append((cols, nuis))
        pairs = None if out[0][0] is None else [(a - b) / (2 * h) for a, b in zip(out[0][0], out[1][0])]
        nuis = None if out[0][1] is None else (out[0][1] - out[1][1]) / (2 * h)
        return pairs, nuis

    def jacobian(self, x: np.ndarray, full: bool = True) -> np.ndarray:
        """Jacobian of `predict(x).residual` by variable projection.

        The derivative of the model at fixed linear coefficients is analytic for
        couplings (eigen-derivatives, physics.derivatives), decay rates and the
        phase delay; Gaussian widths, field components, gyromagnetic ratios and
        nuisance parameters use central differences of their columns only. With A the directions the linear
        solve absorbs (component amplitudes, the shared phase, background and
        nuisance amplitudes) and r the residual,

            J = P_perp dA c - (A^+)^T dA^T r        (Golub and Pereyra),

        where the first term alone is Kaufman's approximation (`full=False`);
        both agree at zero residual.
        """
        x = np.asarray(x, float)
        state = getattr(self, "_last", None)
        values = self.p.values(x)
        if state is None or state["values"] != values:
            self.predict(x)
            state = self._last
        free = self.p.free_names
        driven = {n: [n] + [f for f, leader in self.p.ties.items() if leader == n] for n in free}
        wanted = [f for f in free if self.jacobian_only is None or f in self.jacobian_only]
        cols, gains, bg_cols = state["cols"], state["gains"], state["bg_cols"]
        n_comp, n_free, n_f = len(cols), len(free), len(self.f)
        n_raw = self.amp_map.shape[0] if self.amp_map is not None else n_comp
        n_nuis = self.nuisance_columns(values).shape[1] if self.p.policy.nuisance else 0
        n_bg = bg_cols.shape[1] - n_nuis
        coef_bg = self._background_real(state["background"], bg_cols.shape[1]) if bg_cols.shape[1] else np.zeros(0)
        # d(pair columns of component c) / d x_i and d(nuisance columns) / d x_i.
        dpair = np.zeros((n_raw, n_free, n_f, 2), complex)
        dnuis = np.zeros((n_free, n_f, n_nuis), complex)
        analytic_kinds = ("coupling", "log_rate", "phase_delay", "log_exchange")
        for c in range(n_raw):
            names = [n for n in self.p.order if self.p.parameters[n].kind in analytic_kinds and
                     self.p.parameters[n].component in (c, -1) and any(n in driven[f] for f in wanted)]
            if not names:
                continue
            for n, col in self._component_derivatives(values, c, names).items():
                for i, f in enumerate(free):
                    if n in driven[f]:
                        dpair[c, i] += col
        for i, f in enumerate(free):
            if self.p.parameters[f].kind not in analytic_kinds and f in wanted:
                pairs, nuis = self._numeric_column_derivatives(values, driven[f])
                if pairs is not None:
                    for c in range(n_raw):
                        dpair[c, i] = pairs[c]
                if nuis is not None:
                    dnuis[i] = nuis
        if self.amp_map is not None:
            dpair = np.tensordot(self.amp_map.T, dpair, axes=1)          # merged columns, like predict
        g_real = np.array([[g.real, g.imag] for g in gains]) if n_comp else np.zeros((0, 2))
        dmodel = sum(_mix(dpair[c], g_real[c]) for c in range(n_comp)).T if n_comp else np.zeros((n_f, n_free), complex)
        if n_nuis:
            dmodel = dmodel + np.einsum("ifn,n->fi", dnuis, coef_bg[n_bg:])
        jac = self._realify(dmodel) / self.norm
        # Absorbed directions and their derivatives.
        absorbed = [bg_cols] if bg_cols.shape[1] else []
        phase_block = None
        d_absorbed = [np.concatenate([np.zeros((n_free, n_f, n_bg), complex), dnuis], axis=2)] if bg_cols.shape[1] \
            else []
        if n_comp:
            if self.gain_model == "complex":
                absorbed.append(np.column_stack(cols))
                d_absorbed.append(np.concatenate([dpair[c] for c in range(n_comp)], axis=2))
            else:
                amp = np.abs(gains)
                phi = float(np.angle(gains[np.argmax(amp)])) if amp.max() > 0 else 0.0
                active = [k for k in range(n_comp) if amp[k] > 0]
                if active:
                    rot = np.array([math.cos(phi), math.sin(phi)])
                    tan = np.array([-math.sin(phi), math.cos(phi)])
                    phase_block = (bg_cols.shape[1], active, amp, rot, tan)
                    absorbed.append(np.column_stack([cols[k] @ rot for k in active]))
                    d_absorbed.append(np.stack([_mix(dpair[k], rot) for k in active], axis=2))
                    absorbed.append(sum(amp[k] * (cols[k] @ tan) for k in active)[:, None])
                    d_absorbed.append(sum(amp[k] * _mix(dpair[k], tan) for k in active)[:, :, None])
        if absorbed:
            a = self._realify(np.column_stack(absorbed)) / self.norm
            if not full:
                u, sv, _ = np.linalg.svd(a, full_matrices=False)
                u = u[:, sv > max(sv.max(initial=0.0), 1e-300) * 1e-10]
                return jac - u @ (u.T @ jac)
            # Exact elimination (implicit function theorem): the eliminated variables l satisfy
            # A^T r = 0, so dl/dx = -H^+ (A^T J_x + dA^T r) with H = A^T A + S, where S holds the
            # residual-weighted second derivatives of the model in l. Only the shared phase is
            # nonlinear: d2 model / d phi2 = -(component model), d2 model / d amp_c d phi = tangent_c.
            r = state["residual"]
            k = a.shape[1]
            s_mat = np.zeros((k, k))
            if phase_block is not None:
                first, active, amp, rot, tan = phase_block
                p_idx = k - 1
                comp_model = sum(amp[c] * (cols[c] @ rot) for c in active)
                s_mat[p_idx, p_idx] = -r @ (self._realify(comp_model) / self.norm)
                for j, c in enumerate(active):
                    v = r @ (self._realify(cols[c] @ tan) / self.norm)
                    s_mat[first + j, p_idx] = s_mat[p_idx, first + j] = v
            da = np.concatenate(d_absorbed, axis=2)                                        # (P, F, K)
            da_real = np.stack([self._realify(da[i]) / self.norm for i in range(n_free)])    # (P, m, K)
            scale = np.linalg.norm(a, axis=0)
            use = scale > max(scale.max(initial=0.0), 1e-300) * 1e-12
            a_s = a[:, use] / scale[use]
            h = a_s.T @ a_s + s_mat[np.ix_(use, use)] / np.outer(scale[use], scale[use])
            rhs = a_s.T @ jac + np.einsum("pmk,m->kp", da_real[:, :, use], r) / scale[use][:, None]
            jac = jac - a_s @ (np.linalg.pinv(h, rcond=1e-12, hermitian=True) @ rhs)
        return jac
