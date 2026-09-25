"""Exact analytic rendering of transition lists through the acquisition operator.

For a transition list with amplitudes A_k, frequencies f_k and rates R_k the
real FID is

    x(t) = sum_k Re(A_k exp((2 pi i f_k - R_k) t)),  t = time_origin + m / fs.

Writing x[m] = sum_j c_j z_j^m over conjugate pairs of modes, the processed
finite-record spectrum (SG mirror baseline subtraction, crop, mean removal,
normalized DTFT) is evaluated in closed form:

* interior samples: SG acts as a per-mode gain g_j = 1 - sum_u w_u z_j^u;
* mirror-edge samples: exact corrections from the reflected continuation;
* mean removal: subtract the retained mean times the DTFT of a constant.

`tests/test_render.py` checks this against time-domain synthesis followed by
`process_record` and `evaluate_spectrum`.

Line broadening: each transition decays exponentially with its rate R
(Lorentzian, homogeneous) and optionally carries a Gaussian frequency
distribution of standard deviation sigma (inhomogeneous), giving the envelope
exp(-R t - (2 pi sigma t)^2 / 2), i.e. a Voigt line. Rates and sigmas are
scalars or per-family tables. Finite-record rendering with sigma > 0 uses the
exact sampled path; the continuous renderer uses the Faddeeva function.
"""
from __future__ import annotations

from typing import Sequence, Union

import numpy as np
from scipy.special import wofz

from ..physics.transitions import TransitionList
from ..timing import Timer
from .acquisition import Acquisition, _fft_length, evaluate_spectrum, process_record
from .nufft import nufft_type1

Rates = Union[float, Sequence[float], np.ndarray]


def _per_transition(values: Rates, transitions: TransitionList, name: str, positive: bool = False) -> np.ndarray:
    """Expand a scalar or per-family table to one value per transition."""
    if np.isscalar(values):
        out = np.full(len(transitions), float(values))
    else:
        out = np.asarray(values, float)[transitions.families]
    if not np.isfinite(out).all() or np.any(out < 0) or (positive and np.any(out <= 0)):
        raise ValueError(f"{name} must be finite and {'positive' if positive else 'nonnegative'}.")
    return out


def gaussian_envelope(times_s: np.ndarray, sigma_hz: float) -> np.ndarray:
    """exp(-(2 pi sigma t)^2 / 2): FID envelope of a Gaussian frequency distribution with std sigma."""
    return np.exp(-0.5 * (2 * np.pi * sigma_hz * np.asarray(times_s)) ** 2)


def _geometric(L: np.ndarray, n: int) -> np.ndarray:
    """sum_{m=0}^{n-1} exp(L m), stable near L = 0."""
    small = np.abs(L) < 1e-12
    safe = np.where(small, 1.0, L)
    with np.errstate(invalid="ignore", over="ignore"):
        value = np.expm1(n * safe) / np.expm1(safe)
    return np.where(small, n + (n * (n - 1) / 2) * L, value)


class Renderer:
    """Analytic renderer bound to one acquisition. NumPy backend, float64.

    When a mode decays by more than `exp(stability_limit)` across half an SG
    window, the mirror-edge algebra cancels huge terms; such calls use the exact
    sampled path (time-domain synthesis followed by the acquisition operator).
    """

    BACKENDS = ("auto", "analytic", "time")

    def __init__(self, acquisition: Acquisition, chunk_frequencies: int = 2048,
                 stability_limit: float = 8.0, backend: str = "auto", timer: Timer = None):
        if backend not in self.BACKENDS:
            raise ValueError(f"backend must be one of {self.BACKENDS}.")
        self.acq = acquisition
        self.chunk = chunk_frequencies
        self.stability_limit = stability_limit
        self.backend = backend
        self.timer = timer or Timer(enabled=False)
        self.sampled_calls = 0
        self.analytic_calls = 0
        self.coef = acquisition.sg_coefficients()
        self.half = acquisition.sg_window // 2
        fs, n, start, stop, N, h = (acquisition.sampling_rate_hz, acquisition.n, acquisition.start_sample,
                                    acquisition.stop_sample, acquisition.points, self.half)
        self.offsets = np.arange(-h, h + 1)
        # Edge samples inside the retained range, and their SG correction structure.
        self.start_edge = np.arange(start, min(stop, h)) if h else np.zeros(0, int)
        self.end_edge = np.arange(max(start, N - h), stop) if h else np.zeros(0, int)

    def _time_is_cheaper(self, n_frequencies: int, n_transitions: int, f: np.ndarray) -> bool:
        """Rough cost model: analytic ~ F x 2T; time path ~ N log N plus processing."""
        if _fft_length(f, self.acq) is None:
            return False
        n = self.acq.points
        analytic = n_frequencies * 2 * n_transitions + 2 * n_transitions * self.half * 4
        time_path = 8 * n * max(1.0, np.log2(n))
        return time_path < analytic

    # -- mode bookkeeping ----------------------------------------------------
    def _modes(self, transitions: TransitionList, rates: Rates, gain: complex, phase_delay_s: float):
        f = transitions.frequencies_hz
        r = _per_transition(rates, transitions, "Decay rates")
        lam = -r + 2j * np.pi * f
        a = transitions.amplitudes * gain * np.exp(2j * np.pi * f * phase_delay_s)
        t0 = self.acq.time_origin_s
        c = a / 2 * np.exp(lam * t0)
        lam_all = np.concatenate([lam, lam.conj()])
        c_all = np.concatenate([c, c.conj()])
        return lam_all / self.acq.sampling_rate_hz, c_all  # log z per sample, coefficients

    # -- core -----------------------------------------------------------------
    def _edge_deviation(self, logz: np.ndarray, c: np.ndarray) -> tuple:
        """Deviation of SG-processed edge samples from the interior gain model."""
        h, N = self.half, self.acq.points
        start_dev = np.zeros(len(self.start_edge), complex)
        end_dev = np.zeros(len(self.end_edge), complex)
        if not h:
            return start_dev, end_dev
        v = np.arange(1, h + 1)
        if len(self.start_edge):
            # G(v) = x_v - x_{-v}
            G = (np.exp(np.outer(v, logz)) - np.exp(-np.outer(v, logz))) @ c
            for idx, m in enumerate(self.start_edge):
                vv = np.arange(1, h - m + 1)
                k = h - m - vv
                start_dev[idx] = -np.dot(self.coef[k], G[vv - 1])
        if len(self.end_edge):
            base = np.exp(logz * (N - 1))
            H = (np.exp(-np.outer(v, logz)) - np.exp(np.outer(v, logz))) @ (c * base)
            for idx, m in enumerate(self.end_edge):
                ww = np.arange(1, m + h - (N - 1) + 1)
                k = h + (N - 1) + ww - m
                end_dev[idx] = -np.dot(self.coef[k], H[ww - 1])
        return start_dev, end_dev

    def _dtft_modes(self, logz: np.ndarray, weights: np.ndarray, frequencies: np.ndarray,
                    apodize: bool = True) -> np.ndarray:
        """(1/n) sum_{m=start}^{stop-1} sum_j weights_j z_j^m w_m exp(-i w (m - start)), w_m the apodization."""
        acq = self.acq
        n = acq.n
        omega = 2 * np.pi * np.asarray(frequencies, float) / acq.sampling_rate_hz
        scaled = weights * np.exp(logz * acq.start_sample)
        damp = acq.apodization_rate_per_s / acq.sampling_rate_hz if apodize else 0.0
        out = np.empty(len(omega), complex)
        for first in range(0, len(omega), self.chunk):
            L = logz[None, :] - damp - 1j * omega[first:first + self.chunk, None]
            out[first:first + self.chunk] = _geometric(L, n) @ scaled
        return out / n

    def _kernel(self, frequencies: np.ndarray, rel: np.ndarray) -> np.ndarray:
        """Cached exp(-2 pi i f rel / fs) for edge-sample DTFTs."""
        key = (frequencies.tobytes(), rel.tobytes())
        cache = getattr(self, "_kernel_cache", None)
        if cache is None:
            cache = self._kernel_cache = {}
        if key not in cache:
            if len(cache) > 8:
                cache.clear()
            cache[key] = np.exp(-2j * np.pi * np.outer(frequencies, rel) / self.acq.sampling_rate_hz)
        return cache[key]

    def _dtft_samples(self, sample_index: np.ndarray, values: np.ndarray, frequencies: np.ndarray,
                      apodize: bool = True) -> np.ndarray:
        if not len(sample_index):
            return np.zeros(len(frequencies) if values.ndim == 1 else (len(frequencies),) + values.shape[1:], complex)
        rel = sample_index - self.acq.start_sample
        if apodize and self.acq.apodization_rate_per_s:
            window = np.exp(-self.acq.apodization_rate_per_s * rel / self.acq.sampling_rate_hz)
            values = values * (window if values.ndim == 1 else window[:, None])
        return self._kernel(np.asarray(frequencies, float), rel) @ values / self.acq.n

    def _constant_dtft(self, frequencies: np.ndarray) -> np.ndarray:
        omega = 2 * np.pi * np.asarray(frequencies, float) / self.acq.sampling_rate_hz
        damp = self.acq.apodization_rate_per_s / self.acq.sampling_rate_hz
        return _geometric(-damp - 1j * omega, self.acq.n) / self.acq.n

    def render(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
               gain: complex = 1.0, phase_delay_s: float = 0.0, gaussian_sigma_hz: Rates = 0.0) -> np.ndarray:
        """Processed complex spectrum of the broadened transitions at `frequencies_hz`."""
        f = np.asarray(frequencies_hz, float)
        if not len(transitions):
            return np.zeros(len(f), complex)
        sigma = _per_transition(gaussian_sigma_hz, transitions, "Gaussian sigma")
        logz, c = self._modes(transitions, rates, gain, phase_delay_s)
        unstable = len(self.coef) and float(np.max(-logz.real)) * self.half > self.stability_limit
        if self.backend == "analytic" and np.any(sigma > 0):
            raise ValueError("Finite-record Gaussian broadening has no closed form here; use the 'time' or "
                             "'auto' backend, or the continuous route.")
        use_time = self.backend == "time" or (self.backend == "auto" and (
            unstable or np.any(sigma > 0) or self._time_is_cheaper(len(f), len(transitions), f)))
        if use_time:
            self.sampled_calls += 1
            with self.timer.section("render.time_domain"):
                fid = self.synthesize(transitions, rates, gain, phase_delay_s, gaussian_sigma_hz=gaussian_sigma_hz)
                return evaluate_spectrum(process_record(fid, self.acq), self.acq, f)
        return self._analytic(transitions, rates, f, [gain], phase_delay_s)[:, 0]

    def _analytic(self, transitions: TransitionList, rates: Rates, f: np.ndarray, gains, phase_delay_s: float):
        """Analytic processed spectra for several gains at once: shape (F, len(gains))."""
        logz, c_unit = self._modes(transitions, rates, 1.0, phase_delay_s)
        half = len(logz) // 2
        coeffs = []
        for g in gains:
            c = c_unit.copy()
            c[:half] *= g
            c[half:] *= np.conj(g)
            coeffs.append(c)
        return self.render_modes(logz, np.column_stack(coeffs), f)

    def render_modes(self, logz: np.ndarray, c: np.ndarray, f: np.ndarray) -> np.ndarray:
        """Processed spectra of x[m] = sum_j c[j, g] exp(logz_j m) for each column g: shape (F, G).

        The caller supplies modes whose sums are real signals (conjugate pairs or
        real exponentials). Used for transitions and for nuisance terms.
        """
        self.analytic_calls += 1
        f = np.asarray(f, float)
        c = np.asarray(c, complex)
        if c.ndim == 1:
            c = c[:, None]
        if len(self.coef) and float(np.max(-logz.real)) * self.half > self.stability_limit:
            m = np.arange(self.acq.points)
            signals = np.stack([(np.exp(np.outer(m, logz)) @ c[:, g]).real for g in range(c.shape[1])])
            return evaluate_spectrum(process_record(signals, self.acq), self.acq, f).T
        if len(self.coef):
            gainz = 1 - np.exp(np.outer(self.offsets, logz)).T @ self.coef
        else:
            gainz = np.ones(len(logz), complex)
        weights = c * gainz[:, None]
        spectrum = self._dtft_modes_multi(logz, weights, f)
        devs = [self._edge_deviation(logz, c[:, k]) for k in range(c.shape[1])]
        start_dev = np.column_stack([d[0] for d in devs]) if len(self.start_edge) else np.zeros((0, c.shape[1]))
        end_dev = np.column_stack([d[1] for d in devs]) if len(self.end_edge) else np.zeros((0, c.shape[1]))
        if len(self.start_edge):
            spectrum += self._dtft_samples(self.start_edge, start_dev, f)
        if len(self.end_edge):
            spectrum += self._dtft_samples(self.end_edge, end_dev, f)
        if self.acq.remove_mean:
            zero = np.zeros(1)
            mean = self._dtft_modes_multi(logz, weights, zero, apodize=False)[0]
            if len(self.start_edge):
                mean = mean + self._dtft_samples(self.start_edge, start_dev, zero, apodize=False)[0]
            if len(self.end_edge):
                mean = mean + self._dtft_samples(self.end_edge, end_dev, zero, apodize=False)[0]
            spectrum -= self._constant_dtft(f)[:, None] * mean[None, :]
        phase = self.acq.reference_phase(f)
        return spectrum * (phase[:, None] if np.ndim(phase) else phase)

    def _dtft_modes_multi(self, logz: np.ndarray, weights: np.ndarray, frequencies: np.ndarray,
                          apodize: bool = True) -> np.ndarray:
        """Like `_dtft_modes` for a (2T, G) weight matrix; returns (F, G)."""
        acq = self.acq
        omega = 2 * np.pi * np.asarray(frequencies, float) / acq.sampling_rate_hz
        scaled = weights * np.exp(logz * acq.start_sample)[:, None]
        damp = acq.apodization_rate_per_s / acq.sampling_rate_hz if apodize else 0.0
        out = np.empty((len(omega), weights.shape[1]), complex)
        for first in range(0, len(omega), self.chunk):
            L = logz[None, :] - damp - 1j * omega[first:first + self.chunk, None]
            out[first:first + self.chunk] = _geometric(L, acq.n) @ scaled
        return out / acq.n

    # -- time-domain reference ---------------------------------------------------
    def synthesize(self, transitions: TransitionList, rates: Rates, gain: complex = 1.0,
                   phase_delay_s: float = 0.0, include_dc: bool = False,
                   gaussian_sigma_hz: Rates = 0.0, method: str = "nufft") -> np.ndarray:
        """Real FID over the full record (real part of `synthesize_complex`)."""
        out = self.synthesize_complex(transitions, rates, gain, phase_delay_s, gaussian_sigma_hz, method).real
        if include_dc:
            out = out + np.real(transitions.dc * gain)
        return out

    def render_pair(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
                    phase_delay_s: float = 0.0, gaussian_sigma_hz: Rates = 0.0) -> np.ndarray:
        """Columns for real and imaginary parts of a complex gain: shape (F, 2).

        spectrum(g) = Re(g) * col0 + Im(g) * col1, used for variable projection.
        The time path synthesizes one complex signal z and processes Re z and -Im z.
        """
        f = np.asarray(frequencies_hz, float)
        sigma = _per_transition(gaussian_sigma_hz, transitions, "Gaussian sigma") if len(transitions) else np.zeros(0)
        if not len(transitions):
            return np.zeros((len(f), 2), complex)
        if self.backend == "analytic":
            if np.any(sigma > 0):
                raise ValueError("Finite-record Gaussian broadening needs the 'time' or 'auto' backend.")
            return self._analytic(transitions, rates, f, [1.0, 1j], phase_delay_s)
        logz, _ = self._modes(transitions, rates, 1.0, phase_delay_s)
        unstable = len(self.coef) and float(np.max(-logz.real)) * self.half > self.stability_limit
        if self.backend == "time" or unstable or np.any(sigma > 0) or (
                self.backend == "auto" and self._time_is_cheaper(len(f), len(transitions), f)):
            z = self.synthesize_complex(transitions, rates, 1.0, phase_delay_s, gaussian_sigma_hz)
            pair = process_record(np.stack([z.real, -z.imag]), self.acq)
            return evaluate_spectrum(pair, self.acq, f).T
        if unstable:
            z = self.synthesize_complex(transitions, rates, 1.0, phase_delay_s, gaussian_sigma_hz)
            pair = process_record(np.stack([z.real, -z.imag]), self.acq)
            return evaluate_spectrum(pair, self.acq, f).T
        return self._analytic(transitions, rates, f, [1.0, 1j], phase_delay_s)

    def synthesize_complex(self, transitions: TransitionList, rates: Rates, gain: complex = 1.0,
                           phase_delay_s: float = 0.0, gaussian_sigma_hz: Rates = 0.0,
                           method: str = "nufft") -> np.ndarray:
        """Complex analytic signal z(t_m) whose real part is the FID.

        x(t_m) = Re(sum_k A_k exp((2 pi i f_k - R_k) t_m - (2 pi sigma_k t_m)^2 / 2)),
        t_m = time_origin + m / fs. Transitions sharing (R, sigma) are summed with
        a type-1 NUFFT (`method="nufft"`, relative error about 1e-12) or directly
        (`method="direct"`, reference).
        """
        acq = self.acq
        out = np.zeros(acq.points, complex)
        if not len(transitions):
            return out
        f = transitions.frequencies_hz
        r = _per_transition(rates, transitions, "Decay rates")
        sigma = _per_transition(gaussian_sigma_hz, transitions, "Gaussian sigma")
        a = transitions.amplitudes * gain * np.exp(2j * np.pi * f * (phase_delay_s + acq.time_origin_s))
        times = acq.times()
        omega = 2 * np.pi * f / acq.sampling_rate_hz
        keys = np.stack([r, sigma], axis=1)
        unique, inverse = np.unique(keys, axis=0, return_inverse=True)
        inverse = np.asarray(inverse).ravel()
        for g, (rate, sig) in enumerate(unique):
            members = inverse == g
            if method == "nufft":
                with self.timer.section("render.nufft"):
                    oscillation = nufft_type1(omega[members], a[members], acq.points)
            elif method == "direct":
                m = np.arange(acq.points)
                oscillation = np.zeros(acq.points, complex)
                step = max(1, 2_000_000 // max(1, int(members.sum())))
                for first in range(0, acq.points, step):
                    block = m[first:first + step]
                    oscillation[first:first + step] = np.exp(1j * np.outer(block, omega[members])) @ a[members]
            else:
                raise ValueError("method must be 'nufft' or 'direct'.")
            envelope = np.exp(-rate * times - 0.5 * (2 * np.pi * sig * times) ** 2)
            out += envelope * oscillation
        return out


class ContinuousRenderer:
    """Pure route: infinite-record spectrum with no sampling or processing.

    X(f) = (1/T) integral_0^inf x(t) exp(-2 pi i f t) dt for the damped real FID,
    i.e. complex Lorentzians at +f_k and their mirror terms at -f_k. `T`
    (`normalization_s`) only sets the scale; with T equal to a record length
    the result approximates the finite-record renderer when the signal has
    decayed within that record. Rates must be positive.
    """

    def __init__(self, normalization_s: float = 1.0):
        if not normalization_s > 0:
            raise ValueError("normalization_s must be positive.")
        self.normalization_s = normalization_s

    @staticmethod
    def _line(a: np.ndarray, r: np.ndarray, sigma: np.ndarray, offset: np.ndarray) -> np.ndarray:
        """integral_0^inf a exp(-(r + 2 pi i offset) t - (2 pi sigma t)^2 / 2) dt, elementwise.

        Columns with sigma = 0 are Lorentzian; others use the Faddeeva function:
        integral = sqrt(pi / 2) / b * w(i k / (sqrt(2) b)), b = 2 pi sigma.
        """
        k = r + 2j * np.pi * offset
        out = np.empty(np.broadcast(a, k).shape, complex)
        voigt = sigma > 0
        lorentz = ~voigt
        out[:, lorentz] = (a[lorentz] / k[:, lorentz])
        if np.any(voigt):
            b = 2 * np.pi * sigma[voigt]
            out[:, voigt] = a[voigt] * np.sqrt(np.pi / 2) / b * wofz(1j * k[:, voigt] / (np.sqrt(2) * b))
        return out

    def render(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
               gain: complex = 1.0, phase_delay_s: float = 0.0, gaussian_sigma_hz: Rates = 0.0) -> np.ndarray:
        f = np.asarray(frequencies_hz, float)
        if not len(transitions):
            return np.zeros(len(f), complex)
        fk = transitions.frequencies_hz
        r = _per_transition(rates, transitions, "Decay rates")
        sigma = _per_transition(gaussian_sigma_hz, transitions, "Gaussian sigma")
        if np.any((r <= 0) & (sigma <= 0)):
            raise ValueError("Each transition needs a positive rate or a positive Gaussian sigma.")
        a = transitions.amplitudes * gain * np.exp(2j * np.pi * fk * phase_delay_s)
        out = np.empty(len(f), complex)
        step = max(1, 2_000_000 // max(1, len(fk)))
        for first in range(0, len(f), step):
            ff = f[first:first + step, None]
            positive = self._line(a / 2, r, sigma, ff - fk)
            negative = self._line(np.conj(a) / 2, r, sigma, ff + fk)
            out[first:first + step] = (positive + negative).sum(axis=1)
        return out / self.normalization_s

    def render_pair(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
                    phase_delay_s: float = 0.0, gaussian_sigma_hz: Rates = 0.0) -> np.ndarray:
        return np.column_stack([self.render(transitions, rates, frequencies_hz, 1.0, phase_delay_s, gaussian_sigma_hz),
                                self.render(transitions, rates, frequencies_hz, 1j, phase_delay_s, gaussian_sigma_hz)])
