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
"""
from __future__ import annotations

from typing import Sequence, Union

import numpy as np

from ..physics.transitions import TransitionList
from .acquisition import Acquisition, evaluate_spectrum, process_record

Rates = Union[float, Sequence[float], np.ndarray]


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

    def __init__(self, acquisition: Acquisition, chunk_frequencies: int = 2048,
                 stability_limit: float = 8.0):
        self.acq = acquisition
        self.chunk = chunk_frequencies
        self.stability_limit = stability_limit
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

    # -- mode bookkeeping ----------------------------------------------------
    def _modes(self, transitions: TransitionList, rates: Rates, gain: complex, phase_delay_s: float):
        f = transitions.frequencies_hz
        if np.isscalar(rates):
            r = np.full(len(f), float(rates))
        else:
            table = np.asarray(rates, float)
            r = table[transitions.families]
        if np.any(r < 0) or not np.isfinite(r).all():
            raise ValueError("Decay rates must be finite and nonnegative.")
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

    def _dtft_modes(self, logz: np.ndarray, weights: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
        """(1/n) sum_{m=start}^{stop-1} sum_j weights_j z_j^m exp(-i w (m - start))."""
        acq = self.acq
        n = acq.n
        omega = 2 * np.pi * np.asarray(frequencies, float) / acq.sampling_rate_hz
        scaled = weights * np.exp(logz * acq.start_sample)
        out = np.empty(len(omega), complex)
        for first in range(0, len(omega), self.chunk):
            L = logz[None, :] - 1j * omega[first:first + self.chunk, None]
            out[first:first + self.chunk] = _geometric(L, n) @ scaled
        return out / n

    def _dtft_samples(self, sample_index: np.ndarray, values: np.ndarray, frequencies: np.ndarray) -> np.ndarray:
        if not len(sample_index):
            return np.zeros(len(frequencies), complex)
        rel = sample_index - self.acq.start_sample
        kernel = np.exp(-2j * np.pi * np.outer(frequencies, rel) / self.acq.sampling_rate_hz)
        return kernel @ values / self.acq.n

    def _constant_dtft(self, frequencies: np.ndarray) -> np.ndarray:
        omega = 2 * np.pi * np.asarray(frequencies, float) / self.acq.sampling_rate_hz
        return _geometric(-1j * omega, self.acq.n) / self.acq.n

    def render(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
               gain: complex = 1.0, phase_delay_s: float = 0.0) -> np.ndarray:
        """Processed complex spectrum of the damped transitions at `frequencies_hz`."""
        f = np.asarray(frequencies_hz, float)
        if not len(transitions):
            return np.zeros(len(f), complex)
        logz, c = self._modes(transitions, rates, gain, phase_delay_s)
        if len(self.coef) and float(np.max(-logz.real)) * self.half > self.stability_limit:
            self.sampled_calls += 1
            fid = self.synthesize(transitions, rates, gain, phase_delay_s)
            return evaluate_spectrum(process_record(fid, self.acq), self.acq, f)
        self.analytic_calls += 1
        if len(self.coef):
            gainz = 1 - np.exp(np.outer(self.offsets, logz)).T @ self.coef
        else:
            gainz = np.ones(len(logz), complex)
        weights = c * gainz
        spectrum = self._dtft_modes(logz, weights, f)
        start_dev, end_dev = self._edge_deviation(logz, c)
        spectrum += self._dtft_samples(self.start_edge, start_dev, f)
        spectrum += self._dtft_samples(self.end_edge, end_dev, f)
        if self.acq.remove_mean:
            zero = np.zeros(1)
            mean = (self._dtft_modes(logz, weights, zero) + self._dtft_samples(self.start_edge, start_dev, zero)
                    + self._dtft_samples(self.end_edge, end_dev, zero))[0]
            spectrum -= mean * self._constant_dtft(f)
        return spectrum

    def render_pair(self, transitions: TransitionList, rates: Rates, frequencies_hz: np.ndarray,
                    phase_delay_s: float = 0.0) -> np.ndarray:
        """Columns for real and imaginary parts of a complex gain: shape (F, 2).

        spectrum(g) = Re(g) * col0 + Im(g) * col1, used for variable projection.
        """
        col0 = self.render(transitions, rates, frequencies_hz, 1.0, phase_delay_s)
        col1 = self.render(transitions, rates, frequencies_hz, 1j, phase_delay_s)
        return np.column_stack([col0, col1])

    # -- time-domain reference ---------------------------------------------------
    def synthesize(self, transitions: TransitionList, rates: Rates, gain: complex = 1.0,
                   phase_delay_s: float = 0.0, include_dc: bool = False) -> np.ndarray:
        """Real FID over the full record (reference and noise injection)."""
        logz, c = self._modes(transitions, rates, gain, phase_delay_s)
        m = np.arange(self.acq.points)
        out = np.zeros(self.acq.points)
        step = max(1, 2_000_000 // max(1, len(logz)))
        for first in range(0, len(m), step):
            block = m[first:first + step]
            out[first:first + step] = (np.exp(np.outer(block, logz)) @ c).real
        if include_dc:
            out += np.real(transitions.dc * gain)
        return out
