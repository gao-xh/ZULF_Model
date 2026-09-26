"""Type-1 non-uniform FFT by Gaussian gridding (Greengard and Lee, SIAM Rev. 2004).

Computes s[m] = sum_k a_k exp(i w_k m) for m = 0 .. N-1 and arbitrary angular
frequencies w_k (radians per sample) in O(K * spread + R N log(R N)) instead of
O(K N). With oversampling R = 2 and `spread` = 12 the relative error is about
1e-12; tests compare against the direct sum.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import ifft, next_fast_len


def nufft_type1(omega: np.ndarray, amplitudes: np.ndarray, n: int, oversampling: int = 2,
                spread: int = 12) -> np.ndarray:
    omega = np.mod(np.asarray(omega, float), 2 * np.pi)
    a = np.asarray(amplitudes, complex)
    if omega.shape != a.shape or omega.ndim != 1:
        raise ValueError("omega and amplitudes must be matching 1D arrays.")
    if not len(a):
        return np.zeros(n, complex)
    # Centre the output range: m = m' + c with m' in [-n/2, n/2).
    c = n // 2
    a = a * np.exp(1j * omega * c)
    # Oversampled grid rounded up to a fast FFT length (record lengths often have large prime factors);
    # the Gaussian width uses the effective oversampling ratio.
    m_os = next_fast_len(int(oversampling * n))
    ratio = m_os / n
    h = 2 * np.pi / m_os
    tau = np.pi * spread / (n ** 2 * ratio * (ratio - 0.5))
    j0 = np.floor(omega / h).astype(int)
    offsets = np.arange(-spread + 1, spread + 1)
    j = j0[:, None] + offsets[None, :]
    distance = j * h - omega[:, None]
    weights = np.exp(-distance ** 2 / (4 * tau)) * a[:, None]
    grid = np.zeros(m_os, complex)
    np.add.at(grid, np.mod(j, m_os).ravel(), weights.ravel())
    transformed = ifft(grid) * m_os  # sum_j grid_j exp(+i w_j m)
    m_prime = np.arange(n) - c
    values = transformed[np.mod(m_prime, m_os)] * h
    return values / (np.sqrt(4 * np.pi * tau) * np.exp(-tau * m_prime.astype(float) ** 2))
