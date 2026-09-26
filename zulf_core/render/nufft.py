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
    """s[m] = sum_k a_k exp(i w_k m); `amplitudes` of shape (K,) or (K, P) gives (n,) or (n, P)."""
    omega = np.mod(np.asarray(omega, float), 2 * np.pi)
    a = np.asarray(amplitudes, complex)
    single = a.ndim == 1
    if single:
        a = a[:, None]
    if omega.ndim != 1 or a.ndim != 2 or a.shape[0] != omega.shape[0]:
        raise ValueError("omega must be 1D and amplitudes (K,) or (K, P) with matching K.")
    if not len(a):
        out = np.zeros((n, a.shape[1]), complex)
        return out[:, 0] if single else out
    # Centre the output range: m = m' + c with m' in [-n/2, n/2).
    c = n // 2
    a = a * np.exp(1j * omega * c)[:, None]
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
    kernel = np.exp(-distance ** 2 / (4 * tau))
    index = np.mod(j, m_os).ravel()
    grid = np.zeros((m_os, a.shape[1]), complex)
    for col in range(a.shape[1]):
        grid[:, col] = np.bincount(index, (kernel * a[:, col:col + 1]).real.ravel(), minlength=m_os) + \
            1j * np.bincount(index, (kernel * a[:, col:col + 1]).imag.ravel(), minlength=m_os)
    transformed = ifft(grid, axis=0) * m_os  # sum_j grid_j exp(+i w_j m)
    m_prime = np.arange(n) - c
    values = transformed[np.mod(m_prime, m_os)] * h
    out = values / (np.sqrt(4 * np.pi * tau) * np.exp(-tau * m_prime.astype(float) ** 2))[:, None]
    return out[:, 0] if single else out
