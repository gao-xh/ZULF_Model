"""Spectrum baselines by asymmetric least squares (AsLS, Eilers and Boelens 2005).

A smooth curve z under a real spectrum y minimises

    sum_i w_i (y_i - z_i)^2 + lam sum_i (z_{i-1} - 2 z_i + z_{i+1})^2,

re-weighting until the weights stop changing. Standard AsLS: w = p where y > z, 1 - p elsewhere (lines above
the baseline barely pull on it). Two-sided: points farther than about k noise sigma from z on either side get a
small weight, for spectra whose lines have both signs.

The smoothness is given as a length in Hz: lam = (smooth_hz / step)^4, so the same setting acts alike on any
grid (the second-difference penalty scales as step^-4). The baseline follows structure longer than about
pi * smooth_hz (the cutoff of the Whittaker smoother); broad features of a few hertz need smooth_hz of 1-2.
Contiguous frequency segments are treated separately.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve


def _segments(f: np.ndarray):
    step = float(np.median(np.diff(f)))
    breaks = np.flatnonzero(np.diff(f) > 1.5 * step) + 1
    return step, np.split(np.arange(len(f)), breaks)


def _penalty(n: int, lam: float):
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    return lam * (d.T @ d)


def asls_baseline(values: np.ndarray, frequencies_hz: np.ndarray, smooth_hz: float = 3.0, p: float = 0.01,
                  iterations: int = 50) -> np.ndarray:
    """Standard AsLS baseline of a real spectrum (lines assumed positive)."""
    y = np.asarray(values, float)
    f = np.asarray(frequencies_hz, float)
    z = y.copy()
    if len(f) < 4:
        return z
    step, segments = _segments(f)
    lam = (smooth_hz / step) ** 4
    for seg in segments:
        if len(seg) < 4:
            continue
        ys = y[seg]
        pen = _penalty(len(seg), lam)
        w = np.ones(len(seg))
        for _ in range(iterations):
            zs = spsolve(sparse.csc_matrix(sparse.diags(w) + pen), w * ys)
            w_new = np.where(ys > zs, p, 1.0 - p)
            if np.array_equal(w_new, w):
                break
            w = w_new
        z[seg] = zs
    return z


def asls_two_sided(values: np.ndarray, frequencies_hz: np.ndarray, sigma: float, smooth_hz: float = 3.0,
                   k: float = 3.0, iterations: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """Two-sided AsLS for lines of either sign: weight 1 / (1 + exp(2 (|y - z| - k sigma) / sigma)). Returns the
    baseline and the final weights."""
    y = np.asarray(values, float)
    f = np.asarray(frequencies_hz, float)
    z = y.copy()
    weights = np.ones(len(y))
    if len(f) < 4:
        return z, weights
    step, segments = _segments(f)
    lam = (smooth_hz / step) ** 4
    for seg in segments:
        if len(seg) < 4:
            continue
        ys = y[seg]
        pen = _penalty(len(seg), lam)
        w = np.ones(len(seg))
        zs = ys
        for _ in range(iterations):
            z_new = spsolve(sparse.csc_matrix(sparse.diags(w) + pen), w * ys)
            r = np.abs(ys - z_new)
            w = 1.0 / (1.0 + np.exp(np.clip(2.0 * (r - k * sigma) / sigma, -50, 50)))
            done = np.max(np.abs(z_new - zs)) < 1e-3 * sigma
            zs = z_new
            if done:
                break
        z[seg] = zs
        weights[seg] = w
    return z, weights
