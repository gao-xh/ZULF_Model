"""Common scoring of refined hypotheses: one region, one noise level, one criterion.

Refinements weight each model differently (their own model cores), so the
ranking uses a fixed yardstick instead:

* region: the data cores of signal weighting (weight >= 0.6 around lines
  found in the data alone), the same points for every hypothesis;
* noise: robust (MAD) scale of the real and imaginary parts outside the
  cores after removing a running median, i.e. per-component complex noise;
* decorrelation: every `zero_fill`-th point only (zero-filled points are
  correlated), real and imaginary parts counted separately (real only for a
  phased real spectrum);
* criterion: BIC = chi2 + k ln N (default) or AIC = chi2 + 2 k, with k the
  free nonlinear parameters plus the free amplitude degrees of freedom and
  the background terms.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.ndimage import median_filter

from zulf_core.solver.forward import signal_regions, signal_weights


@dataclass
class Yardstick:
    mask: np.ndarray        # points used (data cores, decimated)
    sigma: float            # per-component complex noise
    n: int                  # independent real numbers compared
    real_only: bool

    def chi2(self, observed_values: np.ndarray, prediction: np.ndarray) -> float:
        r = (observed_values - prediction)[self.mask]
        if self.real_only:
            return float(np.sum(r.real ** 2) / self.sigma ** 2)
        return float(np.sum(r.real ** 2 + r.imag ** 2) / self.sigma ** 2)


def yardstick(observed, threshold: float = 4.0, taper_hz: float = 2.0, outside: float = 0.2,
              baseline_hz: float = 8.0) -> Yardstick:
    sel = observed.selected
    f, y, band = observed.frequencies_hz[sel], observed.values[sel], observed.band_index[sel]
    cores, _ = signal_regions(f, y, band, threshold, baseline_hz)
    weights = signal_weights(f, band, cores, outside, taper_hz) if cores.any() else np.full(len(f), outside)
    region = weights >= 0.6
    quiet = weights < 0.3
    zero_fill = int((observed.metadata or {}).get("zero_fill", 1)) or 1
    spacing = float(np.median(np.diff(f))) if len(f) > 1 else 1.0
    width = max(3, int(round(baseline_hz / spacing)) | 1)
    parts = [y.real] if observed.real_only else [y.real, y.imag]
    residuals = []
    for part in parts:
        for b in np.unique(band):
            m = band == b
            smooth = median_filter(part[m], size=min(width, int(m.sum()) | 1), mode="nearest")
            residuals.append((part[m] - smooth)[quiet[m]])
    pooled = np.concatenate(residuals) if residuals else np.zeros(1)
    sigma = max(1.4826 * float(np.median(np.abs(pooled - np.median(pooled)))), 1e-300)
    decimate = np.zeros(len(f), bool)
    decimate[::zero_fill] = True
    mask = region & decimate
    n = int(mask.sum()) * (1 if observed.real_only else 2)
    return Yardstick(mask, sigma, n, bool(observed.real_only))


def free_parameter_count(model, settings) -> int:
    """Free nonlinear parameters + amplitude degrees of freedom + background terms."""
    param = settings.parameterize(model.interpretation)
    n_nonlinear = len(param.free_names)
    if settings.amplitude_map is not None:
        n_amp = len(settings.amplitude_map[0])
    elif settings.amplitude_ratios is not None:
        n_amp = 1
    else:
        n_amp = len(model.component_labels)
    amp_dof = 2 * n_amp if settings.gain_model == "complex" else n_amp + 1
    return n_nonlinear + amp_dof


def criterion(chi2: float, k: int, n: int, kind: str = "bic") -> float:
    if kind == "aic":
        return chi2 + 2.0 * k
    return chi2 + k * math.log(max(n, 2))
