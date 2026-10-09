"""Display-only baseline of a spectrum and its model (figures, ZULF Studio's "baseline corrected" view).

Never used before a fit: the crop gives every line oscillating wings and a rolling baseline, the fit renders the
model through the same processing and solves a background polynomial per band (WORKFLOW W2; AsLS on the data
before a fit biased 1J(C4,H4) of pyridine by -3.5 Hz). For display the same steps act on data and model alike:
(1) a cubic spline through anchor points away from the model lines, the power-line harmonics and narrow data
peaks near lines; (2) within lift_hz of the model lines an AsLS baseline (1.5 Hz) that lifts the negative crop
wings between close lines. The residual data - model changes by the difference of the two baselines only.
"""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

MAINS_HZ = 60.06


def display_baseline(f, y, m, lines, protect_hz: float = 1.0, knots_hz: float = 1.5,
                     lift_hz: float = 2.5) -> Tuple[np.ndarray, np.ndarray, list]:
    """Baseline-corrected data and model. `lines`: [(component, frequencies (Hz), relative amplitudes), ...]
    (lines below 3 % of their component and below 20 Hz are not protected). Returns (data, model, mains)."""
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    from . import anchor_spline_baseline, asls_baseline, line_mask
    f, y, m = (np.asarray(v, float) for v in (f, y, m))
    model_lines = [v for _, fr, a in lines for v in np.asarray(fr)[np.asarray(a) >= 0.03] if v > 20.0]
    step = float(np.median(np.diff(f)))
    quiet = (f > f.max() - 10) & (f <= f.max())
    noise = 1.4826 * np.median(np.abs(y[quiet] - np.median(y[quiet])))
    peaks, _ = find_peaks(np.abs(y), prominence=6 * noise, width=(None, 1.5 / step))
    mains = [MAINS_HZ * k for k in range(1, int(f.max() / MAINS_HZ) + 1)]
    protect = (line_mask(f, model_lines, protect_hz) | line_mask(f, mains, 0.6)
               | (line_mask(f, f[peaks], 0.8) & line_mask(f, model_lines, 5.0)))
    yc = y - anchor_spline_baseline(y, f, protect, knot_spacing_hz=knots_hz)
    mc = m - anchor_spline_baseline(m, f, protect, knot_spacing_hz=knots_hz, k_sigma=np.inf)
    taper = np.clip(gaussian_filter1d(line_mask(f, model_lines, lift_hz).astype(float), 1.0 / step), 0, 1)
    yc = yc - taper * asls_baseline(yc, f, smooth_hz=1.5, p=0.01)
    mc = mc - taper * asls_baseline(mc, f, smooth_hz=1.5, p=0.01)
    return yc, mc, mains
