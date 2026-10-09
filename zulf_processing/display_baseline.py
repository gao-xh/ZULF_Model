"""Display-only baseline of a spectrum and its model (figures, ZULF Studio's "baseline corrected" view).

Never used before a fit: the crop gives every line oscillating wings and a rolling baseline, the fit renders the
model through the same processing and solves a background polynomial per band (WORKFLOW W2; AsLS on the data
before a fit biased 1J(C4,H4) of pyridine by -3.5 Hz). For display the same steps act on data and model alike:
(1) a cubic spline through anchor points away from the model lines, the power-line harmonics and narrow data
peaks near lines; (2) within lift_hz of the model lines an AsLS baseline (1.5 Hz) that lifts the negative crop
wings between close lines. The residual data - model changes by the difference of the two baselines only.

method (owner, 2026-10-09: the separate baselines made data and model differ where they agreed before):
- "separate": the steps above on data and on model, each its own baseline (the earlier figures);
- "shared" (default since 2026-10-09): one baseline from the data, subtracted from data and model alike: the
  residual is unchanged (ethanol 118-125 Hz: "separate" left the crop's rolling baseline in the data and took it
  out of the model, which looked like a missing line; with "shared" the two agree, as without a baseline);
- "residual": one baseline from the residual data - model (spline through anchors away from the lines),
  subtracted from the data only: the model is shown as the fit made it, the data lose only what varies slowly
  where the model has no line;
- "none": data and model as they are.
"""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

MAINS_HZ = 60.06


def display_baseline(f, y, m, lines, protect_hz: float = 1.0, knots_hz: float = 1.5,
                     lift_hz: float = 2.5, method: str = "shared") -> Tuple[np.ndarray, np.ndarray, list]:
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
    if method == "none":
        return y.copy(), m.copy(), mains
    if method == "residual":
        b = anchor_spline_baseline(y - m, f, protect, knot_spacing_hz=knots_hz, k_sigma=np.inf)
        return y - b, m.copy(), mains
    taper = np.clip(gaussian_filter1d(line_mask(f, model_lines, lift_hz).astype(float), 1.0 / step), 0, 1)
    b = anchor_spline_baseline(y, f, protect, knot_spacing_hz=knots_hz)
    b = b + taper * asls_baseline(y - b, f, smooth_hz=1.5, p=0.01)
    if method == "shared":
        return y - b, m - b, mains
    if method != "separate":
        raise ValueError(f"unknown baseline method {method!r} (separate, shared, residual, none)")
    yc = y - b
    mc = m - anchor_spline_baseline(m, f, protect, knot_spacing_hz=knots_hz, k_sigma=np.inf)
    mc = mc - taper * asls_baseline(mc, f, smooth_hz=1.5, p=0.01)
    return yc, mc, mains
