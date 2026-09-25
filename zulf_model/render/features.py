"""Model input features from complex spectra."""
from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

CHANNELS = ("real", "imag", "magnitude")


def spectrum_scale(values: np.ndarray, method: str = "rms") -> float:
    v = np.abs(np.asarray(values))
    if method == "rms":
        s = float(np.sqrt(np.mean(v ** 2)))
    elif method == "max":
        s = float(v.max())
    else:
        raise ValueError("scale method must be 'rms' or 'max'.")
    return s if s > 0 else 1.0


def spectrum_features(values: np.ndarray, channels: Sequence[str] = ("real", "imag"),
                      scale_method: str = "rms") -> Tuple[np.ndarray, float]:
    """Return (float32 array of shape (channels, F), scale)."""
    values = np.asarray(values, complex)
    scale = spectrum_scale(values, scale_method)
    z = values / scale
    out = []
    for name in channels:
        if name == "real":
            out.append(z.real)
        elif name == "imag":
            out.append(z.imag)
        elif name == "magnitude":
            out.append(np.abs(z))
        else:
            raise ValueError(f"Unknown channel '{name}'. Known: {CHANNELS}")
    return np.stack(out).astype(np.float32), scale
