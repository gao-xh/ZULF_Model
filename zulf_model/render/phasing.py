"""Zero- and first-order phase correction of complex spectra.

Convention: a line with complex amplitude a (signal Re(a exp(2 pi i f t)))
appears in the processed spectrum multiplied by exp(i (phase0 + 2 pi f delay)),
where `delay` is the physical delay between the spin evolution origin and the
acquisition time origin plus the Fourier reference shift of the crop
(`reference_delay_s`). Correction multiplies by the inverse phasor, after which
the real part is the absorption spectrum with the physical line signs.

For experimental data the operator supplies `phase0_rad` and the instrument
`delay_s` (the first-order phase expressed as a time); the crop contribution is
added from the acquisition, so the same values apply to any crop.
"""
from __future__ import annotations

import numpy as np

from .acquisition import Acquisition


def reference_delay_s(acquisition: Acquisition) -> float:
    """First-order phase (as a time) introduced by the Fourier reference of the processed record."""
    if acquisition.phase_reference == "crop_start":
        return acquisition.time_origin_s + acquisition.start_sample / acquisition.sampling_rate_hz
    return 0.0


def correction_phasor(frequencies_hz: np.ndarray, phase0_rad: float, delay_s: float) -> np.ndarray:
    """exp(-i (phase0 + 2 pi f delay)); multiply a spectrum by this to remove that phase."""
    f = np.asarray(frequencies_hz, float)
    return np.exp(-1j * (float(phase0_rad) + 2 * np.pi * f * float(delay_s)))


def phase_correct(values: np.ndarray, frequencies_hz: np.ndarray, phase0_rad: float, delay_s: float,
                  acquisition: Acquisition = None) -> np.ndarray:
    """Apply a zero/first-order correction; includes the crop reference when an acquisition is given."""
    total = float(delay_s) + (reference_delay_s(acquisition) if acquisition is not None else 0.0)
    return np.asarray(values, complex) * correction_phasor(frequencies_hz, phase0_rad, total)


def phase1_to_delay_s(phase1_deg: float, spectral_width_hz: float) -> float:
    """Convert a first-order phase given in degrees across a spectral width into a delay in seconds."""
    return float(phase1_deg) / 360.0 / float(spectral_width_hz)
