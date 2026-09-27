"""Raw-FID diagnostics for one dataset (before any processing choice)."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional, Tuple

import numpy as np

from zulf_core.diagnostics import diagnose_fid


def switching_edge(fid: np.ndarray, sampling_rate_hz: float, search_s: Tuple[float, float] = (0.0005, 0.02),
                   settle_samples: Tuple[int, int] = (2, 10)) -> dict:
    """Time of the field switch-off edge in the raw FID (half height between the start plateau and the first
    extremum inside `search_s`, linear interpolation between samples).

    NMRduino with standard_zf_4000Hz_no_dead: plateau to 2.75 ms, edge at 3.41-3.51 ms in all confirmed
    datasets, then ringing to about 50 ms. The zero-field evolution starts at the edge, so minus the edge time is
    the delay of the processed spectrum in the convention of `zulf_core.render.phasing.phase_correct`. Checked
    against complex fits (ANALYSIS_LOG): within 0.05 ms for alanine, triethylamine and N-ethylmethylamine; lactic
    acid fits equally well with the edge value; fitted delays of narrow-band spectra absorb model error."""
    x = np.asarray(fid, float)
    fs = float(sampling_rate_hz)
    a, b = settle_samples
    plateau = float(np.median(x[a:b]))
    lo_i, hi_i = int(search_s[0] * fs), int(search_s[1] * fs)
    window = x[lo_i:hi_i]
    extreme = int(np.argmax(np.abs(window - plateau)))
    target = float(window[extreme])
    half = plateau + 0.5 * (target - plateau)
    sign = np.sign(target - plateau)
    cross = next((k for k in range(max(lo_i, 1), lo_i + extreme + 1) if sign * (x[k] - half) >= 0), None)
    if cross is None:
        raise ValueError("No switching edge found in the raw FID.")
    t = (cross - 1 + (half - x[cross - 1]) / (x[cross] - x[cross - 1])) / fs
    return {"edge_time_s": float(t), "plateau": plateau, "edge_amplitude": target - plateau,
            "method": "half height of the field switch-off edge"}


@dataclass
class RawDiagnostics:
    sampling_rate_hz: float
    points: int
    edge_time_s: Optional[float]
    plateau_end_s: float
    ringing_end_s: float
    late_noise_rms: float
    first_point_anomaly: bool
    notes: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def diagnose_raw(fid: np.ndarray, sampling_rate_hz: float) -> RawDiagnostics:
    """Switching edge plus the zulf_core FID diagnostics (plateau, ringing end, late noise)."""
    # the multi-exponential fits of diagnose_fid are not used here and cost minutes per FID: skipped
    d = diagnose_fid(fid, sampling_rate_hz, max_exponentials=0)
    notes = list(d.notes)
    try:
        edge = switching_edge(fid, sampling_rate_hz)["edge_time_s"]
    except ValueError as exc:
        edge = None
        notes.append(str(exc))
    return RawDiagnostics(float(sampling_rate_hz), len(fid), edge, float(d.plateau_end_s), float(d.ringing_end_s),
                          float(d.late_noise_rms), bool(d.first_point_anomaly), notes)
