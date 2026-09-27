"""Raw-FID diagnostics for one dataset (before any processing choice)."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional, Sequence, Tuple

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


def signal_extent(fid: np.ndarray, sampling_rate_hz: float, start_sample: int, sg_window: int = 201,
                  band_hz: Tuple[float, float] = (60.0, 330.0), exclude_hz: Sequence[float] = (),
                  block_s: float = 0.25, floor_fraction: float = 0.3, sustain_s: float = 1.0) -> dict:
    """How long the signal lasts in this dataset, from the data alone.

    The record after `start_sample` is high-passed with the same SG window as the processing and cut into blocks;
    in each block the power inside the signal bands (found in the first second, instrument lines excluded) is
    compared with the same bands in the last `floor_fraction` of the record (noise floor, same spectral colour).
    Returns `end_s`, the signal end: the first block after which the power stays below 3x the floor for
    `sustain_s` (single noise blocks of a few bins fluctuate up to about 2.5x); the times (s after the crop start) at which the cumulative excess signal energy reaches 90 / 95 / 99 %,
    the last block with power above 3x the floor, the energy fraction inside 1 s, and the per-block ratios.
    Frequency information grows with t^2, so the late few percent of the energy matter for narrow lines and small
    couplings (ANALYSIS_LOG, window length)."""
    from scipy.signal import savgol_filter
    from .phase import signal_regions_hz
    x = np.asarray(fid, float)
    fs = float(sampling_rate_hz)
    y = x - savgol_filter(x, sg_window, 2, mode="mirror") if sg_window else x - x.mean()
    y = y[int(start_sample):]
    block = max(8, int(round(block_s * fs)))
    nb = len(y) // block
    if nb < 8:
        raise ValueError("Record too short for a signal-extent estimate.")
    first = y[:min(len(y), int(fs))]
    ff = np.fft.rfftfreq(4 * len(first), 1 / fs)
    spec = np.fft.rfft(first, 4 * len(first))
    m = (ff > band_hz[0]) & (ff < band_hz[1])
    regions = signal_regions_hz(spec[m], ff[m], exclude_hz=exclude_hz)
    f = np.fft.rfftfreq(block, 1 / fs)
    sig = np.zeros(len(f), bool)
    for lo, hi in regions:
        sig |= (f >= lo - 1.0) & (f <= hi + 1.0)
    for line in exclude_hz:
        sig &= np.abs(f - line) > 2.0
    if not sig.any():
        raise ValueError("No signal bands for the signal-extent estimate.")
    power = np.array([np.sum(np.abs(np.fft.rfft(y[k * block:(k + 1) * block])[sig]) ** 2) for k in range(nb)])
    floor = float(np.median(power[int(nb * (1 - floor_fraction)):]))
    excess = power - floor              # not clipped: noise-only blocks average to zero instead of adding up
    total = excess.sum()
    cum = np.cumsum(excess) / total if total > 0 else np.ones(nb)
    at = lambda q: float((np.searchsorted(cum, q) + 1) * block / fs)
    above = [k for k in range(nb) if power[k] > 3 * floor]
    hold = max(1, int(round(sustain_s * fs / block)))
    end = next((k for k in range(nb - hold + 1) if np.all(power[k:k + hold] <= 3 * floor)), nb)
    return {"end_s": float(end * block / fs),"t90_s": at(0.90), "t95_s": at(0.95), "t99_s": at(0.99),
            "last_above_3x_floor_s": float((max(above) + 1) * block / fs) if above else 0.0,
            "energy_within_1s": float(cum[max(0, int(round(fs / block)) - 1)]),
            "block_s": block / fs, "ratio_to_floor": [round(float(p / floor), 2) for p in power],
            "bands_hz": regions}


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
