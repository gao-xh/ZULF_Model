"""An averaged FID to the phased complex spectrum of a series entry (scripts/make_series_entry.py, ZULF Studio's
Process mode): the plan of the recipe, zulf_processing.process_dataset (the one processing operator), and the
phase of the instrument calibration at this FID's switching edge, or a phase and delay given by hand."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def series_spectrum(fid: np.ndarray, sampling_rate_hz: float, crop_s: float = 0.1, record_s: float = 8.0,
                    apodization_per_s: float = 0.3, zero_fill: int = 3, grid: Sequence[float] = (20.0, 380.0),
                    sg_window_s: Optional[float] = None, phase0_deg: Optional[float] = None,
                    delay_ms: Optional[float] = None) -> dict:
    """Processed spectrum of an averaged FID. Returns f (Hz), spectrum (complex, phase corrected), acquisition
    (dict), phasing {phase0_rad, delay_s}, edge_s (switching edge + calibration offset) and the plan.

    Phase: zero order `phase0_deg` (default: the calibration of configs/confirmed_samples.json), first order the
    delay `delay_ms` (default: minus the switching edge with the calibration offset). `sg_window_s` sets the drift
    filter length (default: the plan's samples)."""
    from zulf_core.render.phasing import correction_phasor, reference_delay_s
    from . import plan_for_dataset, process_dataset
    from .diagnostics import switching_edge
    cal = json.load(open(ROOT / "configs" / "confirmed_samples.json"))["processing"]["phase_calibration"]
    y = np.asarray(fid, float)
    fs = float(sampling_rate_hz)
    edge = switching_edge(y, fs)["edge_time_s"] + cal["delay_offset_s"]
    start = int(round(crop_s * fs))
    stop = min(start + int(round(record_s * fs)), len(y))
    settings = {"start_sample": start, "stop_sample": stop, "zero_fill": int(zero_fill),
                "apodization_rate_per_s": float(apodization_per_s), "ranges": [list(map(float, grid))]}
    if sg_window_s:
        settings["sg_window"] = max(int(round(sg_window_s * fs)) // 2 * 2 + 1, 5)
    plan = plan_for_dataset(len(y), fs, None, settings)
    ds = process_dataset(y, fs, plan=plan, phase_criterion=None)
    acq = plan.acquisition()
    f = ds.frequencies_hz
    phi0 = np.radians(cal["phase0_deg"] if phase0_deg is None else phase0_deg)
    delay = -(edge + acq.time_origin_s) if delay_ms is None else 1e-3 * float(delay_ms)
    r = ds.spectrum * correction_phasor(f, phi0, delay + reference_delay_s(acq))
    return {"f": np.asarray(f), "spectrum": np.asarray(r), "acquisition": acq.to_dict(), "plan": plan,
            "phasing": {"phase0_rad": float(phi0), "delay_s": float(delay)}, "edge_s": float(edge)}
