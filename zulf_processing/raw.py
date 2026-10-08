"""Reading raw data: averaged FIDs and instrument settings."""
from __future__ import annotations

import configparser
from pathlib import Path
from typing import Optional

import numpy as np


def load_fid(path) -> np.ndarray:
    """An averaged FID as float64 (.npy; other formats are added here as they appear)."""
    p = Path(path)
    if p.suffix == ".npy":
        return np.asarray(np.load(p), dtype=float)
    raise ValueError(f"Unsupported FID format: {p.suffix}")


def read_settings(path) -> dict:
    """NMRduino settings (.ini): device sample rate and length, pulse sequence, date, data directory. Missing keys
    are left out; the caller states any assumption it makes instead."""
    cp = configparser.ConfigParser(interpolation=None, strict=False)
    cp.optionxform = str
    cp.read(path)
    out = {}
    if cp.has_section("General"):
        out["date_time"] = cp.get("General", "DateTime", fallback=None)
    for section in cp.sections():
        if cp.get(section, "DeviceType", fallback="").startswith("Teensy") or section == "NMRduino":
            out["sampling_rate_hz"] = cp.getfloat(section, "SampleRate", fallback=None)
            out["points"] = cp.getint(section, "NumberOfSamples", fallback=None)
        if cp.get(section, "TaskType", fallback="") == "Pulse and Acquire":
            out["pulse_sequence"] = cp.get(section, "PulseSequencePath", fallback=None)
            out["scans"] = cp.getint(section, "NumberOfScans", fallback=None)
            out["data_directory"] = cp.get(section, "DataDirectory", fallback=None)
    return {k: v for k, v in out.items() if v is not None}


def find_sampling_rate(fid_path, default_hz: Optional[float] = None) -> tuple:
    """(sampling rate in Hz, source) of an averaged FID: scans.json of scripts/average_scans.py next to it, else the
    first .ini next to it, else `default_hz` (source "default"). Raises when nothing is found and there is no
    default. The NMRduino sequences differ (2000, 4000, 8333 Hz); a wrong rate scales every frequency."""
    folder = Path(fid_path).expanduser().resolve().parent
    record = folder / "scans.json"
    if record.exists():
        import json
        rate = json.loads(record.read_text()).get("sampling_rate_hz")
        if rate:
            return float(rate), str(record)
    for ini in sorted(folder.glob("*.ini")):
        rate = read_settings(ini).get("sampling_rate_hz")
        if rate:
            return float(rate), str(ini)
    if default_hz is None:
        raise ValueError(f"No sampling rate found next to {fid_path}; pass it explicitly.")
    return float(default_hz), "default"


def scale_sample_settings(settings: dict, sampling_rate_hz: float) -> dict:
    """Copy of processing defaults given in samples at settings["sampling_rate_hz"], rescaled to another rate so the
    same times result: start_sample, stop_sample, and sg_window (kept odd, so the drift filter keeps its cutoff in Hz)."""
    out = dict(settings)
    ratio = float(sampling_rate_hz) / float(settings["sampling_rate_hz"])
    out["sampling_rate_hz"] = float(sampling_rate_hz)
    if abs(ratio - 1.0) < 1e-9:
        return out
    for key in ("start_sample", "stop_sample"):
        if key in out:
            out[key] = int(round(out[key] * ratio))
    if "sg_window" in out:
        out["sg_window"] = int(round(out["sg_window"] * ratio)) // 2 * 2 + 1
    return out
