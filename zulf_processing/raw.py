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
