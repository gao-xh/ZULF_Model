"""Loading experimental FIDs.

The NMRduino decoder reproduces the legacy convention of the earlier projects
(reverse bytes, little-endian int16, keep [20:-2], reverse samples). It is a
reproduction of that convention, not a verified hardware specification.
Experimental folders are read-only; nothing here writes next to the data.
"""
from __future__ import annotations

import configparser
import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np

DECODER = "legacy_nmrduino_v1: reverse bytes; little-endian int16; [20:-2]; reverse samples"


def decode_dat(payload: bytes) -> np.ndarray:
    if len(payload) % 2 or len(payload) < 2 * 54:
        raise ValueError("Invalid DAT length.")
    return np.frombuffer(payload[::-1], dtype="<i2")[20:-2][::-1].astype(np.float64)


def read_ini(path) -> dict:
    """Return sampling rate, declared sample count and the parsed NMRduino section."""
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string(text)
    section = next((s for s in parser.sections() if s.lower() == "nmrduino"), None)
    if section is None:
        raise ValueError(f"{path}: missing NMRduino section.")
    rate = float(parser.get(section, "SampleRate"))
    declared = parser.getint(section, "NumberOfSamples", fallback=0)
    return {"sampling_rate_hz": rate, "declared_points": declared, "section": dict(parser.items(section)),
            "sha256": hashlib.sha256(raw).hexdigest()}


@dataclass
class ExperimentFID:
    fid: np.ndarray
    sampling_rate_hz: float
    scans: int
    source: dict = field(default_factory=dict)

    @property
    def points(self) -> int:
        return len(self.fid)

    @property
    def duration_s(self) -> float:
        return self.points / self.sampling_rate_hz


def load_average(npy_path, ini_path, scans: int = 0) -> ExperimentFID:
    fid = np.load(npy_path, allow_pickle=False).astype(float)
    if fid.ndim != 1 or not np.isfinite(fid).all():
        raise ValueError("Averaged FID must be a finite 1D array.")
    info = read_ini(ini_path)
    return ExperimentFID(fid, info["sampling_rate_hz"], scans,
                         {"npy": str(npy_path), "npy_sha256": hashlib.sha256(Path(npy_path).read_bytes()).hexdigest(),
                          "ini": str(ini_path), "ini_sha256": info["sha256"], "declared_points": info["declared_points"],
                          "decoder": DECODER})


def average_folder(folder, scan_ids: Optional[Sequence[int]] = None, groups: Optional[Sequence[Sequence[int]]] = None):
    """Stream numbered DAT files into one mean or several disjoint group means."""
    folder = Path(folder)
    info = read_ini(folder / "0.ini")
    files = {int(p.stem): p for p in folder.glob("*.dat") if p.stem.isdecimal()}
    if groups is None:
        groups = [sorted(files) if scan_ids is None else list(scan_ids)]
    flat = [i for g in groups for i in g]
    if len(set(flat)) != len(flat):
        raise ValueError("Groups must be disjoint.")
    means: List[np.ndarray] = []
    for group in groups:
        total, count = None, 0
        for scan in group:
            values = decode_dat(files[scan].read_bytes())
            total = values.copy() if total is None else total + values
            count += 1
        means.append(total / count)
    return [ExperimentFID(m, info["sampling_rate_hz"], len(g), {"folder": str(folder), "decoder": DECODER})
            for m, g in zip(means, groups)]
