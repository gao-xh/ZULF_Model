"""Single source of problem dimensions.

Every module that needs a spin count, nucleus list, component limit, J bin
layout or spectral range reads it from a `ProblemSpec`. Extending the problem
(more spins, other nuclei, more components) is a configuration change:
edit `configs/problem_v1.json` or pass a different file.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Tuple

import numpy as np

from .nuclei import get_registry


@dataclass(frozen=True)
class JBinSpec:
    """Symmetric nonuniform J bins.

    `segments` is a list of `[abs_upper_hz, width_hz]` with increasing upper
    limits. Bins are fine near zero and coarse at large |J|. The outermost
    limit is the largest representable |J|.
    """
    segments: Tuple[Tuple[float, float], ...] = ((20.0, 0.25), (60.0, 0.5), (200.0, 1.0), (300.0, 2.0))

    def edges(self) -> np.ndarray:
        positive = [0.0]
        for upper, width in self.segments:
            if width <= 0 or upper <= positive[-1]:
                raise ValueError("J bin segments need increasing limits and positive widths.")
            count = int(round((upper - positive[-1]) / width))
            if count < 1 or not np.isclose(positive[-1] + count * width, upper):
                raise ValueError("Each J segment length must be an integer multiple of its width.")
            positive.extend((positive[-1] + width * np.arange(1, count + 1)).tolist())
        positive = np.asarray(positive)
        return np.concatenate([-positive[:0:-1], positive])

    @property
    def max_abs_hz(self) -> float:
        return float(self.segments[-1][0])


@dataclass(frozen=True)
class GridSpec:
    """Model input spectral range. Resolution comes from the acquisition."""
    f_min_hz: float = 1.0
    f_max_hz: float = 400.0
    zero_fill: int = 1
    channels: Tuple[str, ...] = ("real", "imag")


@dataclass(frozen=True)
class ContributionSpec:
    """Log-spaced bins for relative component contributions (tokens)."""
    log10_min: float = -4.0
    log10_max: float = 0.0
    bins: int = 17


@dataclass(frozen=True)
class ProblemSpec:
    nuclei: Tuple[str, ...] = ("1H", "13C", "15N")
    spin_counts: Tuple[int, ...] = (8,)
    max_components: int = 4
    max_group_size: int = 6
    j_bins: JBinSpec = field(default_factory=JBinSpec)
    grid: GridSpec = field(default_factory=GridSpec)
    contribution: ContributionSpec = field(default_factory=ContributionSpec)
    name: str = "problem-v1"

    def __post_init__(self):
        registry = get_registry()
        if not self.nuclei or len(set(self.nuclei)) != len(self.nuclei):
            raise ValueError("nuclei must be a nonempty list of distinct symbols.")
        for symbol in self.nuclei:
            registry[symbol]
        if not self.spin_counts or any(type(n) is not int or n < 1 for n in self.spin_counts):
            raise ValueError("spin_counts must contain positive integers.")
        if type(self.max_components) is not int or self.max_components < 1:
            raise ValueError("max_components must be a positive integer.")
        if type(self.max_group_size) is not int or self.max_group_size < 1:
            raise ValueError("max_group_size must be a positive integer.")
        if not 0 <= self.grid.f_min_hz < self.grid.f_max_hz:
            raise ValueError("Grid needs 0 <= f_min < f_max.")
        if type(self.grid.zero_fill) is not int or self.grid.zero_fill < 1:
            raise ValueError("zero_fill must be a positive integer.")
        self.j_bins.edges()

    @property
    def max_spins(self) -> int:
        return max(self.spin_counts)

    @property
    def max_pairs(self) -> int:
        return self.max_spins * (self.max_spins - 1) // 2

    def nucleus_index(self, symbol: str) -> int:
        return self.nuclei.index(symbol)

    def allows(self, isotopes) -> bool:
        return len(isotopes) in self.spin_counts and all(s in self.nuclei for s in isotopes)

    def to_dict(self) -> dict:
        return asdict(self)

    def digest(self) -> str:
        return hashlib.sha256(json.dumps(self.to_dict(), sort_keys=True).encode()).hexdigest()[:16]

    @classmethod
    def from_dict(cls, data: dict) -> "ProblemSpec":
        data = dict(data)
        allowed = set(cls.__dataclass_fields__)
        if set(data) - allowed:
            raise ValueError(f"Unknown ProblemSpec keys: {sorted(set(data) - allowed)}")
        if "j_bins" in data:
            data["j_bins"] = JBinSpec(tuple(tuple(s) for s in data["j_bins"]["segments"]))
        if "grid" in data:
            g = dict(data["grid"])
            if "channels" in g:
                g["channels"] = tuple(g["channels"])
            data["grid"] = GridSpec(**g)
        if "contribution" in data:
            data["contribution"] = ContributionSpec(**data["contribution"])
        for key in ("nuclei", "spin_counts"):
            if key in data:
                data[key] = tuple(data[key])
        return cls(**data)

    @classmethod
    def load(cls, path) -> "ProblemSpec":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    def save(self, path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
