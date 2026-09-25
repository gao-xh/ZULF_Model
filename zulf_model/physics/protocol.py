"""Observation protocol: preparation, optional ideal pulses, detection.

v1 uses `sudden_drop`: gamma-weighted longitudinal preparation, sudden removal
of the field, no pulse, gamma-weighted detection. See docs/CONVENTIONS.md.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Tuple, Union

from ..nuclei import get_registry

Weights = Union[str, Dict[str, float]]


@dataclass(frozen=True)
class Pulse:
    """Ideal instantaneous rotation of all spins of the listed nuclei."""
    nuclei: Tuple[str, ...]
    angle_rad: float
    axis: str = "x"

    def __post_init__(self):
        if self.axis not in ("x", "y"):
            raise ValueError("Pulse axis must be 'x' or 'y'.")
        object.__setattr__(self, "nuclei", tuple(self.nuclei))


@dataclass(frozen=True)
class Protocol:
    name: str = "sudden_drop"
    preparation: Weights = "gamma"
    detection: Weights = "gamma"
    pulses: Tuple[Pulse, ...] = field(default_factory=tuple)
    normalize_by_dimension: bool = True
    time_origin_s: float = 0.0

    def __post_init__(self):
        object.__setattr__(self, "pulses", tuple(self.pulses))
        for weights in (self.preparation, self.detection):
            if isinstance(weights, str) and weights != "gamma":
                raise ValueError("Weights must be 'gamma' or a mapping nucleus -> weight.")

    @staticmethod
    def _weight(weights: Weights, symbol: str) -> float:
        if weights == "gamma":
            return get_registry().gamma(symbol)
        return float(weights.get(symbol, 0.0))

    def preparation_weight(self, symbol: str) -> float:
        return self._weight(self.preparation, symbol)

    def detection_weight(self, symbol: str) -> float:
        return self._weight(self.detection, symbol)

    @property
    def is_real(self) -> bool:
        """True when every operator is real symmetric (no pulses)."""
        return not self.pulses

    def to_dict(self) -> dict:
        data = asdict(self)
        data["pulses"] = [asdict(p) for p in self.pulses]
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "Protocol":
        data = dict(data)
        data.pop("frozen", None)
        data.pop("notes", None)
        data["pulses"] = tuple(Pulse(**p) for p in data.get("pulses", ()))
        return cls(**data)

    @classmethod
    def load(cls, path) -> "Protocol":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


SUDDEN_DROP = Protocol()
