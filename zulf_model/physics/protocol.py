"""Observation protocol: preparation, optional ideal pulses, detection.

v1 uses `sudden_drop`: gamma-weighted longitudinal preparation, sudden removal
of the field, no pulse, gamma-weighted detection. See docs/CONVENTIONS.md.

`field_ut` is an optional static field (Bx, By, Bz) in microtesla present
during free evolution, for example the residual field inside a shield. It adds
the Zeeman term -sum gamma_n B . I_n (Hz) to the J Hamiltonian; the default
(0, 0, 0) is exact zero field.
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
    field_ut: Tuple[float, float, float] = (0.0, 0.0, 0.0)

    def __post_init__(self):
        object.__setattr__(self, "pulses", tuple(self.pulses))
        field_ut = tuple(float(v) for v in self.field_ut)
        if len(field_ut) != 3:
            raise ValueError("field_ut must have three components (Bx, By, Bz).")
        object.__setattr__(self, "field_ut", field_ut)
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
        """True when every operator is real symmetric (no pulses, no y field)."""
        return not self.pulses and self.field_ut[1] == 0.0

    @property
    def has_field(self) -> bool:
        return any(v != 0.0 for v in self.field_ut)

    def with_field(self, field_ut) -> "Protocol":
        data = self.to_dict()
        data["field_ut"] = tuple(float(v) for v in field_ut)
        return Protocol.from_dict(data)

    def to_dict(self) -> dict:
        data = asdict(self)
        data["pulses"] = [asdict(p) for p in self.pulses]
        data["field_ut"] = list(self.field_ut)
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "Protocol":
        data = dict(data)
        data.pop("frozen", None)
        data.pop("notes", None)
        data["pulses"] = tuple(Pulse(**p) for p in data.get("pulses", ()))
        if "field_ut" in data:
            data["field_ut"] = tuple(data["field_ut"])
        return cls(**data)

    @classmethod
    def load(cls, path) -> "Protocol":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


SUDDEN_DROP = Protocol()
