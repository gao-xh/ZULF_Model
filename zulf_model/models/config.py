"""Model hyperparameters, loaded from JSON. No problem dimensions live here;
those come from `ProblemSpec`."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List, Tuple


@dataclass(frozen=True)
class EncoderConfig:
    in_channels: int = 2
    d_model: int = 128
    stages: Tuple[Tuple[int, int, int], ...] = ((32, 9, 2), (64, 9, 2), (96, 7, 2), (128, 7, 2), (128, 5, 2))
    dropout: float = 0.1
    frequency_features: int = 32
    frequency_min_period_hz: float = 0.5
    frequency_max_period_hz: float = 2000.0
    transformer_layers: int = 2
    heads: int = 4
    feedforward: int = 256
    max_tokens: int = 1024


@dataclass(frozen=True)
class SetHeadConfig:
    decoder_layers: int = 3
    heads: int = 4
    feedforward: int = 256
    j_scale_hz: float = 100.0
    weight_exists: float = 1.0
    weight_isotope: float = 1.0
    weight_j: float = 1.0
    weight_contribution: float = 0.2


@dataclass(frozen=True)
class SequenceHeadConfig:
    decoder_layers: int = 3
    heads: int = 4
    feedforward: int = 256
    label_smoothing: float = 0.05
    j_neighbor_sigma_bins: float = 1.0
    weight_offset: float = 0.5


@dataclass(frozen=True)
class ModelConfig:
    kind: str = "cnn_set"  # "cnn_set" or "cnn_transformer"
    encoder: EncoderConfig = field(default_factory=EncoderConfig)
    set_head: SetHeadConfig = field(default_factory=SetHeadConfig)
    sequence_head: SequenceHeadConfig = field(default_factory=SequenceHeadConfig)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "ModelConfig":
        enc = dict(data.get("encoder", {}))
        if "stages" in enc:
            enc["stages"] = tuple(tuple(s) for s in enc["stages"])
        return cls(kind=data.get("kind", "cnn_set"), encoder=EncoderConfig(**enc),
                   set_head=SetHeadConfig(**data.get("set_head", {})),
                   sequence_head=SequenceHeadConfig(**data.get("sequence_head", {})))

    @classmethod
    def load(cls, path) -> "ModelConfig":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
