"""Candidate models (PyTorch)."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..spec import ProblemSpec
from .base import CandidateModel
from .config import EncoderConfig, ModelConfig, SequenceHeadConfig, SetHeadConfig

MODEL_KINDS = {}


def register_model(kind: str):
    def decorator(cls):
        MODEL_KINDS[kind] = cls
        return cls
    return decorator


def build_model(spec: ProblemSpec, config: ModelConfig) -> CandidateModel:
    from .set_model import CNNSetModel
    from .transformer import CNNTransformerModel
    MODEL_KINDS.setdefault("cnn_set", CNNSetModel)
    MODEL_KINDS.setdefault("cnn_transformer", CNNTransformerModel)
    if config.kind not in MODEL_KINDS:
        raise ValueError(f"Unknown model kind '{config.kind}'. Known: {sorted(MODEL_KINDS)}")
    return MODEL_KINDS[config.kind](spec, config)


def save_model(path, model: CandidateModel, extra: Optional[dict] = None) -> None:
    """Device-agnostic checkpoint with spec and config for reconstruction."""
    import torch
    state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
    torch.save({"state_dict": state, "spec": model.spec.to_dict(), "model_config": model.model_config.to_dict(),
                "extra": extra or {}}, Path(path))


def load_model(path, device: str = "cpu") -> CandidateModel:
    import torch
    blob = torch.load(Path(path), map_location="cpu", weights_only=False)
    model = build_model(ProblemSpec.from_dict(blob["spec"]), ModelConfig.from_dict(blob["model_config"]))
    model.load_state_dict(blob["state_dict"])
    return model.to(device)


__all__ = ["CandidateModel", "ModelConfig", "EncoderConfig", "SetHeadConfig", "SequenceHeadConfig", "build_model",
           "save_model", "load_model", "register_model", "MODEL_KINDS"]
