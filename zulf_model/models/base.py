"""Common interface for candidate-generating models.

Any model that implements `loss` and `propose` can be trained by
`training.Trainer`, wrapped by `evaluation.ModelProposer` and exposed through
the agent tools.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

import torch
from torch import nn

from ..codec import InterpretationCodec
from ..spec import ProblemSpec
from zulf_core.spinsystem import Interpretation
from .config import ModelConfig


class CandidateModel(nn.Module, ABC):
    kind: str = "base"

    def __init__(self, spec: ProblemSpec, config: ModelConfig):
        super().__init__()
        self.spec = spec
        self.model_config = config
        self.codec = InterpretationCodec(spec)

    @abstractmethod
    def forward(self, batch: Dict[str, torch.Tensor]) -> Dict[str, torch.Tensor]:
        """Training forward pass on a collated batch."""

    @abstractmethod
    def loss(self, outputs: Dict[str, torch.Tensor], batch: Dict[str, torch.Tensor]) -> Tuple[torch.Tensor, Dict[str, float]]:
        """Scalar loss and a dictionary of loggable floats."""

    @abstractmethod
    def propose(self, features: torch.Tensor, frequency_hz: torch.Tensor, k: int = 5) -> List[List[Interpretation]]:
        """Top-k interpretations per spectrum, best first. Must not raise on unusual inputs."""

    def parameter_count(self) -> int:
        return sum(p.numel() for p in self.parameters())
