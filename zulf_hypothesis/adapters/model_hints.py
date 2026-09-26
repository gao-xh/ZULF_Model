"""Neural models (zulf_model checkpoints) as hint providers for `zulf_hypothesis`.

`model_hints(checkpoint, acquisition)` loads a trained checkpoint and wraps
its `zulf_model` `ModelProposer` in `ProposerHints`: the model's top-k interpretations
become group hints and extra candidates. They steer the hypothesis search
and appear in its insight report; they never decide the ranking.
"""
from __future__ import annotations

from typing import Optional

from zulf_core.render.grid import SpectrumGrid

from ..hints import ProposerHints, register_hint_provider


def model_hints(checkpoint, acquisition, phasing: Optional[dict] = None, k: int = 5, name: Optional[str] = None,
                register: bool = False) -> ProposerHints:
    """Hint provider from a checkpoint; `phasing` as for ModelProposer (default: automatic)."""
    # zulf_model (torch) is imported only here, so the hypothesis package runs without it.
    from zulf_model.evaluation.proposers import ModelProposer
    from zulf_model.models import load_model
    model = load_model(checkpoint).eval()
    grid = SpectrumGrid.from_spec(model.spec.grid, acquisition)
    proposer = ModelProposer(model, grid, model.spec.grid.channels, phasing=phasing or {"auto": True})
    provider = ProposerHints(proposer, name=name or f"model:{checkpoint}", k=k)
    if register:
        register_hint_provider(provider)
    return provider
