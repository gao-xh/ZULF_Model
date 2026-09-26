"""Candidate proposers sharing one interface.

A proposer turns an observed spectrum into ranked interpretations. The model
proposer wraps a trained network; the baselines use the same generator prior
without learning, so comparisons at equal simulation budget are fair.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional

import numpy as np

from ..generator.sampler import MixtureSampler
from zulf_core.physics.protocol import SUDDEN_DROP, Protocol
from zulf_core.physics.transitions import TransitionCache
from zulf_core.render.acquisition import evaluate_spectrum, process_record
from ..render.features import spectrum_features
from zulf_core.render.grid import SpectrumGrid
from zulf_core.render.phasing import estimate_phase, phase_correct
from zulf_core.render.renderer import ContinuousRenderer, Renderer
from zulf_core.solver.observed import ObservedSpectrum
from ..spec import ProblemSpec
from zulf_core.spinsystem import Interpretation


class CandidateProposer(ABC):
    name = "proposer"

    @abstractmethod
    def propose(self, observed: ObservedSpectrum, k: int) -> List[Interpretation]:
        """Return up to k interpretations, best first."""


class ModelProposer(CandidateProposer):
    """Wraps a trained CandidateModel. The observed spectrum is re-evaluated on the
    model grid through its acquisition when a stored FID is available.

    Models trained with `grid.phasing = "corrected"` need the operator's phase
    correction: `phasing = {"phase0_rad": ..., "delay_s": ...}` (see render.phasing);
    the crop reference of the acquisition is added automatically. `{"auto": True}`
    estimates both from the spectrum (`render.phasing.estimate_phase`).
    """
    name = "model"

    def __init__(self, model, grid: SpectrumGrid, channels=("real", "imag"), device: str = "cpu",
                 phasing: Optional[dict] = None, **propose_kwargs):
        self.model = model.to(device).eval()
        self.grid = grid
        self.channels = channels
        self.device = device
        self.phasing = phasing
        self.kwargs = propose_kwargs
        if getattr(model.spec.grid, "phasing", "none") == "corrected" and not phasing:
            raise ValueError("This model expects phase-corrected spectra; pass phasing={'phase0_rad', 'delay_s'}.")

    def features(self, observed: ObservedSpectrum) -> np.ndarray:
        if observed.reprocessable:
            acq = observed.acquisition
            values = evaluate_spectrum(process_record(observed.fid, acq), acq, self.grid.frequencies_hz)
        else:
            values = np.interp(self.grid.frequencies_hz, observed.frequencies_hz, observed.values.real) + 1j * np.interp(
                self.grid.frequencies_hz, observed.frequencies_hz, observed.values.imag)
        if self.phasing:
            phasing = dict(self.phasing)
            if phasing.get("auto"):
                estimate = estimate_phase(values, self.grid.frequencies_hz, observed.acquisition,
                                          **{k: v for k, v in phasing.items() if k != "auto"})
                phasing = {"phase0_rad": estimate["phase0_rad"], "delay_s": estimate["delay_s"]}
                self.last_phasing = estimate
            values = phase_correct(values, self.grid.frequencies_hz, float(phasing.get("phase0_rad", 0.0)),
                                   float(phasing.get("delay_s", 0.0)), observed.acquisition)
        return spectrum_features(values, self.channels)[0]

    def propose(self, observed: ObservedSpectrum, k: int) -> List[Interpretation]:
        import torch
        x = torch.as_tensor(self.features(observed)).unsqueeze(0).to(self.device)
        freq = torch.as_tensor(self.grid.frequencies_hz, dtype=torch.float32, device=self.device)
        return self.model.propose(x, freq, k=k, **self.kwargs)[0]


def _render_candidate(interp: Interpretation, observed: ObservedSpectrum, cache: TransitionCache,
                      protocol: Protocol, rate: float) -> np.ndarray:
    renderer = Renderer(observed.acquisition) if observed.acquisition is not None else ContinuousRenderer()
    f = observed.frequencies_hz[observed.selected]
    cols = []
    for c in interp.components:
        cols.append(renderer.render_pair(cache.get(c.system, protocol), rate, f))
    return np.hstack(cols)


def quick_score(interp: Interpretation, observed: ObservedSpectrum, cache: TransitionCache,
                protocol: Protocol = SUDDEN_DROP, rate: float = 1.0) -> float:
    """Relative residual after solving only linear complex gains (no nonlinear refinement)."""
    design = _render_candidate(interp, observed, cache, protocol, rate)
    y = observed.values[observed.selected]
    a = np.vstack([design.real, design.imag])
    b = np.r_[y.real, y.imag]
    coef = np.linalg.lstsq(a, b, rcond=1e-10)[0]
    return float(np.linalg.norm(a @ coef - b) / max(np.linalg.norm(b), 1e-30))


class RandomPriorProposer(CandidateProposer):
    """Random multistart baseline: k independent draws from the generator prior."""
    name = "random_prior"

    def __init__(self, sampler: MixtureSampler, seed: int = 0):
        self.sampler = sampler
        self.rng = np.random.default_rng(seed)

    def propose(self, observed: ObservedSpectrum, k: int) -> List[Interpretation]:
        return [self.sampler.draw(self.rng).interpretation for _ in range(k)]


class PriorSearchProposer(CandidateProposer):
    """Non-learned search: draw `pool` prior samples, rank by gain-only fit, keep top k.

    Its simulation budget is `pool` forward evaluations per spectrum.
    """
    name = "prior_search"

    def __init__(self, sampler: MixtureSampler, pool: int = 200, seed: int = 0, protocol: Protocol = SUDDEN_DROP,
                 rate: float = 1.0):
        self.sampler, self.pool, self.protocol, self.rate = sampler, pool, protocol, rate
        self.rng = np.random.default_rng(seed)
        self.cache = TransitionCache(4096)

    def propose(self, observed: ObservedSpectrum, k: int) -> List[Interpretation]:
        scored = []
        for _ in range(self.pool):
            interp = self.sampler.draw(self.rng).interpretation
            try:
                scored.append((quick_score(interp, observed, self.cache, self.protocol, self.rate), interp))
            except ValueError:
                continue
        scored.sort(key=lambda s: s[0])
        return [Interpretation(i.components, -score, {"source": self.name}) for score, i in scored[:k]]
