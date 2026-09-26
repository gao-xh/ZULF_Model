"""Active-learning rounds: evaluate on development data, mine failures, fine-tune.

Each round records the generator version, the failure histogram and metrics,
so successive rounds are comparable. The loop never touches frozen test data.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence

import numpy as np

from ..generator.sampler import GENERATOR_VERSION, MixtureConfig, MixtureSampler, Sample
from zulf_core.solver.refine import RefineSettings, refine_candidates
from .failures import FailureRecord, classify_failure
from .focused import FocusedSource


@dataclass(frozen=True)
class LoopConfig:
    rounds: int = 3
    k: int = 5
    focused_fraction: float = 0.3
    steps_per_round: int = 1000
    learning_rate_scale: float = 0.3
    jitter_hz: float = 1.0
    refine: bool = True
    tolerance_hz: float = 0.1
    output_dir: str = "runs/active_learning"


class ActiveLearningLoop:
    """Orchestrates rounds around injectable components.

    `observe(sample) -> (train_observation, held_out_list)` renders a development
    sample; `proposer` returns candidates; `fine_tune(sampler, steps, lr_scale)`
    trains the model on a sampler mixing base and focused sources.
    """

    def __init__(self, dev_samples: Sequence[Sample], proposer, observe: Callable, fine_tune: Callable,
                 base_sampler: MixtureSampler, config: LoopConfig = LoopConfig(),
                 refine_settings: RefineSettings = RefineSettings(starts=1)):
        self.dev = list(dev_samples)
        self.proposer, self.observe, self.fine_tune = proposer, observe, fine_tune
        self.base = base_sampler
        self.config = config
        self.refine_settings = refine_settings
        self.output = Path(config.output_dir)
        self.output.mkdir(parents=True, exist_ok=True)

    def evaluate(self) -> List[FailureRecord]:
        records = []
        for sample in self.dev:
            train, held = self.observe(sample)
            candidates = self.proposer.propose(train, self.config.k)
            refined = refine_candidates(candidates, train, held, self.refine_settings) if self.config.refine else None
            records.append(classify_failure(sample, candidates, refined, self.config.tolerance_hz))
        return records

    def mixed_sampler(self, failures: Sequence[FailureRecord]) -> MixtureSampler:
        focused = FocusedSource(self.base.spec, failures, self.config.jitter_hz, self.base.split_config)
        sources = list(self.base.sources.values()) + [focused]
        weights = dict(self.base.config.source_weights)
        total = sum(weights.values())
        weights = {k: v / total * (1 - self.config.focused_fraction) for k, v in weights.items()}
        weights["focused"] = self.config.focused_fraction if focused.failures else 0.0
        config = MixtureConfig(**{**asdict(self.base.config), "source_weights": weights})
        return MixtureSampler(self.base.spec, sources, config, self.base.split_config)

    def run(self) -> List[dict]:
        history = []
        for round_index in range(self.config.rounds):
            records = self.evaluate()
            histogram = Counter(r.category for r in records)
            entry = {"round": round_index, "generator_version": GENERATOR_VERSION,
                     "failures": dict(histogram), "success_rate": histogram.get("success", 0) / max(1, len(records))}
            history.append(entry)
            with (self.output / "rounds.jsonl").open("a", encoding="utf-8") as handle:
                handle.write(json.dumps(entry) + "\n")
            if round_index == self.config.rounds - 1:
                break
            sampler = self.mixed_sampler(records)
            self.fine_tune(sampler, self.config.steps_per_round, self.config.learning_rate_scale)
        return history
