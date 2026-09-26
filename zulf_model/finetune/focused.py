"""Sample source concentrated around development failures.

For each failure it re-draws couplings near the failing system (same spin
graph, J jittered by `jitter_hz`) or, for structure misses, re-samples
isotopologues of the same molecule graph when available. Families from frozen
test splits are refused.
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Set

import numpy as np

from zulf_core.evaluation.identifiability import perturb_couplings
from ..generator.couplings import CouplingRules
from ..generator.graphs import MoleculeGraph
from ..generator.isotopologues import IsotopologueConfig, natural_isotopologues
from ..generator.sampler import Sample, SampleSource, SplitConfig, split_of
from ..spec import ProblemSpec
from zulf_core.spinsystem import Component, Interpretation
from .failures import FailureRecord

FROZEN_SPLITS = ("test_unseen", "test_seen")


class FocusedSource(SampleSource):
    name = "focused"

    def __init__(self, spec: ProblemSpec, failures: Sequence[FailureRecord], jitter_hz: float = 1.0,
                 split_config: SplitConfig = SplitConfig(), rules: Optional[CouplingRules] = None,
                 iso_config: IsotopologueConfig = IsotopologueConfig()):
        self.spec = spec
        self.jitter_hz = jitter_hz
        self.rules = rules or CouplingRules()
        self.iso_config = iso_config
        self.failures = []
        for f in failures:
            if f.category in ("success", "unidentifiable"):
                continue
            if split_of(f.sample.family_id, "", split_config) in FROZEN_SPLITS:
                raise ValueError("Refusing to build focused data from a frozen test family.")
            self.failures.append(f)

    def draw(self, rng: np.random.Generator) -> Optional[Sample]:
        if not self.failures:
            return None
        failure = self.failures[int(rng.integers(len(self.failures)))]
        sample = failure.sample
        if failure.category == "structure_miss" and sample.graph and "elements" in sample.graph:
            graph = MoleculeGraph.from_dict(sample.graph)
            comps = natural_isotopologues(rng, graph, self.rules, self.spec, self.iso_config)
            if comps and len(comps) <= self.spec.max_components:
                return Sample(Interpretation(tuple(comps)), sample.family_id, self.name, [], graph=sample.graph,
                              metadata={"focused_on": failure.category})
        comps = tuple(Component(perturb_couplings(c.system, self.jitter_hz, rng), c.contribution, c.label)
                      for c in sample.interpretation.components)
        return Sample(Interpretation(comps), sample.family_id, self.name, [], graph=sample.graph,
                      metadata={"focused_on": failure.category})
