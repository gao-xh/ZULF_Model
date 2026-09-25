"""Sample sources, mixtures, deterministic family splits and generator versioning.

A `Sample` is a ground-truth `Interpretation` plus provenance. Splits are
decided by hashing the family identifier (the label-free molecular topology,
or a random-system identifier), so every render of one underlying system and
every isotopologue of one molecule land in the same partition.
"""
from __future__ import annotations

import hashlib
import json
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence

import numpy as np

from ..spec import ProblemSpec
from ..spinsystem import Component, Interpretation
from .couplings import CouplingRules
from .graphs import GraphConfig, MoleculeGraph, random_graph
from .isotopologues import IsotopologueConfig, labeled_isotopologue, natural_isotopologues, visible_hydrogens
from .random_j import RandomSystemConfig, random_system

GENERATOR_VERSION = "gen-v1"
SPLITS = ("train", "val", "test_unseen", "test_seen")


@dataclass
class Sample:
    interpretation: Interpretation
    family_id: str
    source: str
    seed: List[int]
    generator_version: str = GENERATOR_VERSION
    graph: Optional[dict] = None
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"interpretation": self.interpretation.to_dict(), "family_id": self.family_id,
                "source": self.source, "seed": list(self.seed), "generator_version": self.generator_version,
                "graph": self.graph, "metadata": self.metadata}

    @classmethod
    def from_dict(cls, data: dict) -> "Sample":
        return cls(Interpretation.from_dict(data["interpretation"]), data["family_id"], data["source"],
                   list(data["seed"]), data.get("generator_version", GENERATOR_VERSION), data.get("graph"),
                   dict(data.get("metadata", {})))


@dataclass(frozen=True)
class SplitConfig:
    val_fraction: float = 0.05
    test_unseen_fraction: float = 0.05
    test_seen_fraction: float = 0.02
    salt: str = "gen-v1"

    @classmethod
    def from_dict(cls, data: dict) -> "SplitConfig":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


def _unit_hash(text: str) -> float:
    return int(hashlib.sha256(text.encode()).hexdigest()[:15], 16) / 16 ** 15


def split_of(family_id: str, sample_key: str, config: SplitConfig) -> str:
    """Deterministic partition. Families decide val/test_unseen; within training
    families a small fraction of individual systems forms test_seen."""
    u = _unit_hash(config.salt + "|family|" + family_id)
    if u < config.test_unseen_fraction:
        return "test_unseen"
    if u < config.test_unseen_fraction + config.val_fraction:
        return "val"
    if _unit_hash(config.salt + "|sample|" + sample_key) < config.test_seen_fraction:
        return "test_seen"
    return "train"


class SampleSource(ABC):
    name: str = "source"

    @abstractmethod
    def draw(self, rng: np.random.Generator) -> Optional[Sample]:
        """Return a sample or None when this attempt was rejected."""


class MoleculeSource(SampleSource):
    name = "molecule"

    def __init__(self, spec: ProblemSpec, graph_config: GraphConfig = GraphConfig(),
                 rules: Optional[CouplingRules] = None, iso_config: IsotopologueConfig = IsotopologueConfig()):
        self.spec = spec
        self.graph_config = graph_config
        self.rules = rules or CouplingRules()
        self.iso_config = iso_config

    def draw(self, rng: np.random.Generator) -> Optional[Sample]:
        target = int(rng.choice(self.spec.spin_counts))
        graph = random_graph(rng, self.graph_config)
        n_h = len(visible_hydrogens(graph, self.iso_config.exchangeable_protons))
        if rng.random() < self.iso_config.natural_probability:
            if n_h + 1 != target:
                return None
            components = natural_isotopologues(rng, graph, self.rules, self.spec, self.iso_config)
            mode = "natural"
        else:
            components = labeled_isotopologue(rng, graph, self.rules, self.spec, self.iso_config, target)
            mode = "labeled"
        if not components or len(components) > self.spec.max_components:
            return None
        return Sample(Interpretation(tuple(components), metadata={"mode": mode}),
                      "mol:" + graph.skeleton_hash(), self.name, [], graph=graph.to_dict(),
                      metadata={"mode": mode})


class RandomJSource(SampleSource):
    name = "random_j"

    def __init__(self, spec: ProblemSpec, config: RandomSystemConfig = RandomSystemConfig(),
                 max_components: Optional[int] = None):
        self.spec = spec
        self.config = config
        self.max_components = max_components or spec.max_components

    def draw(self, rng: np.random.Generator) -> Optional[Sample]:
        count = int(rng.integers(1, self.max_components + 1))
        components = []
        for _ in range(count):
            system = random_system(rng, self.spec, self.config)
            components.append(Component(system, float(10 ** rng.uniform(*self.config.contribution_log10_range)),
                                        "random"))
        key = hashlib.sha256(json.dumps([c.system.to_dict() for c in components]).encode()).hexdigest()[:16]
        return Sample(Interpretation(tuple(components), metadata={"mode": "random_j"}), "rnd:" + key,
                      self.name, [], metadata={"mode": "random_j"})


@dataclass(frozen=True)
class MixtureConfig:
    source_weights: Dict[str, float] = field(default_factory=lambda: {"molecule": 0.85, "random_j": 0.15})
    extra_molecule_probability: float = 0.15
    extra_molecule_log10_ratio: tuple = (-1.0, 1.0)
    max_attempts: int = 2000

    @classmethod
    def from_dict(cls, data: dict) -> "MixtureConfig":
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items()
                      if k in cls.__dataclass_fields__})


class MixtureSampler:
    """Draws samples from weighted sources, optionally combining molecules.

    All families of a combined sample must fall in the requested split.
    """

    def __init__(self, spec: ProblemSpec, sources: Sequence[SampleSource], config: MixtureConfig = MixtureConfig(),
                 split_config: SplitConfig = SplitConfig()):
        self.spec = spec
        self.sources = {s.name: s for s in sources}
        missing = set(config.source_weights) - set(self.sources)
        if missing:
            raise ValueError(f"No source registered for {sorted(missing)}")
        self.config = config
        self.split_config = split_config
        self.attempts = 0
        self.accepted = 0

    def _draw_source(self, rng, name: str, split: Optional[str]) -> Optional[Sample]:
        sample = self.sources[name].draw(rng)
        if sample is None:
            return None
        if split is not None:
            key = sample.family_id + "|" + json.dumps(sample.interpretation.to_dict(), sort_keys=True)[:4096]
            if split_of(sample.family_id, hashlib.sha256(key.encode()).hexdigest(), self.split_config) != split:
                return None
        return sample

    def draw(self, rng: np.random.Generator, split: Optional[str] = None) -> Sample:
        names = list(self.config.source_weights)
        p = np.array([self.config.source_weights[n] for n in names], float)
        p /= p.sum()
        for _ in range(self.config.max_attempts):
            self.attempts += 1
            name = str(rng.choice(names, p=p))
            sample = self._draw_source(rng, name, split)
            if sample is None:
                continue
            if name == "molecule" and rng.random() < self.config.extra_molecule_probability:
                other = None
                for _ in range(50):
                    other = self._draw_source(rng, "molecule", split)
                    if other is not None:
                        break
                if other is not None:
                    ratio = 10 ** rng.uniform(*self.config.extra_molecule_log10_ratio)
                    merged = list(sample.interpretation.components) + [
                        Component(c.system, c.contribution * ratio, c.label, c.metadata)
                        for c in other.interpretation.components]
                    if len(merged) <= self.spec.max_components:
                        sample = Sample(Interpretation(tuple(merged), metadata={"mode": "mixture"}),
                                        sample.family_id + "+" + other.family_id, "mixture", [],
                                        metadata={"families": [sample.family_id, other.family_id],
                                                  "graphs": [sample.graph, other.graph]})
            self.accepted += 1
            return sample
        raise RuntimeError("Sampler exceeded max_attempts; relax the spec or generator settings.")

    def generate(self, count: int, seed: int, split: Optional[str] = None) -> Iterator[Sample]:
        """Reproducible stream: sample k uses rng seeded by (seed, k)."""
        for k in range(count):
            rng = np.random.default_rng([seed, k])
            sample = self.draw(rng, split)
            sample.seed = [seed, k]
            yield sample

    @property
    def acceptance_rate(self) -> float:
        return self.accepted / self.attempts if self.attempts else 0.0


def build_default_sampler(spec: ProblemSpec, config: Optional[dict] = None) -> MixtureSampler:
    """Construct the standard sampler from a generator configuration dictionary."""
    config = config or {}
    rules = CouplingRules.load(config["couplings_path"]) if config.get("couplings_path") else CouplingRules()
    sources = [MoleculeSource(spec, GraphConfig.from_dict(config.get("graph", {})), rules,
                              IsotopologueConfig.from_dict(config.get("isotopologues", {}))),
               RandomJSource(spec, RandomSystemConfig.from_dict(config.get("random_j", {})))]
    return MixtureSampler(spec, sources, MixtureConfig.from_dict(config.get("mixture", {})),
                          SplitConfig.from_dict(config.get("splits", {})))


def load_generator_config(path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))
