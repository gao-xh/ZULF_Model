"""Chemistry-free random spin systems (range-constrained J only).

Used as a fraction of training data and as a separate evaluation subset that
measures how much the model relies on the chemical rules.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from ..spec import ProblemSpec
from zulf_core.spinsystem import Component, SpinSystem
from .isotopologues import has_observable_heteronuclear_coupling


@dataclass(frozen=True)
class RandomSystemConfig:
    heteronuclei: Tuple[int, int] = (1, 2)
    group_probability: float = 0.5
    strong_heteronuclear_probability: float = 0.5
    strong_range_hz: Tuple[float, float] = (60.0, 250.0)
    medium_range_hz: Tuple[float, float] = (-20.0, 20.0)
    small_laplace_hz: float = 1.0
    negative_strong_probability: float = 0.1
    contribution_log10_range: Tuple[float, float] = (-2.0, 0.0)

    @classmethod
    def from_dict(cls, data: dict) -> "RandomSystemConfig":
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items()
                      if k in cls.__dataclass_fields__})


def _partition(rng, count: int, max_size: int, probability: float) -> List[int]:
    sizes, left = [], count
    while left:
        size = 1
        if left > 1 and rng.random() < probability:
            size = int(rng.integers(2, min(max_size, left) + 1))
        sizes.append(size)
        left -= size
    return sizes


def random_system(rng: np.random.Generator, spec: ProblemSpec, config: RandomSystemConfig,
                  n_spins: Optional[int] = None, max_attempts: int = 50) -> SpinSystem:
    hetero_symbols = [s for s in spec.nuclei if s != "1H"]
    if not hetero_symbols or "1H" not in spec.nuclei:
        raise ValueError("Random systems need 1H and at least one heteronucleus in the spec.")
    for _ in range(max_attempts):
        n = int(rng.choice(spec.spin_counts)) if n_spins is None else n_spins
        n_het = int(rng.integers(config.heteronuclei[0], min(config.heteronuclei[1], n - 1) + 1))
        hetero = [str(rng.choice(hetero_symbols)) for _ in range(n_het)]
        h_sizes = _partition(rng, n - n_het, spec.max_group_size, config.group_probability)
        group_iso = hetero + ["1H"] * len(h_sizes)
        sizes = [1] * n_het + h_sizes
        g = len(sizes)
        gj = np.zeros((g, g))
        for a in range(g):
            for b in range(a + 1, g):
                if group_iso[a] != group_iso[b] and rng.random() < config.strong_heteronuclear_probability:
                    value = rng.uniform(*config.strong_range_hz)
                    if rng.random() < config.negative_strong_probability:
                        value = -value
                elif rng.random() < 0.5:
                    value = rng.uniform(*config.medium_range_hz)
                else:
                    value = rng.laplace(0, config.small_laplace_hz)
                gj[a, b] = gj[b, a] = value
        system = SpinSystem.from_group_couplings(group_iso, sizes, gj)
        if spec.allows(system.isotopes) and has_observable_heteronuclear_coupling(system):
            return system
    raise RuntimeError("Failed to draw a valid random system.")


def random_components(rng: np.random.Generator, spec: ProblemSpec, config: RandomSystemConfig,
                      count: int) -> List[Component]:
    return [Component(random_system(rng, spec, config), float(10 ** rng.uniform(*config.contribution_log10_range)),
                      "random") for _ in range(count)]
