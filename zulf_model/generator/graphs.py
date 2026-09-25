"""Procedural molecule-like graphs obeying valence rules.

Graphs contain explicit hydrogens. Heavy atoms carry an element, and after
generation a hybridization label derived from bond orders. Isotope labels
(13C, 15N) are attached later by `isotopologues`. Graphs are plausible
skeletons for sampling spin systems, not validated chemical structures.
"""
from __future__ import annotations

import hashlib
import json
from collections import deque
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

VALENCE = {"C": 4, "N": 3, "O": 2, "H": 1, "F": 1, "P": 3}


@dataclass(frozen=True)
class GraphConfig:
    heavy_atoms: Tuple[int, int] = (2, 7)
    element_weights: Dict[str, float] = field(default_factory=lambda: {"C": 0.75, "N": 0.15, "O": 0.10})
    ring_closure_probability: float = 0.25
    aromatic_ring_probability: float = 0.15
    aromatic_nitrogen_probability: float = 0.2
    double_bond_probability: float = 0.2
    triple_bond_probability: float = 0.03
    max_attempts: int = 50

    @classmethod
    def from_dict(cls, data: dict) -> "GraphConfig":
        data = {k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items() if k in cls.__dataclass_fields__}
        return cls(**data)


@dataclass
class MoleculeGraph:
    elements: List[str]
    bonds: Dict[Tuple[int, int], int]
    aromatic: List[bool]
    labels: Dict[int, str] = field(default_factory=dict)  # atom index -> isotope symbol

    # -- structure --------------------------------------------------------
    def n_atoms(self) -> int:
        return len(self.elements)

    def neighbors(self, i: int) -> List[int]:
        return [b if a == i else a for (a, b) in self.bonds if i in (a, b)]

    def bond_order(self, i: int, j: int) -> int:
        return self.bonds.get((min(i, j), max(i, j)), 0)

    def used_valence(self, i: int) -> int:
        return sum(o for (a, b), o in self.bonds.items() if i in (a, b))

    def free_valence(self, i: int) -> int:
        return VALENCE[self.elements[i]] - self.used_valence(i)

    def add_bond(self, i: int, j: int, order: int = 1) -> None:
        key = (min(i, j), max(i, j))
        self.bonds[key] = self.bonds.get(key, 0) + order

    def heavy_indices(self) -> List[int]:
        return [i for i, e in enumerate(self.elements) if e != "H"]

    def hydrogens_on(self, i: int) -> List[int]:
        return [k for k in self.neighbors(i) if self.elements[k] == "H"]

    def parent(self, h: int) -> int:
        return self.neighbors(h)[0]

    def hybridization(self, i: int) -> str:
        if self.elements[i] == "H":
            return "s"
        if self.aromatic[i]:
            return "ar"
        orders = [self.bond_order(i, k) for k in self.neighbors(i)]
        if 3 in orders or orders.count(2) >= 2:
            return "sp"
        if 2 in orders:
            return "sp2"
        return "sp3"

    def distances(self) -> np.ndarray:
        n = self.n_atoms()
        adjacency = [self.neighbors(i) for i in range(n)]
        out = np.full((n, n), 10**6, dtype=int)
        for s in range(n):
            out[s, s] = 0
            queue = deque([s])
            while queue:
                u = queue.popleft()
                for v in adjacency[u]:
                    if out[s, v] > out[s, u] + 1:
                        out[s, v] = out[s, u] + 1
                        queue.append(v)
        return out

    def shortest_path(self, i: int, j: int) -> List[int]:
        previous = {i: None}
        queue = deque([i])
        while queue:
            u = queue.popleft()
            if u == j:
                break
            for v in self.neighbors(u):
                if v not in previous:
                    previous[v] = u
                    queue.append(v)
        path = [j]
        while previous[path[-1]] is not None:
            path.append(previous[path[-1]])
        return path[::-1]

    def in_ring(self, i: int, j: int) -> bool:
        """True when bond i-j lies on a cycle."""
        key = (min(i, j), max(i, j))
        order = self.bonds.pop(key)
        try:
            seen = {i}
            queue = deque([i])
            while queue:
                u = queue.popleft()
                for v in self.neighbors(u):
                    if v not in seen:
                        seen.add(v)
                        queue.append(v)
            return j in seen
        finally:
            self.bonds[key] = order

    def validate(self) -> None:
        for i in range(self.n_atoms()):
            if self.free_valence(i) != 0:
                raise ValueError(f"Atom {i} ({self.elements[i]}) has unsatisfied valence.")
        if len(self.distances()[0][self.distances()[0] < 10**6]) != self.n_atoms():
            raise ValueError("Graph is disconnected.")

    # -- symmetry ---------------------------------------------------------
    def wl_colors(self, rounds: Optional[int] = None) -> List[int]:
        """Weisfeiler-Lehman color classes including isotope labels and bond orders."""
        n = self.n_atoms()
        colors = [hash((self.elements[i], self.labels.get(i, ""), self.hybridization(i), self.aromatic[i]))
                  for i in range(n)]
        colors = _compress(colors)
        for _ in range(rounds or n):
            new = [hash((colors[i], tuple(sorted((colors[k], self.bond_order(i, k)) for k in self.neighbors(i)))))
                   for i in range(n)]
            new = _compress(new)
            if len(set(new)) == len(set(colors)):
                colors = new
                break
            colors = new
        return colors

    def skeleton_hash(self) -> str:
        """Label-free topology identifier used for family splits."""
        saved = self.labels
        self.labels = {}
        try:
            colors = self.wl_colors()
            summary = sorted((self.elements[i], colors[i]) for i in range(self.n_atoms()))
            multiset = sorted(colors)
        finally:
            self.labels = saved
        text = json.dumps([summary, multiset, len(self.bonds)])
        return hashlib.sha256(text.encode()).hexdigest()[:16]

    def to_dict(self) -> dict:
        return {"elements": self.elements, "bonds": [[a, b, o] for (a, b), o in sorted(self.bonds.items())],
                "aromatic": self.aromatic, "labels": {str(k): v for k, v in self.labels.items()}}

    @classmethod
    def from_dict(cls, data: dict) -> "MoleculeGraph":
        return cls(list(data["elements"]), {(a, b): o for a, b, o in data["bonds"]}, list(data["aromatic"]),
                   {int(k): v for k, v in data.get("labels", {}).items()})

    def copy(self) -> "MoleculeGraph":
        return MoleculeGraph(list(self.elements), dict(self.bonds), list(self.aromatic), dict(self.labels))


def _compress(values: Sequence[int]) -> List[int]:
    ranks = {v: r for r, v in enumerate(sorted(set(values)))}
    return [ranks[v] for v in values]


def _choose(rng, weights: Dict[str, float]) -> str:
    keys = list(weights)
    p = np.array([weights[k] for k in keys], float)
    return str(rng.choice(keys, p=p / p.sum()))


def _add_aromatic_ring(graph: MoleculeGraph, rng, config: GraphConfig) -> List[int]:
    start = len(graph.elements)
    for k in range(6):
        element = "N" if (k > 0 and rng.random() < config.aromatic_nitrogen_probability / 6 * 3) else "C"
        graph.elements.append(element)
        graph.aromatic.append(True)
    ring = list(range(start, start + 6))
    for k in range(6):
        # Kekule form keeps valence bookkeeping exact.
        graph.add_bond(ring[k], ring[(k + 1) % 6], 2 if k % 2 == 0 else 1)
    return ring


def random_graph(rng: np.random.Generator, config: GraphConfig) -> MoleculeGraph:
    """Draw a connected, valence-satisfied graph with explicit hydrogens."""
    for _ in range(config.max_attempts):
        graph = MoleculeGraph([], {}, [])
        target = int(rng.integers(config.heavy_atoms[0], config.heavy_atoms[1] + 1))
        if target >= 6 and rng.random() < config.aromatic_ring_probability:
            _add_aromatic_ring(graph, rng, config)
        if not graph.elements:
            graph.elements.append("C")
            graph.aromatic.append(False)
        while len(graph.elements) < target:
            candidates = [i for i in graph.heavy_indices() if graph.free_valence(i) > 0]
            if not candidates:
                break
            anchor = int(rng.choice(candidates))
            element = _choose(rng, config.element_weights)
            graph.elements.append(element)
            graph.aromatic.append(False)
            graph.add_bond(anchor, len(graph.elements) - 1, 1)
        heavy = graph.heavy_indices()
        distances = graph.distances()
        if rng.random() < config.ring_closure_probability:
            pairs = [(i, j) for i in heavy for j in heavy if i < j and 2 <= distances[i, j] <= 6
                     and graph.free_valence(i) > 0 and graph.free_valence(j) > 0
                     and not (graph.aromatic[i] and graph.aromatic[j])]
            if pairs:
                i, j = pairs[int(rng.integers(len(pairs)))]
                graph.add_bond(i, j, 1)
        for (i, j), order in list(graph.bonds.items()):
            if graph.aromatic[i] or graph.aromatic[j] or order != 1:
                continue
            if graph.free_valence(i) >= 1 and graph.free_valence(j) >= 1:
                u = rng.random()
                if u < config.triple_bond_probability and graph.free_valence(i) >= 2 and graph.free_valence(j) >= 2 \
                        and not graph.in_ring(i, j) and {graph.elements[i], graph.elements[j]} <= {"C", "N"}:
                    graph.add_bond(i, j, 2)
                elif u < config.triple_bond_probability + config.double_bond_probability:
                    graph.add_bond(i, j, 1)
        for i in list(graph.heavy_indices()):
            for _ in range(graph.free_valence(i)):
                graph.elements.append("H")
                graph.aromatic.append(False)
                graph.add_bond(i, len(graph.elements) - 1, 1)
        try:
            graph.validate()
        except ValueError:
            continue
        if not any(e == "C" for e in graph.elements):
            continue
        return graph
    raise RuntimeError("Failed to generate a valid graph within max_attempts.")
