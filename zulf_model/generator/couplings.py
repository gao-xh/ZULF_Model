"""Rule-based J assignment for spin atoms of a labelled molecule graph.

Rules are keyed by bond-path length and the nuclei involved; numeric ranges
come from configuration (`configs/couplings_v1.json`) and are approximate
literature-style ranges to be verified (open question Q3). After sampling, J
values are averaged over pairs that share Weisfeiler-Lehman color classes, so
rotors and symmetry-related atoms become magnetically equivalent.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .graphs import MoleculeGraph

DEFAULT_RULES: Dict[str, dict] = {
    "1J_CH": {"sp3": [120.0, 135.0], "sp2": [150.0, 170.0], "ar": [155.0, 168.0], "sp": [245.0, 255.0],
              "electronegative_shift": [3.0, 8.0]},
    "1J_NH": {"range": [-95.0, -70.0]},
    "1J_CC": {"sp3-sp3": [32.0, 40.0], "sp2": [50.0, 75.0], "ar": [52.0, 60.0], "sp": [60.0, 80.0]},
    "1J_CN": {"range": [-15.0, -4.0]},
    "2J_HH": {"sp3": [-16.0, -8.0], "sp2": [-3.0, 3.0]},
    "2J_CH": {"range": [-6.0, 6.0]},
    "2J_NH": {"range": [-5.0, 5.0]},
    "2J_heavy": {"range": [-3.0, 3.0]},
    "3J_HH": {"karplus": [7.76, -1.10, 1.40], "cis": [6.0, 12.0], "trans": [12.0, 18.0],
              "ortho": [6.5, 8.5]},
    "3J_CH": {"karplus": [4.5, -0.9, 0.9], "range": [0.0, 10.0]},
    "3J_NH": {"range": [-5.0, 5.0]},
    "3J_heavy": {"range": [-2.0, 2.0]},
    "4J_HH_meta": {"range": [1.0, 3.0]},
    "long": {"laplace_scale": 0.4, "clip": 3.0},
    "para": {"range": [0.0, 1.0]},
}


@dataclass
class CouplingRules:
    rules: Dict[str, dict] = field(default_factory=lambda: json.loads(json.dumps(DEFAULT_RULES)))
    rotor_averaging: bool = True
    long_range_zero_probability: float = 0.5
    perturbation_relative: float = 0.0

    @classmethod
    def load(cls, path) -> "CouplingRules":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        rules = json.loads(json.dumps(DEFAULT_RULES))
        rules.update(data.get("rules", {}))
        return cls(rules, data.get("rotor_averaging", True), data.get("long_range_zero_probability", 0.5),
                   data.get("perturbation_relative", 0.0))

    # ------------------------------------------------------------------
    def _uniform(self, rng, key, sub="range") -> float:
        lo, hi = self.rules[key][sub]
        return float(rng.uniform(lo, hi))

    @staticmethod
    def _karplus(coeffs, phi) -> float:
        a, b, c = coeffs
        return a * np.cos(phi) ** 2 + b * np.cos(phi) + c

    def _vicinal(self, rng, graph: MoleculeGraph, path: List[int], key: str) -> float:
        """3J across the central bond path[1]-path[2]."""
        b, c = path[1], path[2]
        order = graph.bond_order(b, c)
        rule = self.rules[key]
        if key == "3J_HH":
            if graph.aromatic[b] and graph.aromatic[c]:
                return float(rng.uniform(*rule["ortho"]))
            if order == 2:
                return float(rng.uniform(*(rule["cis"] if rng.random() < 0.5 else rule["trans"])))
        coeffs = rule.get("karplus")
        if coeffs is None:
            return self._uniform(rng, key)
        if graph.in_ring(b, c):
            return self._karplus(coeffs, rng.uniform(0, np.pi))
        # Free rotation: population-weighted average over staggered rotamers.
        phase = rng.uniform(0, 2 * np.pi / 3)
        populations = rng.dirichlet(np.ones(3) * 3)
        return float(sum(p * self._karplus(coeffs, phase + k * 2 * np.pi / 3) for k, p in enumerate(populations)))

    def pair_coupling(self, rng, graph: MoleculeGraph, i: int, j: int, distance: int) -> float:
        el = {i: graph.elements[i], j: graph.elements[j]}
        nuclei = tuple(sorted((el[i], el[j])))
        if distance == 1:
            if nuclei == ("C", "H"):
                c = i if el[i] == "C" else j
                hyb = graph.hybridization(c)
                value = float(rng.uniform(*self.rules["1J_CH"][hyb]))
                electronegative = sum(1 for k in graph.neighbors(c) if graph.elements[k] in ("N", "O", "F"))
                return value + electronegative * float(rng.uniform(*self.rules["1J_CH"]["electronegative_shift"]))
            if nuclei == ("H", "N"):
                return self._uniform(rng, "1J_NH")
            if nuclei == ("C", "C"):
                hi, hj = graph.hybridization(i), graph.hybridization(j)
                key = "ar" if "ar" in (hi, hj) else "sp" if "sp" in (hi, hj) else "sp2" if "sp2" in (hi, hj) else "sp3-sp3"
                return float(rng.uniform(*self.rules["1J_CC"][key]))
            if nuclei == ("C", "N"):
                return self._uniform(rng, "1J_CN")
            return 0.0
        path = graph.shortest_path(i, j)
        if distance == 2:
            center = path[1]
            if nuclei == ("H", "H"):
                hyb = graph.hybridization(center)
                return float(rng.uniform(*self.rules["2J_HH"]["sp3" if hyb == "sp3" else "sp2"]))
            if nuclei == ("C", "H"):
                return self._uniform(rng, "2J_CH")
            if nuclei == ("H", "N"):
                return self._uniform(rng, "2J_NH")
            return self._uniform(rng, "2J_heavy")
        if distance == 3:
            if nuclei == ("H", "H"):
                return self._vicinal(rng, graph, path, "3J_HH")
            if nuclei == ("C", "H"):
                return self._vicinal(rng, graph, path, "3J_CH")
            if nuclei == ("H", "N"):
                return self._uniform(rng, "3J_NH")
            return self._uniform(rng, "3J_heavy")
        if distance == 4 and nuclei == ("H", "H") and all(graph.aromatic[k] for k in path[1:-1]):
            return self._uniform(rng, "4J_HH_meta")
        if distance == 5 and nuclei == ("H", "H") and all(graph.aromatic[k] for k in path[1:-1]):
            return self._uniform(rng, "para")
        if rng.random() < self.long_range_zero_probability or distance > 5:
            return 0.0
        scale, clip = self.rules["long"]["laplace_scale"], self.rules["long"]["clip"]
        return float(np.clip(rng.laplace(0, scale), -clip, clip))

    def assign(self, rng: np.random.Generator, graph: MoleculeGraph, spin_atoms: Sequence[int]) -> np.ndarray:
        """Symmetric J matrix (Hz) over `spin_atoms`, averaged over WL pair classes."""
        atoms = list(spin_atoms)
        distances = graph.distances()
        n = len(atoms)
        j = np.zeros((n, n))
        for a in range(n):
            for b in range(a + 1, n):
                value = self.pair_coupling(rng, graph, atoms[a], atoms[b], int(distances[atoms[a], atoms[b]]))
                if self.perturbation_relative:
                    value *= 1 + rng.normal(0, self.perturbation_relative)
                j[a, b] = j[b, a] = value
        if self.rotor_averaging:
            colors = graph.wl_colors()
            classes: Dict[Tuple, List[Tuple[int, int]]] = {}
            for a in range(n):
                for b in range(a + 1, n):
                    ca, cb = colors[atoms[a]], colors[atoms[b]]
                    key = (min(ca, cb), max(ca, cb), int(distances[atoms[a], atoms[b]]))
                    classes.setdefault(key, []).append((a, b))
            for pairs in classes.values():
                mean = float(np.mean([j[a, b] for a, b in pairs]))
                for a, b in pairs:
                    j[a, b] = j[b, a] = mean
        return j
