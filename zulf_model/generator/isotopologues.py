"""Isotopologue enumeration and spin-system construction from molecule graphs.

J values are sampled once per molecule over every atom that can carry a spin
(hydrogens plus all C and N), symmetrized over the label-free graph symmetry.
Each isotopologue takes the sub-matrix of its spin atoms, so isotopologues of
one molecule share H-H couplings and remain consistent with each other.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from zulf_core.nuclei import get_registry
from ..spec import ProblemSpec
from zulf_core.spinsystem import Component, SpinSystem
from .couplings import CouplingRules
from .graphs import MoleculeGraph

LABEL_OF_ELEMENT = {"C": "13C", "N": "15N", "F": "19F", "P": "31P"}


@dataclass(frozen=True)
class IsotopologueConfig:
    natural_probability: float = 0.7
    include_natural_n15: bool = True
    labeled_max_labels: int = 3
    exchangeable_protons: str = "drop"  # "drop" or "keep"
    min_relative_contribution: float = 0.0

    @classmethod
    def from_dict(cls, data: dict) -> "IsotopologueConfig":
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class MoleculeCouplings:
    graph: MoleculeGraph
    atoms: List[int]           # graph atom indices that may carry spins
    couplings_hz: np.ndarray   # J over `atoms` (heavy atoms as if labelled)

    def submatrix(self, spin_atoms: Sequence[int]) -> np.ndarray:
        index = [self.atoms.index(a) for a in spin_atoms]
        return self.couplings_hz[np.ix_(index, index)]


def visible_hydrogens(graph: MoleculeGraph, exchangeable: str) -> List[int]:
    out = []
    for i, e in enumerate(graph.elements):
        if e != "H":
            continue
        if exchangeable == "drop" and graph.elements[graph.parent(i)] in ("N", "O"):
            continue
        out.append(i)
    return out


def molecule_couplings(rng: np.random.Generator, graph: MoleculeGraph, rules: CouplingRules,
                       exchangeable: str, label_elements: Sequence[str]) -> MoleculeCouplings:
    unlabelled = graph.copy()
    unlabelled.labels = {}
    atoms = visible_hydrogens(unlabelled, exchangeable)
    atoms += [i for i, e in enumerate(unlabelled.elements) if e in label_elements]
    return MoleculeCouplings(unlabelled, atoms, rules.assign(rng, unlabelled, atoms))


def has_observable_heteronuclear_coupling(system: SpinSystem) -> bool:
    j = system.observable_couplings()
    iso = np.array(system.isotopes, dtype=object)
    return bool(np.any((iso[:, None] != iso[None, :]) & (j != 0)))


def build_system(mc: MoleculeCouplings, labels: Dict[int, str], exchangeable: str,
                 spec: ProblemSpec) -> Optional[SpinSystem]:
    """Spin system of one isotopologue, or None when outside the problem spec."""
    spin_atoms = visible_hydrogens(mc.graph, exchangeable) + sorted(labels)
    isotopes = tuple("1H" if mc.graph.elements[a] == "H" else labels[a] for a in spin_atoms)
    if not spec.allows(isotopes):
        return None
    system = SpinSystem(isotopes, mc.submatrix(spin_atoms))
    if max(len(g) for g in system.groups) > spec.max_group_size:
        return None
    if not has_observable_heteronuclear_coupling(system):
        return None
    return system


def label_classes(graph: MoleculeGraph, element: str) -> List[List[int]]:
    """Symmetry classes (label-free WL colors) of atoms of one element."""
    saved = graph.labels
    graph.labels = {}
    try:
        colors = graph.wl_colors()
    finally:
        graph.labels = saved
    classes: Dict[int, List[int]] = {}
    for i, e in enumerate(graph.elements):
        if e == element:
            classes.setdefault(colors[i], []).append(i)
    return list(classes.values())


def natural_isotopologues(rng: np.random.Generator, graph: MoleculeGraph, rules: CouplingRules,
                          spec: ProblemSpec, config: IsotopologueConfig) -> List[Component]:
    """Singly labelled isotopologues at natural abundance (one component per symmetry class)."""
    registry = get_registry()
    elements = [e for e, sym in LABEL_OF_ELEMENT.items() if sym in spec.nuclei]
    if "N" in elements and not config.include_natural_n15:
        elements.remove("N")
    mc = molecule_couplings(rng, graph, rules, config.exchangeable_protons, elements)
    counts = {e: sum(1 for x in graph.elements if x == e) for e in elements}
    components = []
    for element in elements:
        symbol = LABEL_OF_ELEMENT[element]
        p = registry[symbol].natural_abundance
        for cls in label_classes(graph, element):
            representative = cls[0]
            system = build_system(mc, {representative: symbol}, config.exchangeable_protons, spec)
            if system is None:
                continue
            others = 1.0
            for e2 in elements:
                q = registry[LABEL_OF_ELEMENT[e2]].natural_abundance
                others *= (1 - q) ** (counts[e2] - (1 if e2 == element else 0))
            contribution = len(cls) * p * others
            components.append(Component(system, contribution, f"{symbol}@{representative}",
                                        {"labels": {str(representative): symbol}, "class_size": len(cls)}))
    if components and config.min_relative_contribution > 0:
        top = max(c.contribution for c in components)
        components = [c for c in components if c.contribution >= config.min_relative_contribution * top]
    return components


def labeled_isotopologue(rng: np.random.Generator, graph: MoleculeGraph, rules: CouplingRules,
                         spec: ProblemSpec, config: IsotopologueConfig,
                         target_spins: Optional[int] = None) -> List[Component]:
    """One selectively labelled isotopologue with 1..labeled_max_labels labels."""
    elements = [e for e, sym in LABEL_OF_ELEMENT.items() if sym in spec.nuclei]
    candidates = [i for i, e in enumerate(graph.elements) if e in elements]
    if not candidates:
        return []
    n_h = len(visible_hydrogens(graph, config.exchangeable_protons))
    if target_spins is not None:
        count = target_spins - n_h
    else:
        count = int(rng.integers(1, config.labeled_max_labels + 1))
    if not 1 <= count <= min(config.labeled_max_labels, len(candidates)):
        return []
    chosen = sorted(rng.choice(candidates, size=count, replace=False).tolist())
    labels = {i: LABEL_OF_ELEMENT[graph.elements[i]] for i in chosen}
    mc = molecule_couplings(rng, graph, rules, config.exchangeable_protons, elements)
    system = build_system(mc, labels, config.exchangeable_protons, spec)
    if system is None:
        return []
    return [Component(system, 1.0, "labeled", {"labels": {str(k): v for k, v in labels.items()}})]
