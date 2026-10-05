"""Paired data for J -> structure: a generator graph and the J network a ZULF fit would observe for it.

`observation_from_graph` turns a `MoleculeGraph` into the observation format of zulf_hypothesis.j_structure
(units = carbons that carry protons, grouped into symmetry classes by Weisfeiler-Lehman colours; one proton group
per protonated carbon; couplings from the 13C of copy 0 of every unit to every proton group, and between proton
groups on different carbons, from `CouplingRules` with Gaussian noise) and the true structure as a bond list over
those carbons and unseen atoms X (heteroatoms and carbons without protons). Graphs that route A cannot represent
(more unseen atoms than allowed, a bond between two unseen atoms, an unseen atom with more than three bonds) are
rejected (None), and so are graphs with exchangeable protons on carbon-free parts only. The observation also
carries "bond_counts" ({coupling: bonds in the true graph}) and "hybridization" ({unit: sp3 | sp2 | sp of its copy-0
carbon, aromatic counted as sp2}), the labels route B learns from.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from ..generator.couplings import CouplingRules
from ..generator.graphs import MoleculeGraph


def observation_from_graph(rng: np.random.Generator, graph: MoleculeGraph, rules: CouplingRules,
                           sigma: float = 0.1, max_unseen: int = 2) -> Optional[Tuple[dict, str]]:
    heavy = graph.heavy_indices()
    carbons = [i for i in heavy if graph.elements[i] == "C" and graph.hydrogens_on(i)]
    if not carbons:
        return None
    unseen = [i for i in heavy if i not in carbons]
    if len(unseen) > max_unseen:
        return None
    for x in unseen:
        nb = [k for k in graph.neighbors(x) if graph.elements[k] != "H"]
        if any(k in unseen for k in nb) or len(nb) > 3 or not nb:
            return None
    colors = graph.wl_colors()
    classes: Dict[int, List[int]] = {}
    for c in carbons:
        classes.setdefault(colors[c], []).append(c)
    units, protons, name_of = {}, {}, {}
    for k, (col, members) in enumerate(sorted(classes.items(), key=lambda kv: min(kv[1]))):
        label = f"C{k + 1}"
        units[label] = {"h": len(graph.hydrogens_on(members[0])), "copies": len(members)}
        for copy, atom in enumerate(sorted(members)):
            name_of[atom] = label + "'" * copy
            protons[f"H{label[1:]}{chr(97 + copy)}"] = [label, copy]
    group_atoms = {g: graph.hydrogens_on(sorted(classes_by_label(classes, units, g, protons))[protons[g][1]])
                   for g in protons}
    # J over the spin atoms: copy-0 carbon of every unit and all protons on carbon
    spin_carbons = [sorted(members)[0] for _, members in sorted(classes.items(), key=lambda kv: min(kv[1]))]
    spin_h = [h for g in protons for h in group_atoms[g]]
    atoms = spin_carbons + spin_h
    j = rules.assign(rng, graph, atoms)
    idx = {a: i for i, a in enumerate(atoms)}
    couplings, bond_counts = {}, {}
    dist = graph.distances()
    labels = list(units)
    for u, c0 in zip(labels, spin_carbons):
        for g, hs in group_atoms.items():
            couplings[f"J({u},{g})"] = float(np.mean([j[idx[c0], idx[h]] for h in hs]) + rng.normal(0, sigma))
            bond_counts[f"J({u},{g})"] = int(dist[c0, graph.parent(hs[0])]) + 1
    names = list(protons)
    for a in range(len(names)):
        for b in range(a + 1, len(names)):
            ha, hb = group_atoms[names[a]], group_atoms[names[b]]
            couplings[f"J({names[a]},{names[b]})"] = float(np.mean([j[idx[x], idx[y]] for x in ha for y in hb])
                                                          + rng.normal(0, sigma))
            bond_counts[f"J({names[a]},{names[b]})"] = int(dist[graph.parent(ha[0]), graph.parent(hb[0])]) + 2
    xname = {x: f"X{k + 1}" for k, x in enumerate(unseen)}
    bonds = []
    for (a, b), order in sorted(graph.bonds.items()):
        if graph.elements[a] == "H" or graph.elements[b] == "H":
            continue
        na = name_of.get(a, xname.get(a))
        nb = name_of.get(b, xname.get(b))
        bonds.append(f"{na}{'-=#'[min(order, 3) - 1]}{nb}")
    hyb = {u: {"ar": "sp2"}.get(graph.hybridization(c0), graph.hybridization(c0)) for u, c0 in zip(labels, spin_carbons)}
    return ({"units": units, "protons": protons, "couplings": couplings, "bond_counts": bond_counts,
             "hybridization": hyb},
            ", ".join(bonds))


def classes_by_label(classes, units, group, protons):
    label = protons[group][0]
    k = list(units).index(label)
    return sorted(classes.items(), key=lambda kv: min(kv[1]))[k][1]
