"""Fragment -> refinable natural-abundance isotopologue set.

One component per symmetry-distinct labelled site and isotope (13C, 15N, ...):
the labelled nucleus plus every included proton group. Proton groups that
the labelled site's stabiliser maps onto each other are merged into one
magnetically equivalent group (isopropyl labelled at the CH: six methyl
protons). Couplings related by symmetry, and every proton-proton coupling
shared between isotopologues, become one tied parameter. Amplitude ratios
are natural abundance times the number of equivalent sites.

Replaces hand-built interpretations and hand-mapped ties (group indices
differ between isotopologues). The result feeds `RefineSettings` directly.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from ..nuclei import get_registry
from ..physics import compute_transitions
from ..physics.protocol import SUDDEN_DROP, Protocol
from ..spinsystem import Component, Interpretation, SpinSystem
from .fragment import Fragment, Pair, orbits_under, pair


@dataclass
class HypothesisModel:
    """A built hypothesis: interpretation plus the bookkeeping the solver and the report need."""
    name: str
    fragment: Fragment
    interpretation: Interpretation
    component_labels: List[str]
    ratios: Tuple[float, ...]                 # natural-abundance amplitude ratios (first component = 1)
    coupling_names: Dict[str, List[str]]      # key name 'J(a,b)' -> solver parameter names
    unspecified: Tuple[str, ...]              # key names with no value in the fragment (built as 0 Hz)
    omitted: List[Tuple[str, str]] = field(default_factory=list)   # (component label, reason)

    def ties(self, shared_rate: bool = False, rate_families: int = 1) -> Tuple[Tuple[str, ...], ...]:
        out = [tuple(names) for names in self.coupling_names.values() if len(names) > 1]
        if shared_rate and len(self.component_labels) > 1:
            for f in range(rate_families):
                out.append(tuple(f"c{c}.log_rate{f}" for c in range(len(self.component_labels))))
        return tuple(out)

    def fixed(self) -> Tuple[str, ...]:
        return tuple(n for key in self.unspecified for n in self.coupling_names[key])

    def settings(self, base=None, fixed_ratios: bool = False, shared_rate: bool = False,
                 fix_unspecified: bool = True):
        """RefineSettings with this model's ties, fixed couplings and (optionally) abundance ratios."""
        from ..solver import RefineSettings
        base = base or RefineSettings()
        families = len(base.policy.family_edges_hz) + 1
        ties = tuple(base.ties) + self.ties(shared_rate, families)
        fixed = tuple(base.fixed) + (self.fixed() if fix_unspecified else ())
        return dataclasses.replace(base, ties=ties, fixed=fixed,
                                   amplitude_ratios=self.ratios if fixed_ratios else base.amplitude_ratios)

    def named_couplings(self, parameters: Mapping[str, float]) -> Dict[str, float]:
        """Refined couplings by key name (first instance of each tied parameter)."""
        return {key: float(parameters[names[0]]) for key, names in self.coupling_names.items()
                if names and names[0] in parameters}

    def summary(self) -> dict:
        return {"name": self.name, "components": self.component_labels, "ratios": list(self.ratios),
                "couplings": {k: v for k, v in self.coupling_names.items()}, "unspecified": list(self.unspecified),
                "omitted": [list(o) for o in self.omitted], "fragment": self.fragment.to_dict()}


def build_model(fragment: Fragment, include_exchangeable: bool = True, ranges: Optional[Sequence[Tuple[float, float]]] = None,
                protocol: Protocol = SUDDEN_DROP, min_relative_amplitude: float = 1e-3,
                name: Optional[str] = None) -> HypothesisModel:
    """Build the natural-abundance isotopologue set of `fragment`.

    `include_exchangeable=False` drops exchangeable protons (fast exchange).
    With `ranges`, components without a line in any range are omitted and
    listed with the reason; so are labelled sites without any coupling to an
    included proton group.
    """
    registry = get_registry()
    elements = fragment.group_elements()
    orbits = fragment.coupling_orbits(elements)
    protons = [p for p in fragment.protons if include_exchangeable or not p.exchangeable]
    proton_labels = [p.label for p in protons]
    site_labels = [s.label for s in fragment.sites]

    components, labels, ratios, omitted = [], [], [], []
    names: Dict[str, List[str]] = {}
    unspecified = set()
    for site_orbit in orbits_under(elements, site_labels):
        rep = site_orbit[0]
        site = fragment.site(rep)
        stabiliser = [g for g in elements if g[rep] == rep]
        merged = orbits_under(stabiliser, proton_labels)
        for isotope in site.label_isotopes():
            label = f"{isotope}@{rep}" + (f" (x{len(site_orbit)})" if len(site_orbit) > 1 else "")
            if all(fragment.coupling(rep, g[0], orbits) is None for g in merged):
                omitted.append((label, "no coupling to an included proton group is given"))
                continue
            spins = [(rep,)] + [tuple(g) for g in merged]
            sizes = [1] + [sum(fragment.proton(x).size for x in g) for g in merged]
            n = len(spins)
            matrix = np.zeros((n, n))
            keys: Dict[Tuple[int, int], Pair] = {}
            for a in range(n):
                for b in range(a + 1, n):
                    reps = {orbits[pair(x, y)] for x in spins[a] for y in spins[b]}
                    if len(reps) != 1:
                        raise ValueError(f"{label}: groups {spins[a]} and {spins[b]} are merged as equivalent but "
                                         "their couplings differ; check the symmetry.")
                    key = reps.pop()
                    value = fragment.coupling(spins[a][0], spins[b][0], orbits)
                    if value is None:
                        unspecified.add(fragment.key_name(key))
                        value = 0.0
                    matrix[a, b] = matrix[b, a] = value
                    keys[(a, b)] = key
            system = SpinSystem.from_group_couplings([isotope] + ["1H"] * (n - 1), sizes, matrix)
            if ranges is not None:
                tl = compute_transitions(system, protocol)
                amp = np.abs(tl.amplitudes)
                strong = amp >= min_relative_amplitude * (amp.max() if len(amp) else 1.0)
                inside = np.zeros(len(tl), bool)
                for lo, hi in ranges:
                    inside |= (tl.frequencies_hz >= lo) & (tl.frequencies_hz <= hi)
                if not np.any(strong & inside):
                    omitted.append((label, "no line in the fitted ranges"))
                    continue
            c = len(components)
            for (a, b), key in keys.items():
                names.setdefault(fragment.key_name(key), []).append(f"c{c}.J{a}-{b}")
            ratio = registry[isotope].natural_abundance * len(site_orbit)
            components.append(system)
            labels.append(label)
            ratios.append(ratio)
    if not components:
        raise ValueError(f"Fragment {fragment.name!r} gives no observable isotopologue.")
    top = ratios[0]
    ratios = tuple(r / top for r in ratios)
    interp = Interpretation(tuple(Component(s, r, l) for s, r, l in zip(components, ratios, labels)))
    # Keys that no retained component uses are dropped (for example couplings of omitted sites).
    names = {k: v for k, v in names.items() if v}
    return HypothesisModel(name or fragment.name, fragment, interp, labels, ratios, names,
                           tuple(sorted(k for k in unspecified if k in names)), omitted)
