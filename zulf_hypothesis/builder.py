"""Fragment -> refinable isotopologue set (natural abundance or enriched, see labeling.py).

One component per symmetry-distinct set of labelled sites and isotopes
(13C, 15N, ...; single labels by default, pairs with `max_labels=2`): the
labelled nuclei plus every included proton group. Groups that the label
set's stabiliser maps onto each other are merged into one magnetically
equivalent group (isopropyl labelled at the CH: six methyl protons; both
methyl carbons labelled: one pair of equivalent 13C). Couplings related by
symmetry, and every coupling shared between isotopologues, become one tied
parameter. Amplitude ratios are the exact label-set probabilities of the
`Labeling` times the number of equivalent label sets.

Minor isotopologues (below `primary_fraction` of the strongest set: doubly
labelled ones at natural abundance, partly labelled ones in an enriched
sample) carry their abundance as a fixed weight relative to a parent
component: they never add a free amplitude or a free decay rate
(`amplitude_map`, rate ties), so they cannot absorb background; and label
sets whose relative amplitude is below `min_ratio` are omitted with the
reason (at natural abundance a 13C-13C isotopologue is about 1 % of its
parent, below the noise at SNR < 100).

Replaces hand-built interpretations and hand-mapped ties (group indices
differ between isotopologues). The result feeds `RefineSettings` directly.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_core.physics.protocol import SUDDEN_DROP, Protocol
from zulf_core.spinsystem import Component, Interpretation, SpinSystem
from .fragment import Fragment, Pair, orbits_under, pair
from .labeling import Labeling


@dataclass
class HypothesisModel:
    """A built hypothesis: interpretation plus the bookkeeping the solver and the report need."""
    name: str
    fragment: Fragment
    interpretation: Interpretation
    component_labels: List[str]
    ratios: Tuple[float, ...]                 # label-set abundance ratios (first primary component = 1)
    coupling_names: Dict[str, List[str]]      # key name 'J(a,b)' -> solver parameter names
    unspecified: Tuple[str, ...]              # key names with no value in the fragment (built as 0 Hz)
    omitted: List[Tuple[str, str]] = field(default_factory=list)   # (component label, reason)
    ratio_blocks: Optional[List[List[int]]] = None   # components whose ratios are fixed together (None: all)
    parts: List["HypothesisModel"] = field(default_factory=list)   # separate molecules of a combined model
    parents: List[Optional[int]] = field(default_factory=list)    # minor isotopologue -> its parent component

    def parent(self, c: int) -> Optional[int]:
        return self.parents[c] if c < len(self.parents) else None

    def primaries(self) -> List[int]:
        return [c for c in range(len(self.component_labels)) if self.parent(c) is None]

    def amplitude_map(self, fixed_ratios: bool) -> Optional[Tuple[Tuple[float, ...], ...]]:
        """Components x free amplitudes. Free: one column per primary component, minor isotopologues at their
        abundance factor on the parent's column. Fixed: one column per ratio block. None when nothing is tied."""
        n = len(self.component_labels)
        if fixed_ratios:
            cols = self.blocks()
            m = np.zeros((n, len(cols)))
            for j, block in enumerate(cols):
                for c in block:
                    m[c, j] = self.ratios[c]
        else:
            prim = self.primaries()
            if len(prim) == n:
                return None
            m = np.zeros((n, len(prim)))
            for j, c in enumerate(prim):
                m[c, j] = 1.0
            for c in range(n):
                p = self.parent(c)
                if p is not None:
                    m[c, prim.index(p)] = self.ratios[c] / self.ratios[p]
        return tuple(tuple(float(x) for x in row) for row in m)

    def blocks(self) -> List[List[int]]:
        return self.ratio_blocks or [list(range(len(self.component_labels)))]

    def ties(self, shared_rate: bool = False, rate_families: int = 1) -> Tuple[Tuple[str, ...], ...]:
        out = [tuple(names) for names in self.coupling_names.values() if len(names) > 1]
        if shared_rate:
            for block in self.blocks():
                if len(block) > 1:
                    for f in range(rate_families):
                        out.append(tuple(f"c{c}.log_rate{f}" for c in block))
        else:
            # Minor isotopologues decay with their parent (no free rate of their own).
            children: Dict[int, List[int]] = {}
            for c in range(len(self.component_labels)):
                if self.parent(c) is not None:
                    children.setdefault(self.parent(c), []).append(c)
            for p, cs in children.items():
                for f in range(rate_families):
                    out.append(tuple(f"c{c}.log_rate{f}" for c in [p] + cs))
        return tuple(out)

    def fixed(self) -> Tuple[str, ...]:
        return tuple(n for key in self.unspecified for n in self.coupling_names[key])

    def settings(self, base=None, fixed_ratios: bool = False, shared_rate: bool = False,
                 fix_unspecified: bool = True):
        """RefineSettings with this model's ties, fixed couplings and (optionally) abundance ratios."""
        from zulf_core.solver import RefineSettings
        base = base or RefineSettings()
        families = len(base.policy.family_edges_hz) + 1
        ties = tuple(base.ties) + self.ties(shared_rate, families)
        fixed = tuple(base.fixed) + (self.fixed() if fix_unspecified else ())
        # Fixed ratios hold within each part (separate molecules keep free relative amounts); minor
        # isotopologues always follow their parent.
        return dataclasses.replace(base, ties=ties, fixed=fixed, amplitude_ratios=None,
                                   amplitude_map=self.amplitude_map(fixed_ratios))

    def named_couplings(self, parameters: Mapping[str, float]) -> Dict[str, float]:
        """Refined couplings by key name (first instance of each tied parameter)."""
        return {key: float(parameters[names[0]]) for key, names in self.coupling_names.items()
                if names and names[0] in parameters}

    def summary(self) -> dict:
        return {"name": self.name, "components": self.component_labels, "ratios": list(self.ratios),
                "ratio_blocks": self.blocks(), "parents": list(self.parents),
                "couplings": {k: v for k, v in self.coupling_names.items()}, "unspecified": list(self.unspecified),
                "omitted": [list(o) for o in self.omitted], "fragment": self.fragment.to_dict()}


def build_model(fragment: Fragment, include_exchangeable: bool = True, ranges: Optional[Sequence[Tuple[float, float]]] = None,
                protocol: Protocol = SUDDEN_DROP, min_relative_amplitude: float = 1e-3,
                name: Optional[str] = None, labeling: Optional[Labeling] = None, max_labels: Optional[int] = None,
                min_ratio: float = 0.0) -> HypothesisModel:
    """Build the isotopologue set of `fragment`.

    `include_exchangeable=False` drops exchangeable protons (fast exchange).
    `labeling`: natural abundance (default, single labels) or enrichment
    (`Labeling.enriched(...)`); `max_labels` overrides its label-set size.
    `min_ratio`: omit label sets whose amplitude relative to the strongest
    set is below this (see `min_ratio_for_snr`).
    With `ranges`, components without a line in any range are omitted; so are
    label sets without any coupling to an included proton group. Every
    omission is listed with its reason.
    """
    labeling = labeling or Labeling.natural()
    elements = fragment.group_elements()
    orbits = fragment.coupling_orbits(elements)
    protons = [p for p in fragment.protons if include_exchangeable or not p.exchangeable]
    proton_labels = [p.label for p in protons]
    labelable = [s.label for s in fragment.sites if s.label_isotopes()]
    limit = max_labels if max_labels is not None else labeling.label_limit(len(labelable))
    for g in elements:          # site-specific enrichment must respect the fragment symmetry
        for x in labelable:
            for iso in fragment.site(x).label_isotopes():
                if abs(labeling.abundance(x, iso) - labeling.abundance(g[x], iso)) > 1e-12:
                    raise ValueError(f"Labelling of {x} and {g[x]} differs but the fragment makes them equivalent; "
                                     "remove that symmetry for this sample.")
    unlabelled = {x: 1.0 - sum(labeling.abundance(x, iso) for iso in fragment.site(x).label_isotopes())
                  for x in labelable}

    # Enumerate symmetry-distinct label sets with their isotopes.
    entries = []   # (label groups [(isotope, (sites...))], equivalent count, raw ratio, representative sites)
    for k in range(1, limit + 1):
        seen = set()
        for combo in combinations(labelable, k):
            orbit = {frozenset(g[x] for x in combo) for g in elements}
            key = min(orbit, key=lambda q: sorted(fragment.order(x) for x in q))
            if key in seen:
                continue
            seen.add(key)
            rep = tuple(sorted(key, key=fragment.order))
            stab = [g for g in elements if frozenset(g[x] for x in rep) == frozenset(rep)]
            label_groups = orbits_under(stab, list(rep))
            for isotopes in product(*(fragment.site(grp[0]).label_isotopes() for grp in label_groups)):
                raw = float(len(orbit))
                for iso, grp in zip(isotopes, label_groups):
                    for x in grp:
                        raw *= labeling.abundance(x, iso)
                for x in labelable:
                    if x not in rep:
                        raw *= unlabelled[x]
                entries.append((list(zip(isotopes, [tuple(g) for g in label_groups])), stab, raw, rep, len(orbit)))
    top_raw = max((e[2] for e in entries), default=0.0)
    if top_raw <= 0:
        raise ValueError(f"Fragment {fragment.name!r} has no labelled site with nonzero abundance.")

    components, labels, ratios, omitted, label_sets = [], [], [], [], []
    names: Dict[str, List[str]] = {}
    unspecified = set()
    for groups, stab, raw, rep, n_equiv in entries:
        label = "+".join(f"{iso}@{','.join(grp)}" for iso, grp in groups) + (f" (x{int(round(n_equiv))})"
                                                                   if round(n_equiv) > 1 else "")
        relative = raw / top_raw
        if relative < min_ratio:
            omitted.append((label, f"relative amplitude {relative:.3g} below min_ratio {min_ratio:.3g}"))
            continue
        merged = orbits_under(stab, proton_labels)
        if all(fragment.coupling(grp[0], g[0], orbits) is None for _, grp in groups for g in merged):
            omitted.append((label, "no coupling to an included proton group is given"))
            continue
        def layout(split: bool):
            label_parts = [(iso, (x,)) for iso, grp in groups for x in grp] if split else list(groups)
            proton_parts = [(x,) for g in merged for x in g] if split else [tuple(g) for g in merged]
            spins_ = [grp for _, grp in label_parts] + proton_parts
            isos_ = [iso for iso, _ in label_parts] + ["1H"] * len(proton_parts)
            sizes_ = [len(grp) for _, grp in label_parts] + [sum(fragment.proton(x).size for x in g)
                                                             for g in proton_parts]
            n_ = len(spins_)
            matrix_ = np.zeros((n_, n_))
            keys_: Dict[Tuple[int, int], Pair] = {}
            missing_ = set()
            for a in range(n_):
                for b in range(a + 1, n_):
                    reps = {orbits[pair(x, y)] for x in spins_[a] for y in spins_[b]}
                    if len(reps) != 1:
                        return None
                    key = reps.pop()
                    value = fragment.coupling(spins_[a][0], spins_[b][0], orbits)
                    if value is None:
                        missing_.add(fragment.key_name(key))
                        value = 0.0
                    matrix_[a, b] = matrix_[b, a] = value
                    keys_[(a, b)] = key
            return isos_, sizes_, matrix_, keys_, missing_

        built = layout(split=False)
        if built is None:
            # Chemically but not magnetically equivalent (e.g. both methyl carbons of isopropyl labelled):
            # keep the spins separate; symmetric couplings still share one parameter through their orbit key.
            built = layout(split=True)
            label += " [not magnetically equivalent]"
        if built is None:
            raise ValueError(f"{label}: inconsistent couplings between symmetry-related groups.")
        isos, sizes, matrix, keys, missing = built
        unspecified.update(missing)
        system = SpinSystem.from_group_couplings(isos, sizes, matrix)
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
        components.append(system)
        labels.append(label)
        ratios.append(raw)
        label_sets.append(frozenset(rep))
    if not components:
        raise ValueError(f"Fragment {fragment.name!r} gives no observable isotopologue.")
    # Primary components carry free amplitudes; minor ones follow the primary sharing most labelled sites.
    strongest = max(ratios)
    primary = [i for i, r in enumerate(ratios) if r >= labeling.primary_fraction * strongest]
    parents: List[Optional[int]] = []
    for i in range(len(components)):
        if i in primary:
            parents.append(None)
            continue
        parents.append(max(primary, key=lambda j: (len(label_sets[i] & label_sets[j]), ratios[j])))
    top = ratios[parents.index(None)]
    ratios = tuple(r / top for r in ratios)
    interp = Interpretation(tuple(Component(s, r, l) for s, r, l in zip(components, ratios, labels)))
    names = {k: v for k, v in names.items() if v}
    return HypothesisModel(name or fragment.name, fragment, interp, labels, ratios, names,
                           tuple(sorted(k for k in unspecified if k in names)), omitted, None, [], parents)


def min_ratio_for_snr(peak_snr: float, significance: float = 2.0) -> float:
    """Smallest relative amplitude whose strongest lines could reach `significance` sigma, given the peak SNR of
    the strongest observed band (used to gate minor isotopologues)."""
    return significance / max(peak_snr, 1e-12)


def combine_models(models: Sequence[HypothesisModel], name: Optional[str] = None) -> HypothesisModel:
    """Separate molecules (or fragments without mutual couplings) as one candidate.

    Components are concatenated and renumbered; coupling keys get a part
    prefix ('P1:J(C1,H1)'); ratios stay fixed only within each part
    (`ratio_blocks`) and decay rates are tied only within a part.
    """
    import re
    if len(models) == 1:
        return models[0]
    comps, labels, ratios, names, unspecified, omitted, blocks, parents = [], [], [], {}, [], [], [], []
    for p, m in enumerate(models):
        shift = len(comps)
        rename = lambda n: re.sub(r"^c(\d+)\.", lambda x: f"c{int(x.group(1)) + shift}.", n)
        comps.extend(m.interpretation.components)
        labels.extend(f"P{p + 1}:{l}" for l in m.component_labels)
        ratios.extend(m.ratios)
        for key, ns in m.coupling_names.items():
            names[f"P{p + 1}:{key}"] = [rename(n) for n in ns]
        unspecified.extend(f"P{p + 1}:{k}" for k in m.unspecified)
        omitted.extend((f"P{p + 1}:{l}", r) for l, r in m.omitted)
        blocks.append(list(range(shift, shift + len(m.component_labels))))
        parents.extend(None if m.parent(c) is None else m.parent(c) + shift for c in range(len(m.component_labels)))
    interp = Interpretation(tuple(comps))
    return HypothesisModel(name or " | ".join(m.name for m in models), models[0].fragment, interp, labels,
                           tuple(ratios), names, tuple(unspecified), omitted, blocks, list(models), parents)
