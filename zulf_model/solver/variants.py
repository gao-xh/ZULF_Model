"""Sign variants of candidate interpretations for refinement.

A global sign change of every J leaves a zero-field spectrum unchanged, but
relative signs are observable whenever the coupling network connects the
spins, and the evidence is in small splittings that networks and local
refinement easily miss. A coupling whose refinement interval does not
contain zero can never change sign during refinement, so the opposite sign
has to be offered as a separate start. `sign_variants` returns candidates in
which one such coupling, or every coupling of one heteronuclear group, has
its sign flipped. Variants equivalent under the global sign are removed.
"""
from __future__ import annotations

from typing import List, Optional, Sequence

import numpy as np

from ..spinsystem import Component, Interpretation, SpinSystem
from .parameterization import ParameterPolicy


def crossing_threshold_hz(policy: ParameterPolicy) -> float:
    """|J| above which the policy interval (value +/- margin) excludes zero."""
    rel = policy.coupling_margin_relative
    return policy.coupling_margin_hz / max(1e-9, 1.0 - rel)


def _flipped(system: SpinSystem, pairs: Sequence[tuple]) -> SpinSystem:
    groups = system.groups
    gj = system.group_couplings().copy()
    for a, b in pairs:
        gj[a, b] = gj[b, a] = -gj[a, b]
    return SpinSystem.from_group_couplings([system.isotopes[g[0]] for g in groups], [len(g) for g in groups], gj)


def sign_variants(interpretation: Interpretation, min_abs_hz: float, max_variants: int = 16,
                  nuclei_order: Sequence[str] = ("13C", "15N", "1H")) -> List[Interpretation]:
    """Variants with one large coupling flipped, or all couplings of one heteronuclear group flipped.

    Large means |J| >= `min_abs_hz` (use `crossing_threshold_hz(policy)`). Ordered by the
    magnitude of the largest flipped coupling; the original is not included.
    """
    proposals = []
    for c, component in enumerate(interpretation.components):
        system = component.system
        gj = system.group_couplings()
        g = len(system.groups)
        for a in range(g):
            for b in range(a + 1, g):
                if abs(gj[a, b]) >= min_abs_hz:
                    proposals.append((abs(gj[a, b]), c, [(a, b)], f"c{c}.J{a}-{b}"))
        for a in range(g):
            if system.isotopes[system.groups[a][0]] == "1H":
                continue
            pairs = [(min(a, b), max(a, b)) for b in range(g) if b != a and gj[a, b] != 0]
            if len(pairs) > 1 and any(abs(gj[p]) >= min_abs_hz for p in pairs):
                proposals.append((max(abs(gj[p]) for p in pairs), c, pairs, f"c{c}.group{a}"))
    proposals.sort(key=lambda item: -item[0])
    out, seen = [], {interpretation.key(nuclei_order)}
    for _, c, pairs, label in proposals:
        comps = list(interpretation.components)
        original = comps[c]
        comps[c] = Component(_flipped(original.system, pairs), original.contribution, original.label,
                             dict(original.metadata))
        variant = Interpretation(tuple(comps), interpretation.score,
                                 dict(interpretation.metadata, sign_variant=label))
        key = variant.key(nuclei_order)
        if key in seen:
            continue
        seen.add(key)
        out.append(variant)
        if len(out) >= max_variants:
            break
    return out
