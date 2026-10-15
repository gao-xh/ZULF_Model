"""Structure specifications as JSON: a registered motif or a chain of groups, and coupling overrides by key.

    {"motif": "pyridine ring", "one_bond": {...}}
    {"compound": "acetonitrile", "chain": {"groups": [["C1", "C", 3], ["C2", "C", 0]], "bonds": [["C1", "C2"]],
                                           "sym": [...]}, "one_bond": {"C1": 136.3}}

Used by the fitting scripts (`--structure`), the regression over the confirmed samples and the studio
(zulf_studio). Coupling keys are 'J(a,b)' with site or proton-group labels of the fragment.
"""
from __future__ import annotations

from typing import Dict

from .fragment import Fragment, pair
from .motifs import MOTIFS, _chain


def fragment_from_spec(spec: dict) -> Fragment:
    """The fragment of a structure specification (a registered motif, or a chain of (site, element, H count)
    groups with bonds and optional symmetry maps)."""
    if "motif" in spec:
        return MOTIFS[spec["motif"]].fragment(spec.get("one_bond", {}))
    c = spec["chain"]
    return _chain(spec.get("compound", "chain"), [tuple(g) for g in c["groups"]], [tuple(b) for b in c["bonds"]],
                  tuple(c.get("sym", ())))(spec.get("one_bond", {}))


def coupling_key(a: str, b: str) -> str:
    return f"J({a},{b})"


def parse_coupling_key(key: str):
    """'J(a,b)' -> (a, b)."""
    if not (key.startswith("J(") and key.endswith(")")) or key.count(",") != 1:
        raise ValueError(f"Coupling key must look like 'J(a,b)', got {key!r}.")
    a, b = key[2:-1].split(",")
    return a.strip(), b.strip()


def override_couplings(fragment: Fragment, couplings: Dict[str, float]) -> Fragment:
    """Set couplings by key 'J(a,b)'; every symmetry image of the pair gets the same value."""
    for key, value in couplings.items():
        a, b = parse_coupling_key(key)
        orbit = {pair(a, b)}
        grown = True
        while grown:
            grown = False
            for g in fragment.symmetry:
                for p in list(orbit):
                    x, y = tuple(p)
                    q = pair(g.get(x, x), g.get(y, y))
                    if q not in orbit:
                        orbit.add(q)
                        grown = True
        for p in orbit:
            fragment.couplings[p] = float(value)
    fragment.validate()
    return fragment
