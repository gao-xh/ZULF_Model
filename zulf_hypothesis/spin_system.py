"""Custom spin systems as fit models (PLAN 8d): isotopes and a J matrix typed in, instead of a structure.

    {"compound": "name", "spin_system": {
        "components": [{"name": "A", "isotopes": ["13C", "1H", "1H", "1H"],
                        "J": [[0, "a", "a", "a"], ["a", 0, 0, 0], ...], "weight": 1.0}, ...],
        "variables": {"a": 125.0}, "fixed_weights": false}}

A J entry is a number (Hz) or a variable name; one name in several places ties those couplings (a fit parameter
shared by every place, also across components). Spins of the same isotope whose rows of entries agree (by name
or number, outside their own pair) form an equivalence group; a coupling between two groups is one parameter
c<k>.J<a>-<b> (the solver's naming) under the key of its entry: the variable name, or 'component:J(i,j)' of the
groups' first spins for a number. Entries that are 0 are no coupling and stay 0. Component weights are fitted
amplitudes (each component its own) unless "fixed_weights" is true (then the weight ratios are held, variant
'ratios'). ZULF Studio uses the same grouping and keys for its sliders.
"""
from __future__ import annotations

import copy
import re
from typing import Dict, List, Optional, Tuple

import numpy as np

from zulf_core.nuclei import get_registry
from zulf_core.spinsystem import Component, Interpretation, SpinSystem

from .builder import HypothesisModel

NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
ENTRY = re.compile(r"^(.*):J\((\d+),(\d+)\)$")


def _token(v):
    """A matrix entry as a float or a variable name."""
    if isinstance(v, str):
        t = v.strip()
        if t == "":
            return 0.0
        try:
            return float(t)
        except ValueError:
            if not NAME.match(t):
                raise ValueError(f"J entry {v!r} is neither a number nor a variable name")
            return t
    return float(v)


def tokens(component: dict) -> List[List[object]]:
    """The symmetric token matrix of a component (upper triangle read, diagonal unused)."""
    iso = component["isotopes"]
    n = len(iso)
    J = component.get("J") or [[0] * n for _ in range(n)]
    if len(J) != n or any(len(row) != n for row in J):
        raise ValueError(f"component {component.get('name')}: J must be {n} x {n} for {n} isotopes")
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            out[i][j] = out[j][i] = _token(J[i][j])
    return out


def groups(isotopes, tok) -> List[Tuple[int, ...]]:
    """Equivalence groups by name: same isotope and the same entries to every other spin."""
    n = len(isotopes)
    out: List[List[int]] = []
    for i in range(n):
        for g in out:
            r = g[0]
            if isotopes[r] == isotopes[i] and all(tok[r][k] == tok[i][k] for k in range(n) if k not in (r, i)) \
                    and all(tok[m][i] == tok[r][i] for m in g[1:]):
                g.append(i)
                break
        else:
            out.append([i])
    return [tuple(g) for g in out]


def entry_key(component: str, i: int, j: int) -> str:
    return f"{component}:J({i + 1},{j + 1})"


def couplings(spec: dict) -> List[dict]:
    """One entry per fit parameter key: variables (with every place they appear), then the numeric couplings
    between equivalence groups. Each: key, value (Hz), variable (bool), places [(component, i, j)]."""
    ss = spec["spin_system"]
    variables = {str(k): float(v) for k, v in (ss.get("variables") or {}).items()}
    var_places: Dict[str, list] = {}
    numeric = []
    for c in ss["components"]:
        name = str(c.get("name", ""))
        tok = tokens(c)
        gs = groups(c["isotopes"], tok)
        for a in range(len(gs)):
            for b in range(a + 1, len(gs)):
                t = tok[gs[a][0]][gs[b][0]]
                if isinstance(t, str):
                    var_places.setdefault(t, []).append((name, gs[a][0], gs[b][0]))
                elif t != 0.0:
                    numeric.append({"key": entry_key(name, gs[a][0], gs[b][0]), "value": t, "variable": False,
                                    "places": [(name, gs[a][0], gs[b][0])]})
    missing = [v for v in var_places if v not in variables]
    if missing:
        raise ValueError(f"variables without a value: {missing}")
    return [{"key": v, "value": variables[v], "variable": True, "places": p} for v, p in var_places.items()] + numeric


def set_values(spec: dict, values: Dict[str, float]) -> dict:
    """A copy with variables and numeric couplings changed; a numeric key changes every entry of its group pair
    (all the equivalent spins), so equivalence is kept."""
    out = copy.deepcopy(spec)
    ss = out["spin_system"]
    ss.setdefault("variables", {})
    comps = {str(c.get("name", "")): c for c in ss["components"]}
    for key, v in values.items():
        m = ENTRY.match(key)
        if key in ss["variables"]:
            ss["variables"][key] = float(v)
        elif m:
            c = comps.get(m.group(1))
            if c is None:
                raise ValueError(f"no component {m.group(1)!r}")
            i, j = int(m.group(2)) - 1, int(m.group(3)) - 1
            gs = groups(c["isotopes"], tokens(c))
            gi = next(g for g in gs if i in g)
            gj = next(g for g in gs if j in g)
            for p in gi:
                for q in gj:
                    if p != q:
                        c["J"][p][q] = c["J"][q][p] = float(v)
        else:
            raise ValueError(f"unknown coupling {key!r} (variables {sorted(ss['variables'])} or 'name:J(i,j)')")
    return out


def numeric_matrix(component: dict, variables: Dict[str, float]) -> np.ndarray:
    tok = tokens(component)
    n = len(tok)
    m = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                t = tok[i][j]
                m[i, j] = variables[t] if isinstance(t, str) else t
    return m


def spin_system_model(spec: dict, overrides: Optional[Dict[str, float]] = None) -> HypothesisModel:
    """The HypothesisModel of a spin-system specification (fit_joint_series --structure, ZULF Studio)."""
    if overrides:
        spec = set_values(spec, overrides)
    ss = spec["spin_system"]
    comps = ss.get("components") or []
    if not comps:
        raise ValueError("a spin system needs at least one component")
    variables = {str(k): float(v) for k, v in (ss.get("variables") or {}).items()}
    registry = get_registry()
    labels, components, weights = [], [], []
    names: Dict[str, List[str]] = {}
    zero_keys = set()
    for k, c in enumerate(comps):
        iso = tuple(str(x).strip() for x in c["isotopes"])
        for x in iso:
            registry[x]                                  # an unknown isotope raises here
        tok = tokens(c)
        gs = groups(iso, tok)
        for t in {t for row in tok for t in row if isinstance(t, str)}:
            if t not in variables:
                raise ValueError(f"variable {t!r} has no value (variables: {sorted(variables)})")
        system = SpinSystem(iso, numeric_matrix(c, variables), groups=gs)
        label = str(c.get("name") or f"system {k + 1}")
        for a in range(len(gs)):
            for b in range(a + 1, len(gs)):
                t = tok[gs[a][0]][gs[b][0]]
                key = t if isinstance(t, str) else entry_key(label, gs[a][0], gs[b][0])
                names.setdefault(key, []).append(f"c{k}.J{a}-{b}")
                if not isinstance(t, str) and t == 0.0:
                    zero_keys.add(key)                   # no coupling: held at 0
        labels.append(label)
        components.append(system)
        weights.append(float(c.get("weight", 1.0)))
    if not any(w > 0 for w in weights):
        raise ValueError("every component weight is 0")
    top = next(w for w in weights if w > 0)
    ratios = tuple(w / top for w in weights)
    interp = Interpretation(tuple(Component(s, max(r, 0.0), lab) for s, r, lab in zip(components, ratios, labels)))
    blocks = None if ss.get("fixed_weights") else [[c] for c in range(len(labels))]
    return HypothesisModel(spec.get("compound", "spin system"), None, interp, labels, ratios, names,
                           tuple(sorted(zero_keys)), [], blocks, [], [None] * len(labels))
