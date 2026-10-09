"""Custom spin systems for ZULF Studio (PLAN 8c): isotopes and a J matrix typed in, instead of a structure.

Specification (also what a session file and save_model store):

    {"compound": "name", "spin_system": {
        "components": [{"name": "A", "isotopes": ["13C", "1H", "1H", "1H"],
                        "J": [[0, "a", "a", "a"], ["a", 0, 0, 0], ...], "weight": 1.0}, ...],
        "variables": {"a": 125.0}}}

J entries are numbers (Hz) or variable names; one name in several places ties those couplings together (the
equivalent protons of a methyl group share one J). The matrix is read from its upper triangle and made symmetric.
Every component is one spin system (zulf_core.SpinSystem: exact diagonalisation, equivalent groups found
automatically) with a weight; the spectrum is the weighted sum. ZULF_NMR_Suite molecule folders (structure.csv:
a row of isotopes, then the n x n matrix in Hz) can be read in.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np

from zulf_core.nuclei import get_registry
from zulf_core.spinsystem import SpinSystem

NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


@dataclass
class _Component:
    system: SpinSystem
    contribution: float


@dataclass
class _Interpretation:
    components: List[_Component]


class CustomModel:
    """The model of a spin-system specification, shaped like a structure model for the session (component labels,
    interpretation.components with .system and .contribution, coupling_names, unspecified)."""

    def __init__(self, spec: dict):
        ss = spec["spin_system"]
        self.variables = {str(k): float(v) for k, v in (ss.get("variables") or {}).items()}
        comps = ss.get("components") or []
        if not comps:
            raise ValueError("a spin system needs at least one component")
        registry = get_registry()
        self.component_labels, parts, self.entries = [], [], []
        used = set()
        for c in comps:
            iso = [str(x).strip() for x in c["isotopes"]]
            for x in iso:
                registry[x]                                   # unknown isotopes raise here
            n = len(iso)
            J = c.get("J") or [[0] * n for _ in range(n)]
            if len(J) != n or any(len(row) != n for row in J):
                raise ValueError(f"component {c.get('name')}: J must be {n} x {n} for {n} isotopes")
            mat = np.zeros((n, n))
            for i in range(n):
                for j in range(i + 1, n):
                    v = J[i][j]
                    if isinstance(v, str) and v.strip() == "":
                        v = 0.0
                    if isinstance(v, str) and not _is_number(v):
                        name = v.strip()
                        if not NAME.match(name):
                            raise ValueError(f"J[{i + 1},{j + 1}] = {v!r} is neither a number nor a variable name")
                        if name not in self.variables:
                            raise ValueError(f"variable {name!r} has no value (variables: {sorted(self.variables)})")
                        used.add(name)
                        mat[i, j] = mat[j, i] = self.variables[name]
                    else:
                        mat[i, j] = mat[j, i] = float(v)
                        if float(v) != 0.0:
                            self.entries.append((str(c.get("name", "")), i, j, float(v)))
            label = str(c.get("name") or f"system {len(parts) + 1}")
            self.component_labels.append(label)
            parts.append(_Component(SpinSystem(tuple(iso), mat), float(c.get("weight", 1.0))))
        self.interpretation = _Interpretation(parts)
        self.coupling_names = [v for v in self.variables if v in used] + [entry_key(c, i, j)
                                                                           for c, i, j, _ in self.entries]
        self.unspecified = set()

    def couplings(self) -> List[dict]:
        out = [{"key": k, "value": self.variables[k], "one_bond": abs(self.variables[k]) >= 50.0,
                "unspecified": False, "overridden": False, "variable": True} for k in self.coupling_names
               if k in self.variables]
        out += [{"key": entry_key(c, i, j), "value": v, "one_bond": abs(v) >= 50.0, "unspecified": False,
                 "overridden": False, "variable": False} for c, i, j, v in self.entries]
        return out


def entry_key(component: str, i: int, j: int) -> str:
    """Slider key of a numeric matrix entry: component name and 1-based spin numbers, e.g. 'A:J(1,2)'."""
    return f"{component}:J({i + 1},{j + 1})"


def _is_number(text: str) -> bool:
    try:
        float(text)
        return True
    except ValueError:
        return False


def set_values(spec: dict, values: Dict[str, float]) -> dict:
    """A copy of the specification with variables and numeric entries changed (keys as in CustomModel)."""
    import copy
    out = copy.deepcopy(spec)
    ss = out["spin_system"]
    ss.setdefault("variables", {})
    comps = {str(c.get("name", "")): c for c in ss["components"]}
    for key, v in values.items():
        m = re.match(r"^(.*):J\((\d+),(\d+)\)$", key)
        if m:
            c = comps.get(m.group(1))
            if c is None:
                raise ValueError(f"no component {m.group(1)!r}")
            i, j = int(m.group(2)) - 1, int(m.group(3)) - 1
            c["J"][i][j] = c["J"][j][i] = float(v)
        elif key in ss["variables"]:
            ss["variables"][key] = float(v)
        else:
            raise ValueError(f"unknown coupling {key!r} (variables {sorted(ss['variables'])} or 'name:J(i,j)')")
    return out


def from_structure_model(model, name: str = "") -> dict:
    """The spin-system specification of a structure model: every isotopologue a component with its abundance
    weight; couplings with the same value (equivalent protons, the same coupling in several isotopologues) share
    one variable J1, J2, ... (largest first), so they move together."""
    values = sorted({round(float(v), 6) for c in model.interpretation.components
                     for v in np.asarray(c.system.couplings_hz, float)[np.triu_indices(len(c.system.isotopes), 1)]
                     if abs(v) > 0}, key=lambda v: -abs(v))
    names = {v: f"J{k + 1}" for k, v in enumerate(values)}
    comps = []
    for lab, c in zip(model.component_labels, model.interpretation.components):
        J = np.asarray(c.system.couplings_hz, float)
        n = len(c.system.isotopes)
        cells = [[names.get(round(float(J[i, j]), 6), 0) if i != j and abs(J[i, j]) > 0 else 0 for j in range(n)]
                 for i in range(n)]
        comps.append({"name": lab, "isotopes": list(c.system.isotopes), "J": cells, "weight": float(c.contribution)})
    return {"compound": name or "custom",
            "spin_system": {"components": comps, "variables": {n: v for v, n in names.items()}}}


def read_suite_molecule(folder) -> dict:
    """A ZULF_NMR_Suite molecule folder (structure.csv: isotopes, then the n x n J matrix in Hz)."""
    folder = Path(folder)
    path = folder / "structure.csv" if folder.is_dir() else folder
    rows = [r for r in csv.reader(open(path)) if any(x.strip() for x in r)]
    iso = [x.strip() for x in rows[0]]
    J = [[float(x) for x in r] for r in rows[1:1 + len(iso)]]
    name = (folder if folder.is_dir() else folder.parent).name
    return {"compound": name, "spin_system": {"components": [{"name": name, "isotopes": iso, "J": J, "weight": 1.0}],
                                              "variables": {}}}
