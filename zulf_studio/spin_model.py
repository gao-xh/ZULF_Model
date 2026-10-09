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
from pathlib import Path
from typing import Dict, List

import numpy as np

from zulf_hypothesis.spin_system import NAME, entry_key  # noqa: F401  (the fit's grouping and keys)


class CustomModel:
    """The model of a spin-system specification for the session: zulf_hypothesis.spin_system.spin_system_model
    (the fit's model: components, equivalence groups, keys) and its couplings for the sliders."""

    def __init__(self, spec: dict):
        from zulf_hypothesis.spin_system import couplings, spin_system_model
        m = spin_system_model(spec)
        self.component_labels, self.interpretation = m.component_labels, m.interpretation
        self.coupling_names, self.unspecified = m.coupling_names, set(m.unspecified)
        self._couplings = couplings(spec)

    def couplings(self) -> List[dict]:
        return [{"key": c["key"], "value": c["value"], "one_bond": abs(c["value"]) >= 50.0, "unspecified": False,
                 "overridden": False, "variable": c["variable"]} for c in self._couplings]


def set_values(spec: dict, values: Dict[str, float]) -> dict:
    from zulf_hypothesis.spin_system import set_values as shared
    return shared(spec, values)


def from_structure_model(model, name: str = "") -> dict:
    """The spin-system specification of a structure model: every isotopologue a component with its abundance
    weight; each coupling of the structure becomes one variable named after its key (J(C1,HC1) -> J_C1_HC1),
    shared by every place it appears (equivalent protons, several isotopologues), so a fit of the spin system
    has the structure's parameters. Couplings that are equal only by value stay separate."""
    param_key = {n: k for k, names in model.coupling_names.items() for n in names}
    variables, comps = {}, []
    for c, (lab, comp) in enumerate(zip(model.component_labels, model.interpretation.components)):
        system = comp.system
        J = np.asarray(system.couplings_hz, float)
        n = len(system.isotopes)
        cells = [[0] * n for _ in range(n)]
        gs = system.groups
        for a in range(len(gs)):
            for b in range(a + 1, len(gs)):
                key = param_key.get(f"c{c}.J{a}-{b}")
                value = float(J[gs[a][0], gs[b][0]])
                if key is None and value == 0.0:
                    continue
                var = re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_") if key else f"J_{lab}_{a}_{b}"
                if not NAME.match(var):
                    var = "J_" + var
                variables.setdefault(var, round(value, 6))
                for p in gs[a]:
                    for q in gs[b]:
                        cells[p][q] = cells[q][p] = var
        comps.append({"name": lab, "isotopes": list(system.isotopes), "J": cells, "weight": float(comp.contribution)})
    # abundance ratios are known for a structure: the weights are held (free amplitudes let a fit from a poor
    # start switch one isotopologue off; ethanol: objective 0.78 instead of 0.12)
    return {"compound": name or "custom",
            "spin_system": {"components": comps, "variables": variables, "fixed_weights": True}}


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
