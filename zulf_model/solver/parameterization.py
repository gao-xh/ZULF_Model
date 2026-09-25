"""General parameterization of candidate interpretations for refinement.

Couplings are parameterized between magnetic-equivalence groups (couplings
inside a group are unobservable). Every parameter can be free, fixed, bounded
or tied to others by name. Ties express shared couplings across components
(for example H-H couplings shared by isotopologues of one molecule). Decay
rates are per component or per transition family; Gaussian widths and a
global phase delay are optional nonlinear parameters. Nothing here depends on
a particular molecule, spin count or nucleus list.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..spinsystem import Component, Interpretation, SpinSystem


@dataclass
class Parameter:
    name: str
    value: float
    lower: float
    upper: float
    free: bool = True
    kind: str = "coupling"          # coupling | log_rate | sigma | phase_delay
    component: int = -1
    detail: tuple = ()               # (group_a, group_b) for couplings, (family,) for rates

    def check(self) -> None:
        if not (math.isfinite(self.value) and self.lower <= self.value <= self.upper and self.lower < self.upper):
            raise ValueError(f"Parameter {self.name}: need lower <= value <= upper and lower < upper.")


@dataclass(frozen=True)
class ParameterPolicy:
    coupling_margin_hz: float = 5.0
    coupling_margin_relative: float = 0.05
    fix_small_couplings_below_hz: float = 0.0     # 0 keeps every coupling free
    rate_bounds_per_s: Tuple[float, float] = (0.05, 50.0)
    initial_rate_per_s: float = 1.0
    family_edges_hz: Tuple[float, ...] = ()        # shared transition-family split, empty = one family
    fit_sigma: bool = False
    sigma_bounds_hz: Tuple[float, float] = (1e-3, 5.0)
    initial_sigma_hz: float = 0.05
    fit_phase_delay: bool = False
    phase_delay_bounds_s: Tuple[float, float] = (-0.01, 0.01)

    @classmethod
    def from_dict(cls, data: dict) -> "ParameterPolicy":
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items()
                      if k in cls.__dataclass_fields__})


@dataclass
class ComponentLayout:
    group_nuclei: List[str]
    group_sizes: List[int]
    label: str
    contribution: float


class Parameterization:
    def __init__(self, layouts: List[ComponentLayout], parameters: List[Parameter], ties: Dict[str, str],
                 policy: ParameterPolicy):
        self.layouts = layouts
        self.parameters = {p.name: p for p in parameters}
        self.order = [p.name for p in parameters]
        self.ties = dict(ties)       # follower name -> leader name
        self.policy = policy
        for p in parameters:
            p.check()

    # -- construction ---------------------------------------------------------------
    @classmethod
    def from_interpretation(cls, interpretation: Interpretation, policy: ParameterPolicy = ParameterPolicy(),
                            rate_guess: Optional[Sequence[float]] = None) -> "Parameterization":
        layouts, params = [], []
        for c_index, component in enumerate(interpretation.components):
            system = component.system
            groups = list(system.groups)
            gj = system.group_couplings()
            layouts.append(ComponentLayout([system.isotopes[g[0]] for g in groups], [len(g) for g in groups],
                                           component.label, component.contribution))
            for a in range(len(groups)):
                for b in range(a + 1, len(groups)):
                    value = float(gj[a, b])
                    margin = policy.coupling_margin_hz + policy.coupling_margin_relative * abs(value)
                    free = abs(value) >= policy.fix_small_couplings_below_hz or policy.fix_small_couplings_below_hz == 0
                    params.append(Parameter(f"c{c_index}.J{a}-{b}", value, value - margin, value + margin, free,
                                            "coupling", c_index, (a, b)))
            families = len(policy.family_edges_hz) + 1
            for f in range(families):
                guess = rate_guess[c_index] if rate_guess is not None else policy.initial_rate_per_s
                lo, hi = policy.rate_bounds_per_s
                guess = float(np.clip(guess, lo, hi))
                params.append(Parameter(f"c{c_index}.log_rate{f}", math.log(guess), math.log(lo), math.log(hi), True,
                                        "log_rate", c_index, (f,)))
            if policy.fit_sigma:
                lo, hi = policy.sigma_bounds_hz
                params.append(Parameter(f"c{c_index}.sigma", float(np.clip(policy.initial_sigma_hz, lo, hi)), lo, hi,
                                        True, "sigma", c_index))
        if policy.fit_phase_delay:
            lo, hi = policy.phase_delay_bounds_s
            params.append(Parameter("phase_delay", 0.0, lo, hi, True, "phase_delay"))
        return cls(layouts, params, {}, policy)

    # -- editing ------------------------------------------------------------------------
    def fix(self, *names: str) -> "Parameterization":
        for n in names:
            self.parameters[n].free = False
        return self

    def release(self, *names: str) -> "Parameterization":
        for n in names:
            self.parameters[n].free = True
        return self

    def set(self, name: str, value: Optional[float] = None, lower: Optional[float] = None,
            upper: Optional[float] = None) -> "Parameterization":
        p = self.parameters[name]
        p.value = p.value if value is None else float(value)
        p.lower = p.lower if lower is None else float(lower)
        p.upper = p.upper if upper is None else float(upper)
        p.check()
        return self

    def tie(self, leader: str, *followers: str) -> "Parameterization":
        """Force followers to equal the leader (shared couplings across components)."""
        for f in followers:
            if f == leader or f not in self.parameters:
                raise ValueError(f"Cannot tie {f} to {leader}.")
            self.ties[f] = leader
            self.parameters[f].value = self.parameters[leader].value
        return self

    def coupling_names(self, component: int) -> List[str]:
        return [n for n in self.order if self.parameters[n].kind == "coupling" and self.parameters[n].component == component]

    # -- vector mapping ------------------------------------------------------------------
    @property
    def free_names(self) -> List[str]:
        return [n for n in self.order if self.parameters[n].free and n not in self.ties]

    def vector(self) -> np.ndarray:
        return np.array([self.parameters[n].value for n in self.free_names])

    def bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        return (np.array([self.parameters[n].lower for n in self.free_names]),
                np.array([self.parameters[n].upper for n in self.free_names]))

    def values(self, x: Optional[np.ndarray] = None) -> Dict[str, float]:
        out = {n: p.value for n, p in self.parameters.items()}
        if x is not None:
            out.update(dict(zip(self.free_names, np.asarray(x, float))))
        for follower, leader in self.ties.items():
            out[follower] = out[leader]
        return out

    def update(self, x: np.ndarray) -> None:
        for n, v in self.values(x).items():
            self.parameters[n].value = v

    # -- model quantities ------------------------------------------------------------------
    def systems(self, values: Dict[str, float]) -> List[SpinSystem]:
        out = []
        for c, layout in enumerate(self.layouts):
            g = len(layout.group_sizes)
            gj = np.zeros((g, g))
            for n in self.coupling_names(c):
                a, b = self.parameters[n].detail
                gj[a, b] = gj[b, a] = values[n]
            out.append(SpinSystem.from_group_couplings(layout.group_nuclei, layout.group_sizes, gj))
        return out

    def rates(self, values: Dict[str, float], component: int) -> np.ndarray:
        families = len(self.policy.family_edges_hz) + 1
        return np.array([math.exp(values[f"c{component}.log_rate{f}"]) for f in range(families)])

    def sigma(self, values: Dict[str, float], component: int) -> float:
        return values.get(f"c{component}.sigma", 0.0)

    def phase_delay(self, values: Dict[str, float]) -> float:
        return values.get("phase_delay", 0.0)

    def interpretation(self, values: Dict[str, float], contributions: Sequence[float]) -> Interpretation:
        comps = [Component(system, float(max(c, 0.0)), layout.label, {"refined": True})
                 for system, layout, c in zip(self.systems(values), self.layouts, contributions)]
        return Interpretation(tuple(comps))

    def boundary_hits(self, x: np.ndarray, fraction: float = 0.01) -> List[str]:
        lo, hi = self.bounds()
        width = hi - lo
        return [n for n, v, l, h, w in zip(self.free_names, x, lo, hi, width) if min(v - l, h - v) < fraction * w]
