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
    kind: str = "coupling"          # coupling | log_rate | sigma | phase_delay | log_exchange | field | nuisance
    component: int = -1
    detail: tuple = ()               # (group_a, group_b) for couplings, (family,) for rates, (group,) for exchange

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
    # First-order phase as a time delay. Fitted by default: every acquisition has one (dead time,
    # switching), and a shared zero-order phase alone forces the couplings to absorb it (D24).
    fit_phase_delay: bool = True
    phase_delay_bounds_s: Tuple[float, float] = (-0.01, 0.01)
    # Static field during evolution (Protocol.field_ut, D19), fitted only when asked; the default is zero field.
    # With preparation and detection along z the signal depends on the field only through its transverse size
    # B_transverse (placed in Bx) and |Bz|, so two parameters >= 0 cover every field ("field_transverse_ut", "field_z_ut").
    # The signal is stationary at zero field (zero gradient), so the starts must be nonzero.
    fit_field: bool = False
    field_axes: Tuple[str, ...] = ("transverse", "z")
    field_bounds_ut: Tuple[float, float] = (0.0, 1.0)
    initial_field_ut: Tuple[float, float] = (0.02, 0.02)     # (transverse, z) starts
    field_fixed: bool = False                                  # hold the field at initial_field_ut (a known field)
    # A second field region (D63): components whose label matches this regex (e.g. "^P2:", the second copy of a
    # combined model) evolve in their own fitted field ("field2_transverse_ut", "field2_z_ut"; same axes and
    # bounds), the others in the first. Models a sample that sees two field regions. Empty: one field.
    second_field_components: str = ""
    initial_field2_ut: Tuple[float, float] = (0.1, 0.05)       # (transverse, z) starts of the second field
    # Fitted gyromagnetic ratios (D53): gamma / (2 pi) in Hz/uT of the named nuclei as free parameters, used by the
    # field term and the gamma weights. Only meaningful with a known (fixed) field: the splittings scale with
    # gamma B, so gamma and a free field trade against each other.
    fit_gamma: Tuple[str, ...] = ()
    gamma_bounds_hz_per_ut: Tuple[float, float] = (-50.0, 50.0)
    initial_gamma: Tuple[Tuple[str, float], ...] = ()          # starts (default: the registry value)
    # Nuisance terms rendered through the same operator; linear amplitudes are real and unconstrained.
    # {"kind": "exponential", "rate_bounds_per_s": [lo, hi], "initial_rate_per_s": r}
    # {"kind": "damped_sinusoid", "frequency_bounds_hz": [lo, hi], "initial_frequency_hz": f,
    #  "rate_bounds_per_s": [lo, hi], "initial_rate_per_s": r}
    # {"kind": "template", "template": [...full-record samples...], "shift_bounds_s": [lo, hi]}
    nuisance: Tuple[dict, ...] = ()

    @classmethod
    def from_dict(cls, data: dict) -> "ParameterPolicy":
        data = {k: v for k, v in data.items() if k in cls.__dataclass_fields__}
        if "nuisance" in data:
            data["nuisance"] = tuple(dict(t) for t in data["nuisance"])
        return cls(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items()})


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
            start = 0.0 if lo <= 0.0 <= hi else 0.5 * (lo + hi)     # an instrument prior may exclude zero
            params.append(Parameter("phase_delay", start, lo, hi, True, "phase_delay"))
        if policy.fit_field:
            lo, hi = policy.field_bounds_ut
            for axis in policy.field_axes:
                if axis not in ("transverse", "z"):
                    raise ValueError(f"Unknown field axis '{axis}' (use 'transverse' and/or 'z').")
                start = float(np.clip(policy.initial_field_ut[0 if axis == "transverse" else 1], lo, hi))
                params.append(Parameter(f"field_{axis}_ut", start, lo, hi, not policy.field_fixed, "field", -1,
                                        (axis,)))
            if policy.second_field_components:
                import re
                pattern = re.compile(policy.second_field_components)
                if not any(pattern.search(lay.label) for lay in layouts):
                    raise ValueError(f"second_field_components {policy.second_field_components!r} matches no "
                                     f"component of {[lay.label for lay in layouts]}")
                for axis in policy.field_axes:
                    start = float(np.clip(policy.initial_field2_ut[0 if axis == "transverse" else 1], lo, hi))
                    params.append(Parameter(f"field2_{axis}_ut", start, lo, hi, True, "field", -1, (axis,)))
        if policy.fit_gamma:
            from ..nuclei import get_registry
            lo, hi = policy.gamma_bounds_hz_per_ut
            starts = dict(policy.initial_gamma)
            for symbol in policy.fit_gamma:
                start = float(np.clip(starts.get(symbol, get_registry().gamma(symbol)), lo, hi))
                params.append(Parameter(f"gamma_{symbol}", start, lo, hi, True, "gamma", -1, (symbol,)))
        for i, term in enumerate(policy.nuisance):
            kind = term.get("kind")
            if kind in ("exponential", "damped_sinusoid"):
                lo, hi = term.get("rate_bounds_per_s", (0.01, 100.0))
                r0 = float(np.clip(term.get("initial_rate_per_s", math.sqrt(lo * hi)), lo, hi))
                params.append(Parameter(f"n{i}.log_rate", math.log(r0), math.log(lo), math.log(hi), True,
                                        "nuisance", -1, (i,)))
            if kind == "damped_sinusoid":
                lo, hi = term["frequency_bounds_hz"]
                f0 = float(np.clip(term.get("initial_frequency_hz", (lo + hi) / 2), lo, hi))
                params.append(Parameter(f"n{i}.frequency", f0, lo, hi, True, "nuisance", -1, (i,)))
            if kind == "template" and term.get("shift_bounds_s"):
                lo, hi = term["shift_bounds_s"]
                params.append(Parameter(f"n{i}.shift_s", 0.0, lo, hi, True, "nuisance", -1, (i,)))
            if kind not in ("exponential", "damped_sinusoid", "template"):
                raise ValueError(f"Unknown nuisance kind '{kind}'.")
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

    def add_exchange(self, component: int, group: int, rate_per_s: float = 10.0,
                     bounds_per_s: Tuple[float, float] = (0.01, 1.0e5), free: bool = True,
                     name: Optional[str] = None) -> str:
        """Let every spin of one group of a component exchange with the solvent at a fitted rate k (1/s),
        parameterized as log k (physics.exchange, D45). Returns the parameter name; tie it across the
        isotopologues of one molecule with `tie`."""
        if not 0 <= component < len(self.layouts) or not 0 <= group < len(self.layouts[component].group_sizes):
            raise ValueError("No such component or group.")
        lo, hi = bounds_per_s
        name = name or f"c{component}.log_kex{group}"
        value = float(np.clip(rate_per_s, lo, hi))
        prm = Parameter(name, math.log(value), math.log(lo), math.log(hi), free, "log_exchange", component, (group,))
        prm.check()
        self.parameters[name] = prm
        self.order.append(name)
        return name

    def exchange_rates(self, values: Dict[str, float], component: int) -> Dict[int, float]:
        """{group: k (1/s)} of the exchanging groups of a component (empty: static spin system)."""
        return {self.parameters[n].detail[0]: math.exp(values[n]) for n in self.order
                if self.parameters[n].kind == "log_exchange" and self.parameters[n].component == component}

    def has_exchange(self, component: int) -> bool:
        return any(p.kind == "log_exchange" and p.component == component for p in self.parameters.values())

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

    def gamma_overrides(self, values: Dict[str, float]) -> Dict[str, float]:
        """{symbol: gamma / (2 pi) in Hz/uT} of the fitted gyromagnetic ratios (empty without them); a value at
        zero is kept 1e-6 away from it (zero gamma has no field term and no signal)."""
        out = {}
        for n, p in self.parameters.items():
            if p.kind == "gamma":
                g = float(values.get(n, p.value))
                out[p.detail[0]] = g if abs(g) >= 1e-6 else (1e-6 if g >= 0 else -1e-6)
        return out

    def has_field(self) -> bool:
        return any(p.kind == "field" for p in self.parameters.values())

    def in_second_field(self, component: Optional[int]) -> bool:
        """Whether a component evolves in the second field region (policy.second_field_components)."""
        pattern = self.policy.second_field_components
        if not pattern or component is None or not 0 <= component < len(self.layouts):
            return False
        import re
        return re.search(pattern, self.layouts[component].label) is not None

    def field_ut(self, values: Dict[str, float], component: Optional[int] = None
                 ) -> Optional[Tuple[float, float, float]]:
        """(Bx, By, Bz) in microtesla from the field parameters (B_transverse along x), or None without them. With
        a second field region, the field of that region for its components."""
        if not self.has_field():
            return None
        pre = "field2_" if self.in_second_field(component) else "field_"
        return (float(values.get(pre + "transverse_ut", 0.0)), 0.0, float(values.get(pre + "z_ut", 0.0)))

    def interpretation(self, values: Dict[str, float], contributions: Sequence[float]) -> Interpretation:
        comps = [Component(system, float(max(c, 0.0)), layout.label, {"refined": True})
                 for system, layout, c in zip(self.systems(values), self.layouts, contributions)]
        return Interpretation(tuple(comps))

    def boundary_hits(self, x: np.ndarray, fraction: float = 0.01) -> List[str]:
        lo, hi = self.bounds()
        width = hi - lo
        return [n for n, v, l, h, w in zip(self.free_names, x, lo, hi, width) if min(v - l, h - v) < fraction * w]
