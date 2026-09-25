"""Local identifiability and solver basin measurements.

`local_identifiability` linearizes the processed spectrum with respect to the
free parameters at a given point (gains re-solved by variable projection) and
reports singular values: near-zero values mark directions the observation
cannot resolve. `basin_of_attraction` starts the solver from the truth plus
perturbations of increasing size and records success rates; these numbers set
tolerable network J errors and bin widths (Phase 2 of the plan).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

import numpy as np

from ..evaluation.matching import coupling_error
from ..physics.protocol import SUDDEN_DROP, Protocol
from ..solver.forward import MixtureForward
from ..solver.observed import ObservedSpectrum
from ..solver.parameterization import Parameterization, ParameterPolicy
from ..solver.refine import RefineSettings, refine
from ..spinsystem import Component, Interpretation, SpinSystem


@dataclass
class IdentifiabilityReport:
    parameter_names: List[str]
    singular_values: List[float]
    weakest_direction: Dict[str, float]
    condition_number: float
    relative_threshold: float
    unresolved_directions: int

    def to_dict(self) -> dict:
        return self.__dict__.copy()


def local_identifiability(interpretation: Interpretation, observed: ObservedSpectrum,
                          policy: ParameterPolicy = ParameterPolicy(), protocol: Protocol = SUDDEN_DROP,
                          step_hz: float = 1e-3, relative_threshold: float = 1e-3,
                          couplings_only: bool = True) -> IdentifiabilityReport:
    param = Parameterization.from_interpretation(interpretation, policy)
    if couplings_only:
        param.fix(*[n for n in param.free_names if param.parameters[n].kind != "coupling"])
    forward = MixtureForward(param, observed, protocol, gain_model="complex")
    x0 = param.vector()
    r0 = forward.predict(x0).residual
    jac = np.zeros((len(r0), len(x0)))
    for i in range(len(x0)):
        x = x0.copy()
        x[i] += step_hz
        jac[:, i] = (forward.predict(x).residual - r0) / step_hz
    u, s, vt = np.linalg.svd(jac, full_matrices=False)
    weakest = dict(zip(param.free_names, vt[-1].tolist())) if len(s) else {}
    cond = float(s[0] / s[-1]) if len(s) and s[-1] > 0 else float("inf")
    unresolved = int(np.sum(s < relative_threshold * s[0])) if len(s) else 0
    return IdentifiabilityReport(param.free_names, s.tolist(), weakest, cond, relative_threshold, unresolved)


def perturb_couplings(system: SpinSystem, scale_hz: float, rng: np.random.Generator) -> SpinSystem:
    groups = system.groups
    j = system.couplings_hz.copy()
    for a in range(len(groups)):
        for b in range(a + 1, len(groups)):
            d = rng.normal(0, scale_hz)
            j[np.ix_(groups[a], groups[b])] += d
            j[np.ix_(groups[b], groups[a])] += d
    return SpinSystem(system.isotopes, j, groups)


def basin_of_attraction(truth: Interpretation, observed: ObservedSpectrum, scales_hz: Sequence[float],
                        trials: int = 5, tolerance_hz: float = 0.05, settings: RefineSettings = RefineSettings(starts=1),
                        protocol: Protocol = SUDDEN_DROP, seed: int = 0) -> List[dict]:
    rng = np.random.default_rng(seed)
    rows = []
    for scale in scales_hz:
        successes, evaluations = 0, []
        for _ in range(trials):
            start = Interpretation(tuple(Component(perturb_couplings(c.system, scale, rng), c.contribution, c.label)
                                         for c in truth.components))
            margin = max(settings.policy.coupling_margin_hz, 4 * scale)
            local = RefineSettings(**{**settings.__dict__,
                                      "policy": ParameterPolicy(**{**settings.policy.__dict__,
                                                                   "coupling_margin_hz": margin})})
            result = refine(start, observed, local, protocol)
            error = max(coupling_error(t.system, r.system)[0]
                        for t, r in zip(truth.components, result.interpretation.components))
            successes += error <= tolerance_hz
            evaluations.append(result.evaluations)
        rows.append({"scale_hz": scale, "success_rate": successes / trials,
                     "median_evaluations": float(np.median(evaluations)), "trials": trials,
                     "tolerance_hz": tolerance_hz})
    return rows
