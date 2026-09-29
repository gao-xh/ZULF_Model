"""Linearised uncertainties of refined couplings: standard errors, correlations and poorly determined directions.

At the optimum x* of a fit, with J the Jacobian of the weighted residual r (variable projection: the linear gains,
shared phase and background are solved at every x, so J already accounts for them),

    cov(x) = s^2 (J^T J)^+,   s^2 = r.r / (n - k)

(s^2 absorbs misfit and correlated noise: the scale is the same quasi-likelihood c_hat the ranking uses). Standard
errors are sqrt(diag(cov)) in the parameter units (Hz for couplings), correlations cov_ij / (se_i se_j). The singular
values of J with columns scaled to unit norm show directions the data do not fix: a small value means a combination
of parameters (its right singular vector) that barely changes the spectrum. Linearisation holds near the optimum
only; a fit in a wrong local minimum has meaningless errors, so check starts and residuals first.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from zulf_core.physics.protocol import SUDDEN_DROP
from zulf_core.solver.forward import MixtureForward
from zulf_core.solver.refine import RefineSettings, _signal_kwargs

from .search import Evaluated, _settings_for


def coupling_uncertainties(e: Evaluated, observed, base: RefineSettings, protocol=SUDDEN_DROP,
                           weak_threshold: float = 1e-3) -> dict:
    """Standard errors (Hz), correlation matrix and weakly determined directions of the free couplings of one fit.

    `base` must be the RefineSettings the search used (fit_settings().base or the settings passed to
    search_hypotheses): weighting, background and ties change the Jacobian."""
    settings = _settings_for(e.model, e.variant, base)
    param = settings.parameterize(e.model.interpretation)
    values = e.summary["parameters"]
    free = param.free_names
    x = np.array([values[n] for n in free], float)
    forward = MixtureForward(param, observed, protocol, settings.gain_model, settings.background_order,
                             settings.band_weighting, signal_extra_hz=tuple(e.summary.get("signal_model_lines_hz") or ()),
                             **_signal_kwargs(settings))
    residual = forward.predict(x).residual
    jac = forward.jacobian(x)
    dof = max(len(residual) - e.k, 1)
    s2 = float(residual @ residual) / dof
    cov = s2 * np.linalg.pinv(jac.T @ jac)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    corr = cov / np.maximum(np.outer(se, se), 1e-300)
    norms = np.maximum(np.linalg.norm(jac, axis=0), 1e-300)
    _, sv, vt = np.linalg.svd(jac / norms, full_matrices=False)
    # key name ('J(a,b)') of each free coupling parameter
    key_of: Dict[str, str] = {}
    for key, names in e.model.coupling_names.items():
        for n in names:
            leader = param.ties.get(n, n)
            if leader in free and leader not in key_of:
                key_of[leader] = key
    idx = [i for i, n in enumerate(free) if n in key_of]
    keys = [key_of[free[i]] for i in idx]
    weak: List[dict] = []
    for s, v in zip(sv / sv[0], vt):
        if s < weak_threshold:
            parts = sorted(((abs(c), free[i], c) for i, c in enumerate(v) if abs(c) > 0.2), reverse=True)
            weak.append({"relative_singular_value": float(s),
                         "direction": {key_of.get(n, n): round(float(c), 3) for _, n, c in parts}})
    return {"hypothesis": e.key, "keys": keys, "values_hz": [float(x[i]) for i in idx],
            "std_hz": [float(se[i]) for i in idx],
            "correlation": [[float(corr[i, j]) for j in idx] for i in idx],
            "relative_singular_values": [float(s) for s in sv / sv[0]], "weak_directions": weak,
            "residual_scale": s2, "dof": dof,
            "note": "linearised at the optimum; valid only if the fit is in the right minimum"}


def uncertainty_markdown(u: dict, strong_correlation: float = 0.8) -> str:
    lines = [f"## Coupling uncertainties: {u['hypothesis']}", "",
             "| coupling | value (Hz) | std (Hz) | strongly correlated with |", "|---|---|---|---|"]
    keys = u["keys"]
    for i, k in enumerate(keys):
        partners = [f"{keys[j]} ({u['correlation'][i][j]:+.2f})" for j in range(len(keys))
                    if j != i and abs(u["correlation"][i][j]) >= strong_correlation]
        lines.append(f"| {k} | {u['values_hz'][i]:.3f} | {u['std_hz'][i]:.3f} | {', '.join(partners)} |")
    if u["weak_directions"]:
        lines += ["", "Weakly determined directions (relative singular value, main components):"]
        for w in u["weak_directions"]:
            lines.append(f"- {w['relative_singular_value']:.1e}: {w['direction']}")
    return "\n".join(lines)


def start_agreement(e: Evaluated, base: RefineSettings, tolerance: float = 0.01) -> dict:
    """How often the best minimum was found: the starts whose final score is within `tolerance` (relative) of the
    best, and the range of every free coupling over those starts and over all starts. Many starts in the best
    minimum with a small range: the fit is reproducible; one start alone, or a wide range among near-best starts:
    several minima fit about equally well."""
    sols = e.summary.get("start_solutions") or []
    if not sols:
        return {"starts": 0}
    settings = _settings_for(e.model, e.variant, base)
    param = settings.parameterize(e.model.interpretation)
    key_of: Dict[str, str] = {}
    for key, names in e.model.coupling_names.items():
        for n in names:
            leader = param.ties.get(n, n)
            if leader in sols[0]["parameters"] and leader not in key_of:
                key_of[leader] = key
    best = min(s["score"] for s in sols)
    near = [s for s in sols if s["score"] <= best * (1 + tolerance)]
    out = {"starts": len(sols), "near_best": len(near), "tolerance": tolerance,
           "scores_relative": sorted(round(s["score"] / best, 4) for s in sols), "couplings": {}}
    for name, key in key_of.items():
        v_near = [s["parameters"][name] for s in near]
        v_all = [s["parameters"][name] for s in sols]
        out["couplings"][key] = {"near_best_range_hz": [min(v_near), max(v_near)],
                                 "all_range_hz": [min(v_all), max(v_all)]}
    return out
