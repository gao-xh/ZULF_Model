"""Bounded multistart refinement of candidate interpretations.

Every refined result is a conditional numerical candidate: it reports its
parameter boundary hits, convergence status, budgets, residuals per band and
(when supplied) frozen predictions on held-out spectra. Ranking uses held-out
prediction when available; otherwise it is marked provisional.
"""
from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import least_squares

from ..physics.protocol import SUDDEN_DROP, Protocol
from ..spinsystem import Interpretation
from ..timing import Timer
from .forward import MixtureForward, Prediction
from .observed import ObservedSpectrum
from .parameterization import Parameterization, ParameterPolicy
from .search import SearchSettings, global_search


class _BudgetReached(Exception):
    pass


@dataclass(frozen=True)
class RefineSettings:
    starts: int = 4
    start_spread_hz: float = 1.0
    max_nfev: int = 200
    max_evaluations: int = 6000     # shared by the direct and continuation paths (D28)
    max_seconds: float = 120.0
    seed: int = 0
    gain_model: str = "shared_phase"
    background_order: int = -1
    band_weighting: str = "equal"
    diff_step: float = 1e-6
    continuation_rates_per_s: tuple = (10.0, 3.0, 1.0, 0.0)
    guard_continuation: bool = True   # also fit directly at full resolution from each start (D28)
    ties: tuple = ()      # ((leader, follower, ...), ...) parameter names; missing names are skipped
    fixed: tuple = ()     # parameter names held at their candidate values
    policy: ParameterPolicy = field(default_factory=ParameterPolicy)
    search: Optional[dict] = None   # SearchSettings fields; when set, a global pattern search supplies the starts
    sign_variants: bool = False     # also refine candidates with large couplings sign-flipped (D25)
    sign_variant_min_hz: Optional[float] = None   # default: |J| whose interval excludes zero
    max_sign_variants: int = 16

    def parameterize(self, candidate: Interpretation) -> Parameterization:
        param = Parameterization.from_interpretation(candidate, self.policy)
        for group in self.ties:
            names = [n for n in group if n in param.parameters]
            if len(names) >= 2:
                param.tie(names[0], *names[1:])
        param.fix(*[n for n in self.fixed if n in param.parameters])
        return param

    @classmethod
    def from_dict(cls, data: dict) -> "RefineSettings":
        data = dict(data)
        policy = ParameterPolicy.from_dict(data.pop("policy", {}))
        data = {k: (tuple(v) if isinstance(v, list) else v) for k, v in data.items() if k in cls.__dataclass_fields__}
        if "ties" in data:
            data["ties"] = tuple(tuple(g) for g in data["ties"])
        return cls(policy=policy, **data)


@dataclass
class RefinementResult:
    interpretation: Interpretation
    parameters: Dict[str, float]
    gains: List[complex]
    background: List[complex]
    relative_residual: float
    signal_relative_residual: float
    band_relative_residuals: List[float]
    boundary_hits: List[str]
    flags: List[str]
    converged: bool
    budget_exhausted: bool
    evaluations: int
    elapsed_s: float
    attempts: List[dict]
    prediction: np.ndarray
    component_spectra: List[np.ndarray]
    validation: List[dict] = field(default_factory=list)
    candidate_index: int = -1
    search: Optional[dict] = None
    note: str = ("Conditional numerical refinement; not an assignment. Inspect boundary hits, residuals, "
                 "component spectra and held-out prediction.")

    @property
    def validation_residual(self) -> Optional[float]:
        if not self.validation:
            return None
        return float(np.mean([v["relative_residual"] for v in self.validation]))

    def summary(self) -> dict:
        return {"candidate_index": self.candidate_index, "relative_residual": self.relative_residual,
                "signal_relative_residual": self.signal_relative_residual,
                "band_relative_residuals": self.band_relative_residuals,
                "validation_relative_residual": self.validation_residual, "flags": self.flags,
                "boundary_hits": self.boundary_hits, "converged": self.converged,
                "budget_exhausted": self.budget_exhausted, "evaluations": self.evaluations,
                "elapsed_s": self.elapsed_s, "parameters": self.parameters,
                "gains": [[g.real, g.imag] for g in self.gains],
                "background": [[b.real, b.imag] for b in self.background],
                "interpretation": self.interpretation.to_dict(), "validation": self.validation,
                "search": self.search, "note": self.note}


def _band_residuals(forward: MixtureForward, prediction: Prediction) -> List[float]:
    out = []
    for b in np.unique(forward.band):
        m = forward.band == b
        out.append(float(np.linalg.norm(prediction.model[m] - forward.y[m]) / max(np.linalg.norm(forward.y[m]), 1e-30)))
    return out


def refine(candidate: Interpretation, observed: ObservedSpectrum, settings: RefineSettings = RefineSettings(),
           protocol: Protocol = SUDDEN_DROP, parameterization: Optional[Parameterization] = None,
           timer: Optional[Timer] = None, initial_points: Optional[Sequence[np.ndarray]] = None) -> RefinementResult:
    """Refine one candidate.

    Starts: `initial_points` (free-parameter vectors) when given; otherwise the
    candidate values plus `settings.starts - 1` random perturbations, or, when
    `settings.search` is set, the distinct starts of a global pattern search.
    """
    timer = timer or Timer()
    param = parameterization or settings.parameterize(candidate)
    search_summary = None
    if initial_points is None and settings.search is not None:
        with timer.section("solver.search"):
            found = global_search(param, observed, SearchSettings.from_dict(settings.search), protocol, timer)
        initial_points = found.points
        search_summary = found.summary(param)

    def make_forward(extra: float) -> MixtureForward:
        """Matched continuation: the same extra apodization is applied to data and model."""
        data = observed
        if extra > 0:
            data = observed.with_acquisition(observed.acquisition.with_processing(
                apodization_rate_per_s=observed.acquisition.apodization_rate_per_s + extra))
        return MixtureForward(param, data, protocol, settings.gain_model, settings.background_order,
                              settings.band_weighting, timer)

    forward = make_forward(0.0)
    base_forward = forward
    x0 = param.vector()
    lower, upper = param.bounds()
    rng = np.random.default_rng(settings.seed)
    start_time = time.perf_counter()
    evaluations = 0
    attempts = []

    final_best = {"score": math.inf}

    def evaluate(x):
        nonlocal evaluations
        if evaluations >= settings.max_evaluations:
            raise _BudgetReached("max_evaluations")
        if evaluations and time.perf_counter() - start_time > settings.max_seconds:
            raise _BudgetReached("max_seconds")
        evaluations += 1
        with timer.section("solver.evaluate"):
            pred = forward.predict(x)
        if forward is base_forward and pred.score < final_best["score"]:
            final_best.update(score=pred.score, x=np.array(x, float), start=current_start)
        return pred.residual

    budget = False
    current_start = 0
    history = []
    schedule = tuple(settings.continuation_rates_per_s) or (0.0,)
    if schedule[-1] != 0.0:
        schedule = schedule + (0.0,)
    continuation_note = None
    if len(schedule) > 1 and not observed.reprocessable:
        schedule = (0.0,)
        continuation_note = "continuation_unavailable_without_fid"
    level_forwards = {extra: (base_forward if extra == 0.0 else make_forward(extra)) for extra in schedule}
    explicit = [np.asarray(p, float) for p in initial_points] if initial_points is not None else None
    if explicit is not None and any(p.shape != x0.shape for p in explicit):
        raise ValueError("initial_points must be vectors over the free parameters.")
    n_starts = len(explicit) if explicit is not None else settings.starts
    for start in range(n_starts):
        current_start = start
        x = x0.copy() if explicit is None else explicit[start].copy()
        if explicit is None and start and len(x):
            couplings = np.array([param.parameters[n].kind == "coupling" for n in param.free_names])
            x = x + np.where(couplings, rng.normal(0, settings.start_spread_hz, len(x)), rng.normal(0, 0.3, len(x)))
        x = np.clip(x, lower + 1e-9 * (upper - lower), upper - 1e-9 * (upper - lower))
        if not len(x):
            evaluate(x)
            attempts.append({"start": start, "status": "no_free_parameters", "evaluations": 1})
            break
        # Guarded continuation (D28): a direct fit at full resolution runs first from the same start, then the
        # broadening schedule; the best full-resolution score of either path wins. Broadening can move the
        # optimum of a crowded spectrum away from a start that is already inside the right basin.
        paths = [schedule]
        if settings.guard_continuation and len(schedule) > 1:
            paths = [(0.0,), schedule]
        x_start = x.copy()
        for path_index, path in enumerate(paths):
            x = x_start.copy()
            for level, extra in enumerate(path):
                forward = level_forwards[extra]
                final_level = level == len(path) - 1
                before = evaluations
                try:
                    sol = least_squares(evaluate, x, bounds=(lower, upper), x_scale="jac", max_nfev=settings.max_nfev,
                                        diff_step=settings.diff_step, ftol=1e-10, xtol=1e-10, gtol=1e-10)
                except _BudgetReached as exc:
                    attempts.append({"start": start, "path": path_index, "extra_rate_per_s": extra,
                                     "status": "budget_exhausted", "reason": str(exc),
                                     "evaluations": evaluations - before})
                    budget = True
                    break
                x = np.clip(sol.x, lower, upper)
                attempts.append({"start": start, "path": path_index, "extra_rate_per_s": extra,
                                 "status": "converged" if sol.success else "stopped", "message": str(sol.message),
                                 "evaluations": evaluations - before, "score": float(sol.cost * 2)})
                if final_level:
                    history.append((start, sol.success, np.array(sol.x), float(sol.cost * 2)))
            if budget:
                break
        if budget:
            break
    forward = base_forward
    best = final_best if math.isfinite(final_best["score"]) else {"x": x0, "start": 0}
    x_best = best.get("x", x0)
    converged = any(s == best.get("start") and ok and np.allclose(xs, x_best, rtol=1e-6, atol=1e-9)
                    for s, ok, xs, _ in history)
    final = forward.predict(x_best)
    values = param.values(x_best)
    hits = param.boundary_hits(x_best)
    flags = []
    if hits:
        flags.append("search_boundary")
    if not converged:
        flags.append("optimizer_not_converged")
    if budget:
        flags.append("search_budget_exhausted")
    if continuation_note:
        flags.append(continuation_note)
    contributions = np.abs(final.gains)
    top = contributions.max() if len(contributions) and contributions.max() > 0 else 1.0
    interp = param.interpretation(values, contributions / top)
    interp.metadata.update(source="solver", flags=flags)
    relative = float(np.linalg.norm(final.model - forward.y) / max(np.linalg.norm(forward.y), 1e-30))
    signal_relative = float(np.linalg.norm(final.model - forward.y) / max(forward.background_only_residual(values), 1e-30))
    return RefinementResult(interp, values, final.gains.tolist(), final.background.tolist(), relative, signal_relative,
                            _band_residuals(forward, final), hits,
                            flags, converged, budget, evaluations, time.perf_counter() - start_time, attempts,
                            final.model, final.component_spectra, search=search_summary)


def frozen_prediction(result: RefinementResult, candidate_param: Parameterization, held_out: ObservedSpectrum,
                      settings: RefineSettings = RefineSettings(), protocol: Protocol = SUDDEN_DROP) -> dict:
    """Evaluate a refined result on a held-out spectrum with every parameter frozen.

    A conditional diagnostic refits only the linear gains; it is reported
    separately and never replaces the frozen score.
    """
    forward = MixtureForward(candidate_param, held_out, protocol, settings.gain_model, settings.background_order,
                             settings.band_weighting)
    frozen = forward.predict(values=result.parameters, fixed_gains=np.array(result.gains),
                             fixed_background=np.array(result.background, complex))
    refit = forward.predict(values=result.parameters)
    norm = max(np.linalg.norm(forward.y), 1e-30)
    return {"label": held_out.label,
            "relative_residual": float(np.linalg.norm(frozen.model - forward.y) / norm),
            "gain_refit_relative_residual": float(np.linalg.norm(refit.model - forward.y) / norm),
            "note": "Frozen score uses training gains and parameters; gain refit is a conditional diagnostic."}


def refine_candidates(candidates: Sequence[Interpretation], observed: ObservedSpectrum,
                      held_out: Sequence[ObservedSpectrum] = (), settings: RefineSettings = RefineSettings(),
                      protocol: Protocol = SUDDEN_DROP, timer: Optional[Timer] = None) -> List[RefinementResult]:
    """Refine every candidate independently and keep all branches.

    Ranking: mean frozen held-out residual when held-out spectra are supplied,
    otherwise training residual (flagged `provisional_ranking`). With
    `settings.sign_variants`, every candidate is accompanied by variants whose
    large couplings have the opposite sign (`solver.variants`); each variant
    keeps the index of its source candidate and records the flip.
    """
    from .variants import crossing_threshold_hz, sign_variants
    expanded = []
    for index, candidate in enumerate(candidates):
        expanded.append((index, candidate))
        if settings.sign_variants:
            threshold = (settings.sign_variant_min_hz if settings.sign_variant_min_hz is not None
                         else crossing_threshold_hz(settings.policy))
            expanded += [(index, v) for v in sign_variants(candidate, threshold, settings.max_sign_variants)]
    results = []
    for index, candidate in expanded:
        param = settings.parameterize(candidate)
        try:
            result = refine(candidate, observed, settings, protocol, param, timer)
        except ValueError as exc:
            results.append(None)
            continue
        result.candidate_index = index
        if "sign_variant" in candidate.metadata:
            result.interpretation.metadata["sign_variant"] = candidate.metadata["sign_variant"]
        param_frozen = Parameterization.from_interpretation(result.interpretation, settings.policy)
        result.validation = [frozen_prediction(result, param_frozen, h, settings, protocol) for h in held_out]
        if not held_out:
            result.flags.append("provisional_ranking")
        results.append(result)
    kept = [r for r in results if r is not None]
    key = (lambda r: r.validation_residual) if held_out else (lambda r: r.relative_residual)
    return sorted(kept, key=key)
