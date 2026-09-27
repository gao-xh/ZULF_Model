"""Stage 3: budgeted search over hypotheses.

    proposals (stage 2) [+ hinted interpretations]
      -> refine each in variants ("free": free amplitudes and rates;
         "fixed": abundance ratios and one shared rate), in parallel processes
      -> common yardstick (scoring.py): chi2 on the data cores, BIC
      -> checks (checks.py) per result, with the other variant as peer
      -> extension round(s): moves (moves.py) triggered by the findings of the
         best clean hypotheses, built like the originals and refined on the
         same data; accepted only if the BIC improves by `accept_delta`
      -> knowledge matching of the leaders (knowledge.py)
      -> ranked table and a step log

The solver is used through `refine` / `RefineSettings` only. Hint sources
never enter the ranking; hinted interpretations are refined as ordinary
candidates and labelled by origin.
"""
from __future__ import annotations

import dataclasses
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np

from zulf_core.solver import RefineSettings, refine

from .builder import HypothesisModel, build_model, model_from_interpretation
from .checks import Finding, run_checks
from .knowledge import KnowledgeBase, knowledge_matches
from .moves import propose_all
from .scoring import Yardstick, criterion, free_parameter_count, yardstick


def default_base() -> RefineSettings:
    return RefineSettings(starts=2, start_spread_hz=0.5, band_weighting="signal", background_order=1,
                          max_seconds=300.0, max_evaluations=20000, continuation_rates_per_s=(0.0,))


@dataclass
class SearchSettings:
    base: RefineSettings = field(default_factory=default_base)
    variants: Tuple[str, ...] = ("free", "fixed")
    top_models: int = 4               # stage-2 group proposals refined in round 0
    top_motifs: int = 3               # motif-scan proposals refined in round 0
    include_hinted: int = 0           # hinted interpretations refined in round 0 (free variant only)
    rounds: int = 1                   # extension rounds
    extend_top: int = 2               # hypotheses extended per round
    extension_starts: int = 8         # starts for extensions (warm start plus perturbations): new couplings are
    extension_spread_hz: float = 1.5  # unknown and have several comparable minima
    accept_delta: float = 6.0         # (quasi-)BIC improvement needed to accept an extension
    clean_margin: float = 10.0        # prefer a hypothesis without warnings if it is within this of the minimum
    overdispersion: bool = True       # scale chi2 by the best reduced chi2 (quasi-likelihood) before ranking
    workers: int = 1                  # parallel refinement processes
    criterion: str = "bic"            # "bic" | "aic"
    knowledge_kinds: Tuple[str, ...] = ("measured",)
    knowledge_exclude: Tuple[str, ...] = ()   # e.g. the sample id: never match references fitted on this sample


@dataclass
class Evaluated:
    name: str
    variant: str
    origin: str
    model: HypothesisModel
    summary: dict
    prediction: np.ndarray = field(repr=False)
    k: int = 0
    chi2: float = float("nan")
    n: int = 0
    score: float = float("nan")
    findings: List[Finding] = field(default_factory=list)
    couplings: dict = field(default_factory=dict)
    knowledge: List[dict] = field(default_factory=list)
    parent: Optional[str] = None
    status: str = "evaluated"         # evaluated | accepted | rejected

    @property
    def clean(self) -> bool:
        return not any(f.severity in ("warn", "reject") for f in self.findings)

    @property
    def key(self) -> str:
        return f"{self.name} [{self.variant}]"

    def row(self, best: float) -> dict:
        amp = [abs(complex(*g)) for g in self.summary["gains"]]
        top = max(amp) if amp and max(amp) > 0 else 1.0
        return {"hypothesis": self.name, "variant": self.variant, "origin": self.origin, "status": self.status,
                "score": round(self.score, 2), "delta": round(self.score - best, 2), "chi2": round(self.chi2, 1),
                "n": self.n, "k": self.k, "data_region_residual": self.summary.get("data_region_residual"),
                "amplitudes": [round(a / top, 3) for a in amp], "expected": [round(r, 3) for r in self.model.ratios]
                if self.model.abundance_known else None,
                "bounds": self.summary.get("boundary_hits"), "findings": [f.code for f in self.findings],
                "couplings": {k: round(v, 2) for k, v in self.couplings.items()}, "knowledge": self.knowledge[:2]}


@dataclass
class SearchResult:
    evaluated: List[Evaluated]
    log: List[dict]
    yardstick: Yardstick = field(repr=False)
    c_hat: float = 1.0                # overdispersion used in the scores
    clean_margin: float = 10.0

    def ranked(self) -> List[Evaluated]:
        return sorted(self.evaluated, key=lambda e: e.score)

    @property
    def best(self) -> Optional[Evaluated]:
        """Lowest score; a hypothesis without warnings is preferred only when it is within `clean_margin`."""
        ranked = self.ranked()
        if not ranked:
            return None
        clean = [e for e in ranked if e.clean and e.score <= ranked[0].score + self.clean_margin]
        return clean[0] if clean else ranked[0]

    def table(self) -> List[dict]:
        ranked = self.ranked()
        best = ranked[0].score if ranked else 0.0
        return [dict(rank=i + 1, **e.row(best)) for i, e in enumerate(ranked)]


def _settings_for(model: HypothesisModel, variant: str, base: RefineSettings) -> RefineSettings:
    if not model.abundance_known:
        return base
    fixed = variant == "fixed"
    return model.settings(base, fixed_ratios=fixed, shared_rate=fixed)


def _refine_job(job):
    model, variant, observed, base = job[:4]
    overrides = dict(job[4]) if len(job) > 4 and job[4] else {}
    alternates = overrides.pop("alternates", None)
    settings = _settings_for(model, variant, base)
    if overrides:
        settings = dataclasses.replace(settings, **overrides)
    start = time.perf_counter()
    initial = None
    if alternates:
        # Several starting models of the same structure (e.g. warm start from the parent and the cold proposal
        # values), each also perturbed; refine keeps the best.
        param = settings.parameterize(model.interpretation)
        lo, hi = param.bounds()
        rng = np.random.default_rng(settings.seed)
        coupling = np.array([param.parameters[n].kind == "coupling" for n in param.free_names])
        initial = []
        for alt in [model] + list(alternates):
            x = np.clip(settings.parameterize(alt.interpretation).vector(), lo, hi)
            initial.append(x)
            for _ in range(max(settings.starts - 1, 0) // (1 + len(alternates))):
                initial.append(np.clip(x + np.where(coupling, rng.normal(0, settings.start_spread_hz, len(x)), 0.0),
                                       lo, hi))
    try:
        result = refine(model.interpretation, observed, settings, initial_points=initial)
    except Exception as exc:              # a failed candidate is reported, not fatal
        return {"error": f"{type(exc).__name__}: {exc}"}
    return {"summary": result.summary(), "prediction": np.asarray(result.prediction),
            "k": free_parameter_count(model, settings), "seconds": time.perf_counter() - start}


class _Runner:
    def __init__(self, observed, settings: SearchSettings, stick: Yardstick):
        self.observed, self.settings, self.stick = observed, settings, stick
        sel = observed.selected
        self.values = observed.values[sel]

    def run(self, jobs: Sequence[tuple], log: List[dict], overrides: Optional[dict] = None) -> List[Evaluated]:
        payload = [(j[0], j[1], self.observed, self.settings.base,
                    dict(overrides or {}, **({"alternates": j[4]} if len(j) > 4 and j[4] else {})))
                   for j in jobs]
        if self.settings.workers > 1 and len(payload) > 1:
            with ProcessPoolExecutor(max_workers=self.settings.workers) as pool:
                outputs = list(pool.map(_refine_job, payload))
        else:
            outputs = [_refine_job(p) for p in payload]
        out = []
        for job, res in zip(jobs, outputs):
            model, variant, origin, parent = job[:4]
            if "error" in res:
                log.append({"step": "refine", "hypothesis": model.name, "variant": variant, "error": res["error"]})
                continue
            chi2 = self.stick.chi2(self.values, res["prediction"])
            score = criterion(chi2, res["k"], self.stick.n, self.settings.criterion)    # rescored later
            couplings = model.named_couplings(res["summary"]["parameters"])
            res["summary"]["reduced_chi2"] = chi2 / max(self.stick.n - res["k"], 1)
            out.append(Evaluated(model.name, variant, origin, model, res["summary"], res["prediction"], res["k"],
                                 chi2, self.stick.n, score, couplings=couplings, parent=parent))
            log.append({"step": "refine", "hypothesis": model.name, "variant": variant, "score": round(score, 2),
                        "seconds": round(res["seconds"], 1)})
        return out


def _fixed_differs(model: HypothesisModel) -> bool:
    """The fixed variant constrains something only if some ratio block has two or more components."""
    return model.abundance_known and any(len(b) > 1 for b in model.blocks())


def _rescore(evaluated: List[Evaluated], settings: SearchSettings, n: int) -> float:
    """Quasi-likelihood scores: chi2 / c_hat + penalty, c_hat = smallest reduced chi2 (>= 1). Real spectra are
    never fitted to the noise, and the unmodelled structure would otherwise inflate every difference."""
    c_hat = 1.0
    if settings.overdispersion and evaluated:
        c_hat = max(1.0, min(e.chi2 / max(n - e.k, 1) for e in evaluated))
    for e in evaluated:
        e.score = criterion(e.chi2 / c_hat, e.k, n, settings.criterion)
    return c_hat


def _decide_extensions(evaluated: List[Evaluated], settings: SearchSettings, log: List[dict]) -> None:
    by_key = {e.key: e for e in evaluated}
    for e in evaluated:
        if e.parent is None:
            continue
        parent = by_key.get(e.parent)
        better = parent is not None and e.score < parent.score - settings.accept_delta
        e.status = "accepted" if better else "rejected"
        log.append({"step": "accept" if better else "reject", "hypothesis": e.key, "parent": e.parent,
                    "delta": round(e.score - parent.score, 2) if parent else None})


def _same_model(a: HypothesisModel, b: HypothesisModel) -> bool:
    ca, cb = a.interpretation.components, b.interpretation.components
    return len(ca) == len(cb) and all(
        x.system.isotopes == y.system.isotopes and x.system.groups == y.system.groups
        and np.allclose(x.system.couplings_hz, y.system.couplings_hz) for x, y in zip(ca, cb))


def _check(evaluated: List[Evaluated]) -> None:
    for e in evaluated:
        peers = [p.summary for p in evaluated if p.name == e.name and p is not e]
        e.findings = run_checks(e.model, e.summary, peers)


def search_hypotheses(observed, proposals, settings: Optional[SearchSettings] = None,
                      knowledge: Optional[KnowledgeBase] = None) -> SearchResult:
    """Refine, check, extend and rank. `proposals`: a ProposalSet (stage 2) or a list of HypothesisModel."""
    settings = settings or SearchSettings()
    models = list(getattr(proposals, "models", proposals))[:settings.top_models]
    build_options = dict(getattr(proposals, "build_options", {}) or {})
    hinted = list(getattr(proposals, "hinted_interpretations", []))[:settings.include_hinted]
    stick = yardstick(observed)
    runner = _Runner(observed, settings, stick)
    log: List[dict] = [{"step": "yardstick", "points": int(stick.mask.sum()), "n": stick.n, "sigma": stick.sigma}]

    motif_models = list(getattr(proposals, "motif_models", []))[:settings.top_motifs]
    jobs = [(m, v, f"proposal {i + 1}", None) for i, m in enumerate(models) for v in settings.variants
            if v == "free" or _fixed_differs(m)]
    for i, m in enumerate(motif_models):
        twin = next((o for o in models if _same_model(m, o)), None)
        if twin is not None:          # the same structure from both routes is refined once
            log.append({"step": "dedupe", "motif": m.name, "same_as": twin.name})
            continue
        jobs += [(m, v, f"motif {i + 1}", None) for v in settings.variants if v == "free" or _fixed_differs(m)]
    jobs += [(model_from_interpretation(interp, f"hint {i + 1}"), "free", "hint", None)
             for i, interp in enumerate(hinted)]
    evaluated = runner.run(jobs, log)
    _rescore(evaluated, settings, stick.n)
    _check(evaluated)

    for round_ in range(settings.rounds):
        ranked = sorted(evaluated, key=lambda e: e.score)
        parents, seen = [], set()
        for e in sorted(ranked, key=lambda e: (not e.clean, e.score)):
            if e.model.fragment is None or e.model.parts or e.name in seen:
                continue
            seen.add(e.name)
            parents.append(e)
            if len(parents) >= settings.extend_top:
                break
        new_jobs = []
        for parent in parents:
            findings = [f for e in evaluated if e.name == parent.name for f in e.findings]
            # Warm start: the parent's refined couplings replace the proposal's starting values, so an
            # extension (a nested model) starts next to the parent's optimum instead of from generic values.
            refined = {k: v for k, v in parent.couplings.items() if k.startswith("J(") and ":" not in k}
            warm = parent.model.fragment.with_couplings(refined) if refined else parent.model.fragment
            cold_fragments = propose_all(parent.model.fragment, findings)
            for n_move, fragment in enumerate(propose_all(warm, findings)):
                # A move hypothesises what it adds (a coupled proton is in slow exchange by definition), so
                # extensions keep exchangeable protons whatever the proposals were built with.
                options = dict(build_options, include_exchangeable=True)
                try:
                    model = build_model(fragment, **options)
                except ValueError as exc:
                    log.append({"step": "extend", "from": parent.name, "proposal": fragment.name, "error": str(exc)})
                    continue
                if _same_model(model, parent.model):
                    log.append({"step": "extend", "from": parent.name, "proposal": fragment.name,
                                "note": "adds nothing observable; skipped"})
                    continue
                # The warm start inherits the parent's compensations (e.g. broadened lines standing in for the
                # missing coupling), so the cold start from the proposal values is refined as well.
                alternates = []
                if n_move < len(cold_fragments):
                    try:
                        alternates = [build_model(cold_fragments[n_move], **options)]
                    except ValueError:
                        alternates = []
                for v in settings.variants:
                    if v == "free" or _fixed_differs(model):
                        new_jobs.append((model, v, f"extension of {parent.name}", f"{parent.name} [{v}]",
                                         alternates))
                log.append({"step": "extend", "round": round_ + 1, "from": parent.name, "proposal": fragment.name,
                            "triggers": sorted({f.code for f in findings})})
        if not new_jobs:
            log.append({"step": "extend", "round": round_ + 1, "note": "no move triggered"})
            break
        added = runner.run(new_jobs, log, {"starts": settings.extension_starts,
                                           "start_spread_hz": settings.extension_spread_hz})
        evaluated.extend(added)
        _rescore(evaluated, settings, stick.n)
        _check(evaluated)

    c_hat = _rescore(evaluated, settings, stick.n)          # final scale; extensions decided on it
    log.append({"step": "score", "c_hat": round(c_hat, 3), "criterion": settings.criterion})
    _decide_extensions(evaluated, settings, log)

    kb = knowledge or KnowledgeBase.packaged()
    for e in sorted(evaluated, key=lambda e: e.score)[:max(3, settings.top_models)]:
        e.knowledge = knowledge_matches(e.model, e.couplings, kb, kinds=settings.knowledge_kinds,
                                        exclude=settings.knowledge_exclude)
    return SearchResult(evaluated, log, stick, c_hat, settings.clean_margin)
