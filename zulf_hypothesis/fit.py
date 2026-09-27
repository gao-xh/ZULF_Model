"""Fit a given structure without hand-picked variants.

`fit_structure(fragment, observed)` is a thin entry point over `search_hypotheses` (the solver is unchanged): it
builds the isotopologue set of the structure in every exchange regime that applies (fast N-H / O-H exchange: those
protons dropped; slow: kept), refines each in the amplitude variants (free, fixed, ratios) from several starts
plus the global pattern search, applies the checks and the triggered extension moves, and ranks everything on
one yardstick. The work that used to be proposed by hand for every sample (e66a4b08: exchange regime, per-
isotopologue rates, more starts) is done here by default.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field, replace
from types import SimpleNamespace
from typing import List, Optional, Sequence, Tuple

from zulf_core.solver import ParameterPolicy, RefineSettings

from .builder import HypothesisModel, build_model
from .fragment import Fragment
from .labeling import Labeling
from .search import SearchResult, SearchSettings, search_hypotheses, warm_fragment

# Protons on these elements exchange with the solvent; their regime (fast or slow) is not known in advance.
EXCHANGE_ELEMENTS = ("N", "O", "S")


def default_fit_base() -> RefineSettings:
    """Refinement settings used for the confirmed samples (complex data, signal weighting, linear background)."""
    return RefineSettings(starts=4, start_spread_hz=0.5, band_weighting="signal", background_order=1,
                          max_seconds=900.0, max_evaluations=20000, continuation_rates_per_s=(0.0,),
                          policy=ParameterPolicy(coupling_margin_hz=8.0, coupling_margin_relative=0.05,
                                                 rate_bounds_per_s=(0.05, 20.0), initial_rate_per_s=2.0))


def exchangeable_groups(fragment: Fragment) -> List[str]:
    """Proton groups flagged exchangeable or sitting on N, O or S."""
    elements = {s.label: s.element for s in fragment.sites}
    return [p.label for p in fragment.protons if p.exchangeable or elements.get(p.site) in EXCHANGE_ELEMENTS]


def without_groups(fragment: Fragment, labels: Sequence[str], name: Optional[str] = None) -> Fragment:
    """The fragment with the given proton groups (and their couplings and symmetry images) removed."""
    drop = set(labels)
    protons = tuple(p for p in fragment.protons if p.label not in drop)
    couplings = {k: v for k, v in fragment.couplings.items() if not (set(k) & drop)}
    symmetry = tuple({a: b for a, b in g.items() if a not in drop and b not in drop} for g in fragment.symmetry)
    return replace(fragment, name=name or fragment.name, protons=protons, couplings=couplings, symmetry=symmetry)


def exchange_variants(fragment: Fragment, mode: str = "auto") -> List[Fragment]:
    """Fragments for the exchange regimes to fit. mode: "auto" (both when the fragment has exchangeable protons),
    "fast", "slow"."""
    groups = exchangeable_groups(fragment)
    if not groups:
        return [fragment]
    fast = without_groups(fragment, groups, f"{fragment.name} (fast exchange)")
    slow = replace(fragment, name=f"{fragment.name} (slow exchange)",
                   protons=tuple(replace(p, exchangeable=False) if p.label in groups else p for p in fragment.protons))
    return {"auto": [fast, slow], "fast": [fast], "slow": [slow]}[mode]


# Moves that refine how a given structure is fitted without changing it (proton counts, added protons would).
REFINING_MOVES = ("free_remote_couplings", "gaussian_line_shape")


class SlowExchange:
    """Model move of one fit_structure call: the structure with its exchangeable protons kept (slow exchange),
    warm-started from the fitted fast-exchange couplings. Slow-exchange models have more spins and minima; from
    generic starts they were start-sensitive (e66a4b08: chi2 25640 with 8 starts, 37976 with 4)."""
    name = "slow_exchange"

    def __init__(self, slow: Fragment, build_options: dict):
        self.slow, self.options = slow, dict(build_options)

    def triggered_by(self):
        return ("misfit", "rate_asymmetry")

    def applies(self, model, findings=()) -> bool:
        return model.fragment is not None and "(fast exchange)" in model.fragment.name and \
            any(f.code in self.triggered_by() for f in findings)

    def propose_model(self, model, findings=()) -> list:
        couplings = dict(self.slow.couplings)
        couplings.update(model.fragment.couplings)               # fitted values of the shared couplings
        name = model.fragment.name.replace("(fast exchange)", "(slow exchange)")
        fragment = replace(self.slow, name=name, couplings=couplings)
        try:
            built = build_model(fragment, include_exchangeable=True, name=model.name.replace(
                "(fast exchange)", "(slow exchange)"), **self.options)
        except ValueError:
            return []
        return [dataclasses.replace(built, line_shape=dict(model.line_shape))]


def fit_settings(workers: int = 4, **changes) -> SearchSettings:
    """Search settings for a given structure: no motif scan, thorough round-0 starts (4 perturbed starts plus a
    60 s global pattern search), one round of the structure-preserving moves in the best variant. Budget: about
    10 min on 4 cores for a 3-isotopologue amine (e66a4b08: 29 min with 8 starts, 240 s searches and all moves)."""
    base = SearchSettings(base=default_fit_base(), top_models=100, top_motifs=0, motif_screen=0, rounds=2,
                          extend_top=1, workers=workers, initial_global_search=True, moves=REFINING_MOVES,
                          extension_starts=4, extension_global_search={"max_seconds": 60.0, "solutions": 3,
                                                                       "popsize": 10, "maxiter": 40,
                                                                       "one_bond_min_hz": 50.0})
    return dataclasses.replace(base, **changes)


@dataclass
class StructureFit:
    """Result of fit_structure: the complex-data fit, and optionally the fit of the phase-corrected real part
    (phase taken from the best complex fit)."""
    result: SearchResult
    phased: Optional[SearchResult] = None
    phased_observed: object = None
    phasing: dict = field(default_factory=dict)

    @property
    def best(self):
        return self.result.best

    def table(self):
        return self.result.table()


def phased_observation(observed, phase0_rad: float, delay_s: float):
    """The phase-corrected real (absorption) part of a complex observation as a real-only ObservedSpectrum; the
    solver applies the same correction to its models (ObservedSpectrum.phasing)."""
    from zulf_core.render.phasing import phase_correct
    from zulf_core.solver import ObservedSpectrum
    sel = observed.selected
    f = observed.frequencies_hz[sel]
    real = phase_correct(observed.values[sel], f, phase0_rad, delay_s, observed.acquisition).real
    ranges = [tuple(r) for r in (observed.metadata or {}).get("ranges", [])] or None
    out = ObservedSpectrum.from_spectrum(f, real, ranges, record=observed.acquisition,
                                         phasing={"phase0_rad": phase0_rad, "delay_s": delay_s})
    out.metadata["zero_fill"] = (observed.metadata or {}).get("zero_fill", 1)
    return out


def fit_structure(fragment: Fragment, observed, settings: Optional[SearchSettings] = None,
                  labeling: Optional[Labeling] = None, exchange: str = "auto", knowledge=None,
                  report: Optional[str] = None, route: str = "complex") -> StructureFit:
    """Fit a structure to an observed spectrum with every applicable variant.

    exchange: "auto" (default) fits fast exchange first and tries slow exchange as a warm-started extension
    (accepted only if it scores better); "both" fits both regimes from the start; "fast" / "slow" one regime.

    route: "complex" (default: complex data, phase and delay fitted with the couplings), or "both": in addition,
    the spectrum is phase-corrected with the best complex fit's phase and the real part alone is refitted from
    the complex optima (delay held; a residual zero-order phase stays free and is reported). A model-free phase is
    not used: for J-spectra with dense multiplets it is not accurate enough (D43).
    `report`: path prefix for the automatic report (report.write_report; "_phased" for the second route)."""
    settings = settings or fit_settings()
    labeling = labeling or Labeling.natural()
    # components without a line in the fitted ranges are omitted (e.g. a 15N isotopologue with only small
    # couplings to carbon-bound protons: low-frequency lines only)
    ranges = [tuple(r) for r in (getattr(observed, "metadata", None) or {}).get("ranges", [])] or None
    variants = exchange_variants(fragment, "auto" if exchange in ("auto", "both") else exchange)
    extra = ()
    if exchange == "auto" and len(variants) == 2:
        extra = (SlowExchange(variants[1], {"labeling": labeling, "ranges": ranges}),)
        variants = variants[:1]
    models: List[HypothesisModel] = [build_model(f, include_exchangeable=True, labeling=labeling, name=f.name,
                                                 ranges=ranges) for f in variants]
    settings = dataclasses.replace(settings, extra_model_moves=tuple(settings.extra_model_moves) + extra)
    # extensions are built with the same labelling scheme as the structure
    proposals = SimpleNamespace(models=models, build_options={"labeling": labeling, "ranges": ranges})
    result = search_hypotheses(observed, proposals,
                               dataclasses.replace(settings, top_models=max(settings.top_models, len(models))), knowledge)
    out = StructureFit(result)
    if report:
        from .report import write_report
        write_report(result, observed, report)
    if route == "both" and result.best is not None:
        from .report import fit_phasing
        phase0, delay = fit_phasing(result.best)
        phased = phased_observation(observed, phase0, delay)
        # warm starts: each structure's best complex optimum; delay held (the data carry the correction)
        warm_models, seen = [], set()
        for e in result.ranked():
            if e.name in seen or e.model.fragment is None or e.status == "rejected" or e.demoted(settings.demote):
                continue
            seen.add(e.name)
            m = build_model(warm_fragment(e), include_exchangeable=True, labeling=labeling, name=e.name,
                            ranges=ranges)
            warm_models.append(dataclasses.replace(m, line_shape=dict(e.model.line_shape)))
        base = dataclasses.replace(settings.base, policy=dataclasses.replace(settings.base.policy,
                                                                              fit_phase_delay=False))
        phased_settings = dataclasses.replace(settings, base=base, rounds=0, initial_global_search=False,
                                              top_models=len(warm_models))
        out.phased = search_hypotheses(phased, SimpleNamespace(models=warm_models,
                                                               build_options={"labeling": labeling, "ranges": ranges}),
                                       phased_settings, knowledge)
        out.phased_observed = phased
        out.phasing = {"phase0_rad": phase0, "delay_s": delay, "source": f"complex fit {result.best.key}"}
        if report:
            write_report(out.phased, phased, report + "_phased")
    return out
