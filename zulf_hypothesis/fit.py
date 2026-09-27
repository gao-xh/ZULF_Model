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
from dataclasses import replace
from types import SimpleNamespace
from typing import List, Optional, Sequence, Tuple

from zulf_core.solver import ParameterPolicy, RefineSettings

from .builder import HypothesisModel, build_model
from .fragment import Fragment
from .labeling import Labeling
from .search import SearchResult, SearchSettings, search_hypotheses

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


def fit_settings(workers: int = 4, **changes) -> SearchSettings:
    """Search settings for a given structure: no motif scan, thorough round-0 starts, one extension round."""
    base = SearchSettings(base=default_fit_base(), top_models=100, top_motifs=0, motif_screen=0, rounds=1,
                          extend_top=1, workers=workers, initial_global_search=True)
    return dataclasses.replace(base, **changes)


def fit_structure(fragment: Fragment, observed, settings: Optional[SearchSettings] = None,
                  labeling: Optional[Labeling] = None, exchange: str = "auto", knowledge=None,
                  report: Optional[str] = None) -> SearchResult:
    """Fit a structure to an observed spectrum with every applicable variant; returns the ranked SearchResult.
    `report`: path prefix for the automatic report (report.write_report)."""
    settings = settings or fit_settings()
    labeling = labeling or Labeling.natural()
    # components without a line in the fitted ranges are omitted (e.g. a 15N isotopologue with only small
    # couplings to carbon-bound protons: low-frequency lines only)
    ranges = [tuple(r) for r in (getattr(observed, "metadata", None) or {}).get("ranges", [])] or None
    models: List[HypothesisModel] = [build_model(f, include_exchangeable=True, labeling=labeling, name=f.name,
                                                 ranges=ranges)
                                     for f in exchange_variants(fragment, exchange)]
    # extensions are built with the same labelling scheme as the structure
    proposals = SimpleNamespace(models=models, build_options={"labeling": labeling, "ranges": ranges})
    result = search_hypotheses(observed, proposals,
                               dataclasses.replace(settings, top_models=max(settings.top_models, len(models))), knowledge)
    if report:
        from .report import write_report
        write_report(result, observed, report)
    return result
