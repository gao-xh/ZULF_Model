"""Hints: outside suggestions (neural models, other tools, a person) that steer
the hypothesis search without deciding it.

A `HintProvider` returns `Hint`s for an observed spectrum. Hints can
* raise the search priority of data-derived group candidates they agree with
  (`apply_hints`: a bounded bonus; nothing is removed),
* add group candidates the data scoring ranked low or missed, marked as
  hint-sourced (they are refined like any other candidate),
* carry whole interpretations, refined as extra candidates.
The final ranking always comes from refinement and checks, never from hint
confidence. `insight_report` lists where hints and data agree or disagree,
which is the useful output of a weak model.

`ProposerHints` adapts any object with `propose(observed, k) ->
List[Interpretation]` (for example `zulf_model.evaluation.ModelProposer`),
so this module needs no torch. Other sources implement `HintProvider` and
are added with `register_hint_provider`.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..spinsystem import Interpretation
from .groups import GroupCandidate, PATTERNS
from .inventory import Inventory


@dataclass
class Hint:
    kind: str                     # "group" | "interpretation" | "note"
    source: str                   # provider name (model checkpoint, tool, person)
    confidence: float             # the provider's own score in [0, 1]; used for search order only
    payload: dict = field(default_factory=dict)
    # group: {"nucleus", "n_h", "j_hz"}; interpretation: {"interpretation": Interpretation}; note: {"text"}

    def describe(self) -> str:
        if self.kind == "group":
            p = self.payload
            return f"{p['nucleus']}H{p['n_h']} J={p['j_hz']:.1f} Hz ({self.source}, {self.confidence:.2f})"
        if self.kind == "note":
            return f"{self.payload.get('text', '')} ({self.source})"
        return f"{self.kind} ({self.source}, {self.confidence:.2f})"


class HintProvider(ABC):
    name: str = "hints"

    @abstractmethod
    def hints(self, observed, context: Optional[dict] = None) -> List[Hint]:
        """Hints for the observed spectrum; `context` may carry the inventory and ranges."""


HINT_PROVIDERS: Dict[str, HintProvider] = {}


def register_hint_provider(provider: HintProvider) -> HintProvider:
    HINT_PROVIDERS[provider.name] = provider
    return provider


def group_hints_from_interpretation(interp: Interpretation, source: str, confidence: float,
                                    one_bond_min_hz: float = 40.0) -> List[Hint]:
    """Per component: every non-1H spin with its directly coupled protons (|J| >= one_bond_min_hz)."""
    out = []
    for comp in interp.components:
        s = comp.system
        j = np.asarray(s.couplings_hz)
        for i, iso in enumerate(s.isotopes):
            if iso == "1H":
                continue
            partners = [k for k, other in enumerate(s.isotopes) if other == "1H" and abs(j[i, k]) >= one_bond_min_hz]
            if not partners:
                continue
            out.append(Hint("group", source, confidence,
                            {"nucleus": iso, "n_h": len(partners), "j_hz": float(np.mean(np.abs(j[i, partners]))),
                             "component_label": comp.label}))
    return out


class ProposerHints(HintProvider):
    """Adapter for candidate proposers (neural or not): top-k interpretations as hints.

    Confidence falls with rank (1, 1/2, 1/3, ...). The interpretations are
    passed on as well, so the search can refine them as extra candidates.
    """

    def __init__(self, proposer, name: Optional[str] = None, k: int = 5):
        self.proposer = proposer
        self.name = name or getattr(proposer, "name", "proposer")
        self.k = k

    def hints(self, observed, context: Optional[dict] = None) -> List[Hint]:
        out = []
        for rank, interp in enumerate(self.proposer.propose(observed, self.k)):
            conf = 1.0 / (rank + 1)
            out.append(Hint("interpretation", self.name, conf, {"interpretation": interp, "rank": rank}))
            out.extend(group_hints_from_interpretation(interp, self.name, conf))
        return out


def apply_hints(candidates: List[GroupCandidate], hints: Sequence[Hint], inventory: Inventory,
                weight: float = 0.1, j_tol_hz: float = 4.0, add_missing: bool = True) -> List[GroupCandidate]:
    """Bonus for data candidates a group hint agrees with (at most `weight`); hint-only candidates added.

    A hint-only candidate is created when a group hint predicts a line inside
    an observed band but no data candidate of that type and J exists; it
    starts with score 0, so it ranks below data-supported alternatives unless
    no data candidate explains that band at all.
    """
    out = list(candidates)
    for h in hints:
        if h.kind != "group":
            continue
        p = h.payload
        match = [c for c in out if c.nucleus == p["nucleus"] and c.n_h == p["n_h"] and abs(c.j_hz - p["j_hz"]) <= j_tol_hz]
        for c in match:
            bonus = min(weight, c.hint_support + weight * h.confidence)
            if bonus > c.hint_support:
                c.hint_support = bonus
            if h.source not in c.sources:
                c.sources.append(h.source)
        if match or not add_missing:
            continue
        pattern = next((q for q in PATTERNS.values() if q.nucleus == p["nucleus"] and q.n_h == p["n_h"]), None)
        if pattern is None:
            continue
        lines = [(r * p["j_hz"], a, "hinted", inventory.support(r * p["j_hz"])) for r, a in pattern.unit_lines]
        bands = sorted({b.index for f, *_ in lines for b in [inventory.band_of(f, 1.5)] if b is not None})
        if not bands:
            continue
        out.append(GroupCandidate(pattern.name, pattern.nucleus, pattern.n_h, float(p["j_hz"]), bands[0], lines,
                                  tuple(bands), 0.0, pattern.prior, weight * h.confidence, [h.source],
                                  ["hint-only: no data candidate of this type and J"]))
    out.sort(key=lambda c: (c.anchor_band,) + c.sort_key())
    return out


def insight_report(inventory: Inventory, candidates: Sequence[GroupCandidate], hints: Sequence[Hint]) -> List[dict]:
    """Per band: the best data candidate, the hints that touch the band, and whether they agree."""
    from .groups import by_band
    groups = by_band([c for c in candidates if "data" in c.sources])
    rows = []
    for band in inventory.bands:
        data_best = groups.get(band.index, [None])[0]
        touching = []
        for h in hints:
            if h.kind != "group":
                continue
            pattern = next((q for q in PATTERNS.values() if q.nucleus == h.payload["nucleus"]
                            and q.n_h == h.payload["n_h"]), None)
            if pattern and any(band.contains(r * h.payload["j_hz"], 1.5) for r, _ in pattern.unit_lines):
                touching.append(h)
        agree = [h for h in touching if data_best is not None and h.payload["nucleus"] == data_best.nucleus
                 and h.payload["n_h"] == data_best.n_h and abs(h.payload["j_hz"] - data_best.j_hz) <= 4.0]
        rows.append({"band": band.index, "peak_hz": round(band.peak.frequency_hz, 2),
                     "data_best": None if data_best is None else f"{data_best.pattern} J={data_best.j_hz:.1f}",
                     "hints": [h.describe() for h in touching],
                     "verdict": ("no hint" if not touching else "agree" if agree else "disagree")})
    silent = [h.describe() for h in hints if h.kind == "group" and not any(
        b.contains(r * h.payload["j_hz"], 1.5) for b in inventory.bands
        for q in PATTERNS.values() if q.nucleus == h.payload["nucleus"] and q.n_h == h.payload["n_h"]
        for r, _ in q.unit_lines)]
    if silent:
        rows.append({"band": None, "peak_hz": None, "data_best": None, "hints": silent,
                     "verdict": "hint predicts lines where the data show none"})
    for source in sorted({h.source for h in hints}):
        mine = [h for h in hints if h.source == source]
        groups_ = [h for h in mine if h.kind == "group"]
        agree = sum(1 for r in rows if r["verdict"] == "agree" for h in groups_ if h.describe() in r["hints"])
        off = sum(1 for h in groups_ if h.describe() in silent)
        rows.append({"band": None, "peak_hz": None, "data_best": None, "hints": [],
                     "verdict": f"source {source}: {len(groups_)} group hints, {agree} agree with the data, "
                                f"{off} predict lines where the data show none"
                                + ("; no group hints (no one-bond couplings proposed)" if not groups_ else "")})
    return rows
