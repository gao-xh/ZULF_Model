"""Stage 2 pipeline: spectrum -> ranked, buildable hypotheses (no refinement yet).

    inventory -> group candidates -> hints (optional) -> fragment enumeration
    -> built isotopologue models (+ interpretations passed on by hint providers)

Stage 3 (planned) refines the models in parallel, runs the checks, applies
triggered extension moves and ranks by an information criterion. Hints
change the search order and add candidates; they never enter the ranking.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from zulf_core.spinsystem import Interpretation
from .builder import HypothesisModel, build_model, combine_models, min_ratio_for_snr
from .enumerate import FragmentProposal, enumerate_fragments
from .groups import GroupCandidate, group_candidates
from .hints import Hint, HintProvider, apply_hints, insight_report
from .inventory import Inventory, band_inventory
from .labeling import Labeling


@dataclass
class ProposalSet:
    inventory: Inventory
    candidates: List[GroupCandidate]
    hints: List[Hint]
    insight: List[dict]
    proposals: List[FragmentProposal]
    models: List[HypothesisModel]
    hinted_interpretations: List[Interpretation] = field(default_factory=list)
    failed: List[dict] = field(default_factory=list)

    def table(self) -> List[dict]:
        return [{"rank": i + 1, "proposal": p.describe(), "score": round(p.score, 3), "sources": p.sources,
                 "components": m.component_labels, "ratios": [round(r, 3) for r in m.ratios],
                 "omitted": m.omitted}
                for i, (p, m) in enumerate(zip(self.proposals, self.models))]


def propose_hypotheses(observed, instrument_hz: Sequence[float] = (), providers: Sequence[HintProvider] = (),
                       hint_weight: float = 0.1, top_per_band: int = 3, max_groups: int = 4,
                       prior: str = "generic_sp3", max_proposals: int = 20, threshold: float = 4.0,
                       include_exchangeable: bool = False, labeling: Optional[Labeling] = None,
                       min_ratio: Optional[float] = None) -> ProposalSet:
    """`labeling`: natural abundance with up to two labels per isotopologue (default) or an enriched scheme
    (`Labeling.enriched`). Minor isotopologues are built at their abundance weight and kept only when their
    lines could reach 2 sigma (`min_ratio` defaults to 2 / peak SNR of the strongest band); omitted ones are
    listed with the reason."""
    labeling = labeling or Labeling.natural(max_labels=2)
    inventory = band_inventory(observed, instrument_hz, threshold=threshold)
    candidates = group_candidates(inventory)
    hints: List[Hint] = []
    for provider in providers:
        hints.extend(provider.hints(observed, {"inventory": inventory}))
    if hints:
        candidates = apply_hints(candidates, hints, inventory, weight=hint_weight)
    insight = insight_report(inventory, candidates, hints) if hints else []
    proposals = enumerate_fragments(candidates, inventory, top_per_band, max_groups, prior,
                                    max_proposals=max_proposals)
    if min_ratio is None:
        peak = max((b.peak.snr for b in inventory.bands), default=1.0)
        min_ratio = min_ratio_for_snr(peak)
    built, kept, failed = [], [], []
    for p in proposals:
        try:
            parts = [build_model(f, include_exchangeable=include_exchangeable, ranges=inventory.ranges,
                                 labeling=labeling, min_ratio=min_ratio)
                     for f in p.fragments]
        except ValueError as exc:
            failed.append({"proposal": p.describe(), "reason": str(exc)})
            continue
        built.append(combine_models(parts, p.describe()))
        kept.append(p)
    interps = [h.payload["interpretation"] for h in hints if h.kind == "interpretation"]
    return ProposalSet(inventory, candidates, hints, insight, kept, built, interps, failed)
