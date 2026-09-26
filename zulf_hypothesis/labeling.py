"""Isotope labelling schemes: natural abundance (default) or enrichment.

A `Labeling` gives the probability that a site carries a spin-1/2 isotope:
the natural abundance from the nucleus registry, overridden uniformly per
isotope (`isotopes={"13C": 0.99}`) or per site (`sites={"Ca": {"13C": 0.99}}`).
The builder weights every label set exactly:

    P(L) = prod_{s in L} a_s(iso_s) * prod_{s labelable, not in L} (1 - sum_iso a_s(iso))

so natural samples (P ~ product of small abundances) and enriched samples
(where the fully labelled molecule dominates and singly labelled ones are
the minor species) use the same code. `max_labels` limits the label-set
size (natural default 1 in `build_model`, 2 in the proposal pipeline;
enriched default: all labelable sites). Components below
`primary_fraction` of the strongest set are minor: they follow a parent
component's amplitude and decay rate instead of getting free ones.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional

from zulf_core.nuclei import get_registry


@dataclass(frozen=True)
class Labeling:
    mode: str = "natural"                                         # "natural" | "enriched"
    isotopes: Mapping[str, float] = field(default_factory=dict)   # uniform overrides per isotope
    sites: Mapping[str, Mapping[str, float]] = field(default_factory=dict)   # per-site overrides
    max_labels: Optional[int] = None
    primary_fraction: float = 0.05

    @classmethod
    def natural(cls, max_labels: int = 1, primary_fraction: float = 0.05) -> "Labeling":
        return cls("natural", {}, {}, max_labels, primary_fraction)

    @classmethod
    def enriched(cls, isotopes: Optional[Mapping[str, float]] = None,
                 sites: Optional[Mapping[str, Mapping[str, float]]] = None, max_labels: Optional[int] = None,
                 primary_fraction: float = 0.05) -> "Labeling":
        if not isotopes and not sites:
            raise ValueError("An enriched labelling needs isotope or site enrichment levels.")
        return cls("enriched", dict(isotopes or {}), {k: dict(v) for k, v in (sites or {}).items()},
                   max_labels, primary_fraction)

    def __post_init__(self):
        if self.mode not in ("natural", "enriched"):
            raise ValueError("Labeling mode must be 'natural' or 'enriched'.")
        for level in list(self.isotopes.values()) + [v for s in self.sites.values() for v in s.values()]:
            if not 0.0 <= float(level) <= 1.0:
                raise ValueError("Isotope fractions must lie in [0, 1].")

    def abundance(self, site: str, isotope: str) -> float:
        if site in self.sites and isotope in self.sites[site]:
            return float(self.sites[site][isotope])
        if isotope in self.isotopes:
            return float(self.isotopes[isotope])
        return float(get_registry()[isotope].natural_abundance)

    def label_limit(self, n_labelable: int) -> int:
        if self.max_labels is not None:
            return max(1, min(self.max_labels, n_labelable))
        return 1 if self.mode == "natural" else max(1, n_labelable)

    def to_dict(self) -> Dict:
        return {"mode": self.mode, "isotopes": dict(self.isotopes), "sites": {k: dict(v) for k, v in self.sites.items()},
                "max_labels": self.max_labels, "primary_fraction": self.primary_fraction}
