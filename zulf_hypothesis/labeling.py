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


# ---- labelling hypotheses for blind samples (D54) ------------------------------------------------------------
LABELING_MODES = ("natural", "15N", "2H-exchange", "unknown")


def deuterate_exchangeable(fragment, name: Optional[str] = None):
    """The fragment with its exchangeable proton groups as 2H (a sample in D2O); couplings to them are converted
    at equal reduced coupling K (J x g_2H / g_1H per deuterated partner, D51). None without exchangeable groups."""
    from dataclasses import replace
    from zulf_core.nuclei import convert_coupling
    from .fit import exchangeable_groups              # flagged, or on N / O / S (the exchange variants' rule)
    deut = set(exchangeable_groups(fragment))
    if not deut:
        return None
    couplings = {}
    for key, value in fragment.couplings.items():
        n = len(set(key) & deut)
        couplings[key] = value if n == 0 else convert_coupling(value, ("1H", "1H"),
                                                               ("2H", "1H") if n == 1 else ("2H", "2H"))
    protons = tuple(replace(p, isotope="2H") if p.label in deut else p for p in fragment.protons)
    return replace(fragment, name=name or f"{fragment.name} [2H exchange]", protons=protons, couplings=couplings)


def labeling_variants(fragment, mode: str = "unknown", enrichment: float = 0.98):
    """[(name, fragment, Labeling)] of the labelling hypotheses to fit for a sample whose labelling is `mode`:
    "natural" (natural abundance), "15N" (N sites enriched to `enrichment`), "2H-exchange" (exchangeable protons
    deuterated), "unknown" (every one that applies to this fragment: 15N needs an N site, 2H an exchangeable
    proton group)."""
    if mode not in LABELING_MODES:
        raise ValueError(f"labelling mode must be one of {LABELING_MODES}")
    out = []
    if mode in ("natural", "unknown"):
        out.append(("natural", fragment, Labeling.natural()))
    has_n = any(s.element == "N" for s in fragment.sites)
    if mode == "15N" or (mode == "unknown" and has_n):
        if not has_n:
            raise ValueError("15N labelling needs a nitrogen site in the fragment")
        out.append(("15N", fragment, Labeling.enriched(isotopes={"15N": enrichment}, max_labels=2)))
    deuterated = deuterate_exchangeable(fragment)
    if mode == "2H-exchange" or (mode == "unknown" and deuterated is not None):
        if deuterated is None:
            raise ValueError("2H exchange labelling needs exchangeable protons in the fragment")
        out.append(("2H-exchange", deuterated, Labeling.natural()))
    return out
