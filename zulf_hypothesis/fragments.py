"""Fragment templates used so far, with generic starting couplings.

Starting values are generic sp3 values (1J(C,H) about 125-150 Hz, 2J(C,H)
about -4.5 Hz, vicinal 3J(H,H) about 7 Hz, 1J(C,C) about 35 Hz for doubly
labelled isotopologues), not molecule-specific literature
couplings; pass measured line positions through `one_bond` to start near the
data. Add templates with `register_template`.
"""
from __future__ import annotations

from typing import Callable, Dict, Optional

from .fragment import Fragment, ProtonGroup, Site, pair

TEMPLATES: Dict[str, Callable[..., Fragment]] = {}


def register_template(name: str):
    def wrap(fn):
        TEMPLATES[name] = fn
        return fn
    return wrap


def template(name: str, **kwargs) -> Fragment:
    return TEMPLATES[name](**kwargs)


@register_template("CH-CH3")
def ch_ch3(one_bond: Optional[Dict[str, float]] = None, substituent: Optional[str] = None,
           heavy_neighbour: Optional[str] = None) -> Fragment:
    """R-CH(X)-CH3. `substituent`: element of X on the CH ('N' adds a site that can carry 15N).
    `heavy_neighbour`: element of R (e.g. 'C' for a carboxyl carbon, labelled 13C)."""
    ob = {"Ca": 146.0, "Cb": 130.0, **(one_bond or {})}
    sites = [Site("Ca", "C"), Site("Cb", "C")]
    bonds = [("Ca", "Cb")]
    protons = [ProtonGroup("Ha", 1, "Ca"), ProtonGroup("Hb", 3, "Cb")]
    couplings = {pair("Ca", "Ha"): ob["Ca"], pair("Cb", "Hb"): ob["Cb"], pair("Ca", "Hb"): -4.5,
                 pair("Cb", "Ha"): -4.5, pair("Ha", "Hb"): 7.0, pair("Ca", "Cb"): 35.0}
    if substituent:
        sites.append(Site("X", substituent))
        bonds.append(("Ca", "X"))
    if heavy_neighbour:
        sites.append(Site("Cr", heavy_neighbour))
        bonds.append(("Ca", "Cr"))
    return Fragment("CH-CH3", tuple(sites), tuple(protons), couplings, (), tuple(bonds),
                    "vicinal CH-CH3 (alanine, lactic acid type)")


@register_template("CH(CH3)2")
def isopropyl(one_bond: Optional[Dict[str, float]] = None) -> Fragment:
    """(CH3)2CH-R with two equivalent methyls."""
    ob = {"Ca": 146.0, "Cb": 130.0, **(one_bond or {})}
    sites = (Site("Ca", "C"), Site("Cb1", "C"), Site("Cb2", "C"))
    protons = (ProtonGroup("Ha", 1, "Ca"), ProtonGroup("Hb1", 3, "Cb1"), ProtonGroup("Hb2", 3, "Cb2"))
    couplings = {pair("Ca", "Ha"): ob["Ca"], pair("Cb1", "Hb1"): ob["Cb"], pair("Ca", "Hb1"): -4.5,
                 pair("Cb1", "Ha"): -4.5, pair("Cb1", "Hb2"): 4.5, pair("Ha", "Hb1"): 7.0,
                 pair("Ca", "Cb1"): 35.0}
    swap = {"Cb1": "Cb2", "Cb2": "Cb1", "Hb1": "Hb2", "Hb2": "Hb1"}
    return Fragment("CH(CH3)2", sites, protons, couplings, (swap,), (("Ca", "Cb1"), ("Ca", "Cb2")),
                    "isopropyl; the 4J between the two methyls is left unspecified (built as 0 Hz, fixed)")


@register_template("CH2-CH2")
def ch2_ch2(one_bond: Optional[Dict[str, float]] = None) -> Fragment:
    """X-CH2-CH2-X with two equivalent CH2 groups (succinate, ethylenediamine, ethylene glycol type)."""
    ob = {"C1": 132.0, **(one_bond or {})}
    swap = {"C1": "C2", "C2": "C1", "H1": "H2", "H2": "H1"}
    return Fragment("CH2-CH2", (Site("C1", "C"), Site("C2", "C")), (ProtonGroup("H1", 2, "C1"), ProtonGroup("H2", 2, "C2")),
                    {pair("C1", "H1"): ob["C1"], pair("C1", "H2"): -4.5, pair("H1", "H2"): 7.0, pair("C1", "C2"): 35.0},
                    (swap,), (("C1", "C2"),), "symmetric CH2-CH2")


@register_template("N-ethyl (Et3N)")
def n_ethyl_et3n(one_bond: Optional[Dict[str, float]] = None) -> Fragment:
    """One ethyl of N(CH2CH3)3, reduced: the four protons of the two other N-CH2 groups form one group HX on a
    never-labelled pseudo-site X, coupled to the labelled CH2 carbon through N (3J(C,N,C,H)). Their couplings to
    their own methyls and the remote methyls are dropped (the full 16-spin molecule is not tractable)."""
    ob = {"C1": 131.0, "C2": 125.0, **(one_bond or {})}
    sites = (Site("C1", "C"), Site("C2", "C"), Site("N1", "N", ()), Site("X", "C", ()))
    protons = (ProtonGroup("H1", 2, "C1"), ProtonGroup("H2", 3, "C2"), ProtonGroup("HX", 4, "X"))
    couplings = {pair("C1", "H1"): ob["C1"], pair("C2", "H2"): ob["C2"], pair("C1", "H2"): -4.5,
                 pair("C2", "H1"): -4.5, pair("H1", "H2"): 7.0, pair("C1", "HX"): 4.5}
    return Fragment("N-ethyl (Et3N)", sites, protons, couplings, (), (("C1", "C2"), ("C1", "N1"), ("N1", "X")),
                    "N-CH2-CH3 of a triethylamine-type N(CH2R)3; other N-CH2 protons as one group")


@register_template("CH3CH2-N-CH3")
def ch3ch2_n_ch3(one_bond: Optional[Dict[str, float]] = None) -> Fragment:
    """N-ethyl-N-methyl amine with the N-H decoupled (fast exchange); same labels as the motif."""
    from .motifs import MOTIFS
    return MOTIFS["CH3CH2-N-CH3"].fragment({"C1": 125.0, "C2": 135.0, "C3": 131.5, **(one_bond or {})})

