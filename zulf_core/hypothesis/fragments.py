"""Fragment templates used so far, with generic starting couplings.

Starting values are generic sp3 values (1J(C,H) about 125-150 Hz, 2J(C,H)
about -4.5 Hz, vicinal 3J(H,H) about 7 Hz), not molecule-specific literature
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
                 pair("Cb", "Ha"): -4.5, pair("Ha", "Hb"): 7.0}
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
                 pair("Cb1", "Ha"): -4.5, pair("Cb1", "Hb2"): 4.5, pair("Ha", "Hb1"): 7.0}
    swap = {"Cb1": "Cb2", "Cb2": "Cb1", "Hb1": "Hb2", "Hb2": "Hb1"}
    return Fragment("CH(CH3)2", sites, protons, couplings, (swap,), (("Ca", "Cb1"), ("Ca", "Cb2")),
                    "isopropyl; the 4J between the two methyls is left unspecified (built as 0 Hz, fixed)")
