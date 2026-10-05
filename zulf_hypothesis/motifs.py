"""Motif library and motif scan: coupled spectra of whole structural motifs.

Stage-2 group candidates assume isolated X-Hn patterns, which fails when
protons couple strongly (CH2 next to CH3, aromatic rings). A motif is a
whole fragment (with symmetry) whose spectrum is computed exactly; the scan
takes the one-bond couplings from band positions, evaluates every
combination with a closed-form linear solve (one complex amplitude, the
isotopologue ratios fixed by the labelling, and the background; no
nonlinear fit), scores it on the common
yardstick and keeps the best combinations per motif for stage 3.

Couplings beyond 1J are generic values (sp3: 2J(C,H) -4.5, 3J(C,H) 4.5,
vicinal 3J(H,H) 7; aromatic: 2J(C,H) 1, 3J(C,H) 7.5, 4J(C,H) -1.2, ortho /
meta / para J(H,H) 7.5 / 1.5 / 0.7; amines with slow N-H exchange: 1J(15N,H)
negative, 2J(N,H) -1, 3J(N,H) 1, H-C-N-H 5.5), not compound-specific values, so a
blind scan does not borrow literature couplings. Add motifs with
`register_motif`; each confirmed structure type should become one.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import product
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from zulf_core.solver import RefineSettings
from zulf_core.solver.forward import MixtureForward

from .builder import HypothesisModel, build_model
from .fragments import n_ethyl_et3n as _et3n_template
from .fragment import Fragment, ProtonGroup, Site, pair
from .groups import _unit_lines
from .inventory import Inventory
from .labeling import Labeling
from .scoring import Yardstick, criterion, free_parameter_count


@dataclass(frozen=True)
class OneBondSite:
    site: str               # representative site label (symmetry partners share its 1J)
    nucleus: str            # labelled isotope, e.g. "13C"
    n_h: int                # protons on the site
    j_range_hz: Tuple[float, float]


@dataclass
class Motif:
    name: str
    factory: Callable[[Dict[str, float]], Fragment]    # one_bond {site: J} -> fragment
    one_bond: Tuple[OneBondSite, ...]
    description: str = ""

    def fragment(self, one_bond: Dict[str, float]) -> Fragment:
        return self.factory(one_bond)


MOTIFS: Dict[str, Motif] = {}


def register_motif(motif: Motif) -> Motif:
    MOTIFS[motif.name] = motif
    return motif


# -- generic couplings ---------------------------------------------------------------------------------
SP3 = {"2JCH": -4.5, "3JCH": 4.5, "3JHH": 7.0}
# Amines (slow N-H exchange): 1J(15N,H) is negative (negative gamma of 15N); generic small couplings.
AMINE = {"2JNH": -1.0, "3JNH": 1.0, "3JHH_N": 5.5}
AROM = {"2JCH": 1.0, "3JCH": 7.5, "4JCH": -1.2, "o": 7.5, "m": 1.5, "p": 0.7}


def _chain(name, groups, bonds, sym=(), extra_sites=(), note=""):
    """Fragment from (site, element, n_h) groups and site bonds; couplings by bond distance (sp3)."""
    def factory(one_bond):
        sites = [Site(s, el) for s, el, _ in groups] + [Site(s, el) for s, el in extra_sites]
        protons = [ProtonGroup(f"H{s}", n, s) for s, _, n in groups if n > 0]
        adj = {s.label: set() for s in sites}
        for a, b in bonds:
            adj[a].add(b)
            adj[b].add(a)

        def dist(a, b):
            seen, frontier, d = {a}, [a], 0
            while frontier:
                if b in frontier:
                    return d
                frontier = [y for x in frontier for y in adj[x] if y not in seen]
                seen.update(frontier)
                d += 1
            return None
        c = {}
        for s, el, n in groups:
            if n == 0:
                continue
            for t, el2, n2 in groups:
                if n2 == 0:
                    continue
                d = dist(s, t)
                nitrogen = el == "N"
                if s == t:
                    if s in one_bond:        # symmetry partners take it from the orbit
                        c[pair(s, f"H{t}")] = -abs(one_bond[s]) if nitrogen else one_bond[s]
                elif d == 1:
                    c[pair(s, f"H{t}")] = AMINE["2JNH"] if nitrogen else SP3["2JCH"]
                elif d == 2:
                    c[pair(s, f"H{t}")] = AMINE["3JNH"] if nitrogen else SP3["3JCH"]
                if s < t and d == 1:
                    c[pair(f"H{s}", f"H{t}")] = AMINE["3JHH_N"] if "N" in (el, el2) else SP3["3JHH"]
        return Fragment(name, tuple(sites), tuple(protons), c, tuple(sym), tuple(bonds), note)
    return factory


def _ring(name, n_atoms, hetero: Optional[str], note=""):
    """Six-membered aromatic ring: C-H at every carbon, optional heteroatom at position 1."""
    labels = [f"A{i}" for i in range(1, n_atoms + 1)]
    elements = ["C"] * n_atoms
    if hetero:
        elements[0] = hetero

    def factory(one_bond):
        sites = tuple(Site(l, e) for l, e in zip(labels, elements))
        carbons = [l for l, e in zip(labels, elements) if e == "C"]
        protons = tuple(ProtonGroup(f"H{l}", 1, l) for l in carbons)
        ring_d = lambda a, b: min(abs(labels.index(a) - labels.index(b)), n_atoms - abs(labels.index(a) - labels.index(b)))
        c = {}
        for s in carbons:
            for t in carbons:
                d = ring_d(s, t)
                if d == 0:
                    c[pair(s, f"H{t}")] = one_bond[_rep(s)]
                elif d == 1:
                    c[pair(s, f"H{t}")] = AROM["2JCH"]
                elif d == 2:
                    c[pair(s, f"H{t}")] = AROM["3JCH"]
                else:
                    c[pair(s, f"H{t}")] = AROM["4JCH"]
                if s < t:
                    c[pair(f"H{s}", f"H{t}")] = {1: AROM["o"], 2: AROM["m"], 3: AROM["p"]}[d]
        rotation = {l: labels[(i + 1) % n_atoms] for i, l in enumerate(labels)}
        mirror = {l: labels[(-i) % n_atoms] for i, l in enumerate(labels)}
        for m in (rotation, mirror):
            for l in list(m):
                if l in carbons:
                    m[f"H{l}"] = f"H{m[l]}"
        sym = (mirror,) if hetero else (rotation, mirror)
        bonds = tuple((labels[i], labels[(i + 1) % n_atoms]) for i in range(n_atoms))
        return Fragment(name, sites, protons, c, sym, bonds, note)

    def _rep(s):
        i = labels.index(s)
        if not hetero:
            return labels[1] if n_atoms else s
        return labels[min(i, n_atoms - i)]
    return factory, _rep


_benzene, _ = _ring("benzene ring", 6, None, "C6H6-type ring, generic aromatic couplings")
_pyridine, _ = _ring("pyridine ring", 6, "N", "C5H5N-type ring, generic aromatic couplings")

register_motif(Motif("methyl", _chain("methyl", [("C1", "C", 3)], []), (OneBondSite("C1", "13C", 3, (115, 150)),),
                     "isolated CH3 (e.g. methanol, acetate)"))
register_motif(Motif("CH-CH3", _chain("CH-CH3", [("C1", "C", 1), ("C2", "C", 3)], [("C1", "C2")]),
                     (OneBondSite("C1", "13C", 1, (125, 175)), OneBondSite("C2", "13C", 3, (115, 150))),
                     "vicinal CH-CH3 (lactate, alanine type)"))
register_motif(Motif("ethyl", _chain("ethyl", [("C1", "C", 2), ("C2", "C", 3)], [("C1", "C2")]),
                     (OneBondSite("C1", "13C", 2, (120, 160)), OneBondSite("C2", "13C", 3, (115, 150))),
                     "CH2-CH3 (ethanol, ethyl ester type)"))
register_motif(Motif("isopropyl", _chain("isopropyl", [("C1", "C", 1), ("C2", "C", 3), ("C3", "C", 3)],
                                         [("C1", "C2"), ("C1", "C3")],
                                         ({"C2": "C3", "C3": "C2", "HC2": "HC3", "HC3": "HC2"},)),
                     (OneBondSite("C1", "13C", 1, (125, 175)), OneBondSite("C2", "13C", 3, (115, 150))),
                     "CH(CH3)2"))
register_motif(Motif("CH2(CH3)2", _chain("CH2(CH3)2", [("C1", "C", 2), ("C2", "C", 3), ("C3", "C", 3)],
                                         [("C1", "C2"), ("C1", "C3")],
                                         ({"C2": "C3", "C3": "C2", "HC2": "HC3", "HC3": "HC2"},)),
                     (OneBondSite("C1", "13C", 2, (115, 160)), OneBondSite("C2", "13C", 3, (115, 150))),
                     "propane-like CH3-CH2-CH3"))
register_motif(Motif("CH2-CH2", _chain("CH2-CH2", [("C1", "C", 2), ("C2", "C", 2)], [("C1", "C2")],
                                       ({"C1": "C2", "C2": "C1", "HC1": "HC2", "HC2": "HC1"},)),
                     (OneBondSite("C1", "13C", 2, (115, 160)),), "symmetric CH2-CH2 (succinate type)"))
register_motif(Motif("CH3-C-CH3", _chain("CH3-C-CH3", [("C1", "C", 3), ("C2", "C", 3), ("C0", "C", 0)],
                                         [("C0", "C1"), ("C0", "C2")],
                                         ({"C1": "C2", "C2": "C1", "HC1": "HC2", "HC2": "HC1"},)),
                     (OneBondSite("C1", "13C", 3, (115, 150)),), "two methyls on one non-protonated carbon"))
register_motif(Motif("benzene ring", _benzene, (OneBondSite("A2", "13C", 1, (145, 180)),),
                     "six equivalent aromatic C-H"))
register_motif(Motif("pyridine ring", _pyridine,
                     (OneBondSite("A2", "13C", 1, (160, 195)), OneBondSite("A3", "13C", 1, (145, 180)),
                      OneBondSite("A4", "13C", 1, (145, 180))), "azine ring, N at position 1"))

# Amines with slow N-H exchange (acidic or non-aqueous samples). With fast exchange the N-H protons decouple and
# the spectrum is that of the carbon motif (methyl, ethyl, isopropyl ...): the 15N isotopologue keeps only small
# couplings with lines below the usual fitted ranges, so no separate fast-exchange motifs are needed.
_ISO_SWAP = {"C2": "C3", "C3": "C2", "HC2": "HC3", "HC3": "HC2"}
register_motif(Motif("CH3-NH3+", _chain("CH3-NH3+", [("C1", "C", 3), ("N1", "N", 3)], [("C1", "N1")]),
                     (OneBondSite("C1", "13C", 3, (125, 155)), OneBondSite("N1", "15N", 3, (65, 90))),
                     "methylammonium, slow N-H exchange"))
register_motif(Motif("CH3-NH2", _chain("CH3-NH2", [("C1", "C", 3), ("N1", "N", 2)], [("C1", "N1")]),
                     (OneBondSite("C1", "13C", 3, (125, 150)), OneBondSite("N1", "15N", 2, (55, 80))),
                     "methylamine, slow N-H exchange"))
register_motif(Motif("CH3CH2-NH3+", _chain("CH3CH2-NH3+", [("C1", "C", 2), ("C2", "C", 3), ("N1", "N", 3)],
                                           [("C1", "C2"), ("C1", "N1")]),
                     (OneBondSite("C1", "13C", 2, (130, 155)), OneBondSite("C2", "13C", 3, (115, 140)),
                      OneBondSite("N1", "15N", 3, (65, 90))), "ethylammonium, slow N-H exchange"))
register_motif(Motif("(CH3)2CH-NH3+", _chain("(CH3)2CH-NH3+", [("C1", "C", 1), ("C2", "C", 3), ("C3", "C", 3),
                                                              ("N1", "N", 3)],
                                             [("C1", "C2"), ("C1", "C3"), ("C1", "N1")], (_ISO_SWAP,)),
                     (OneBondSite("C1", "13C", 1, (130, 160)), OneBondSite("C2", "13C", 3, (115, 140)),
                      OneBondSite("N1", "15N", 3, (65, 90))), "isopropylammonium, slow N-H exchange"))
register_motif(Motif("H2N-CH2-CH2-NH2", _chain("H2N-CH2-CH2-NH2", [("C1", "C", 2), ("C2", "C", 2), ("N1", "N", 2),
                                                                    ("N2", "N", 2)],
                                                [("C1", "C2"), ("C1", "N1"), ("C2", "N2")],
                                                ({"C1": "C2", "C2": "C1", "N1": "N2", "N2": "N1", "HC1": "HC2",
                                                  "HC2": "HC1", "HN1": "HN2", "HN2": "HN1"},)),
                     (OneBondSite("C1", "13C", 2, (125, 150)), OneBondSite("N1", "15N", 2, (55, 80))),
                     "ethylenediamine, slow N-H exchange (fast exchange: use CH2-CH2)"))
def _xch2ch2x(one_bond):
    """X-CH2-CH2-X with fast X-H exchange, the CH2 protons kept individually (AA'BB'): the two vicinal couplings
    differ, J between protons on the same side (HC1a-HC2a, HC1b-HC2b) and J' across (HC1a-HC2b, HC1b-HC2a), so
    the protons of one CH2 are not magnetically equivalent and their geminal coupling matters. Symmetry: the two
    halves (C1 <-> C2) and the two sides (a <-> b)."""
    sites = (Site("C1", "C"), Site("C2", "C"), Site("N1", "N"), Site("N2", "N"))
    protons = tuple(ProtonGroup(f"HC{i}{x}", 1, f"C{i}") for i in (1, 2) for x in "ab")
    j1 = one_bond.get("C1", 132.0)
    c = {}
    for i, k in ((1, 2), (2, 1)):
        for x in "ab":
            c[pair(f"C{i}", f"HC{i}{x}")] = j1
            c[pair(f"C{i}", f"HC{k}{x}")] = SP3["2JCH"]
    c[pair("HC1a", "HC2a")] = c[pair("HC1b", "HC2b")] = 4.0          # J (same side)
    c[pair("HC1a", "HC2b")] = c[pair("HC1b", "HC2a")] = 8.7          # J' (across)
    c[pair("HC1a", "HC1b")] = c[pair("HC2a", "HC2b")] = -12.0        # geminal
    sym = ({"C1": "C2", "C2": "C1", "N1": "N2", "N2": "N1", "HC1a": "HC2a", "HC2a": "HC1a", "HC1b": "HC2b",
            "HC2b": "HC1b"},
           {"HC1a": "HC1b", "HC1b": "HC1a", "HC2a": "HC2b", "HC2b": "HC2a"})
    return Fragment("X-CH2-CH2-X (AA'BB')", sites, protons, c, sym, (("C1", "C2"), ("C1", "N1"), ("C2", "N2")),
                    "AA'BB' ethylene unit, fast X-H exchange")


def _h2nch2ch2nh2_aabb(one_bond):
    """Ethylenediamine with the AA'BB' ethylene unit and the N-H protons kept (slow exchange): each NH2 is one
    group (its two protons equivalent by fast rotation and inversion), coupled to the 13C (2J, 3J), to the CH2
    protons (3J(H,C,N,H) to the near pair, 4J to the far pair) and, in the 15N isotopologue, to 15N (1J)."""
    base = _xch2ch2x(one_bond)
    sites = base.sites
    protons = base.protons + (ProtonGroup("HN1", 2, "N1"), ProtonGroup("HN2", 2, "N2"))
    c = dict(base.couplings)
    jn = -abs(one_bond.get("N1", 65.0))
    for i, k in ((1, 2), (2, 1)):
        c[pair(f"N{i}", f"HN{i}")] = jn
        c[pair(f"C{i}", f"HN{i}")] = -3.0          # 2J(C,N,H)
        c[pair(f"C{i}", f"HN{k}")] = 1.5           # 3J(C,C,N,H)
        for x in "ab":
            c[pair(f"HC{i}{x}", f"HN{i}")] = 5.5   # 3J(H,C,N,H)
            c[pair(f"HC{k}{x}", f"HN{i}")] = 0.0   # 4J
            c[pair(f"N{i}", f"HC{i}{x}")] = -1.0   # 2J(N,C,H)
            c[pair(f"N{i}", f"HC{k}{x}")] = 1.0    # 3J(N,C,C,H)
    c[pair("HN1", "HN2")] = 0.0
    sym = (dict(base.symmetry[0], HN1="HN2", HN2="HN1"), base.symmetry[1])
    return Fragment("H2N-CH2-CH2-NH2 (AA'BB')", sites, protons, c, sym, base.bonds,
                    "AA'BB' ethylene unit with the N-H protons kept (slow exchange)")


register_motif(Motif("H2N-CH2-CH2-NH2 (AA'BB')", _h2nch2ch2nh2_aabb,
                     (OneBondSite("C1", "13C", 2, (125, 150)), OneBondSite("N1", "15N", 2, (55, 90))),
                     "ethylenediamine, AA'BB' ethylene unit, slow N-H exchange (fast: use X-CH2-CH2-X (AA'BB'))"))
register_motif(Motif("X-CH2-CH2-X (AA'BB')", _xch2ch2x, (OneBondSite("C1", "13C", 2, (125, 150)),),
                     "1,2-disubstituted ethane with fast X-H exchange (ethylenediamine): J and J' separate"))
register_motif(Motif("(CH3)2CH-NH2", _chain("(CH3)2CH-NH2", [("C1", "C", 1), ("C2", "C", 3), ("C3", "C", 3),
                                                            ("N1", "N", 2)],
                                            [("C1", "C2"), ("C1", "C3"), ("C1", "N1")], (_ISO_SWAP,)),
                     (OneBondSite("C1", "13C", 1, (125, 155)), OneBondSite("C2", "13C", 3, (115, 140)),
                      OneBondSite("N1", "15N", 2, (55, 80))), "isopropylamine, slow N-H exchange"))

# Tertiary amine without N-H: the ethyl carbons see the other N-CH2 protons through N (3J(C,N,C,H)); the reduced
# fragment keeps them as one 4H group (template "N-ethyl (Et3N)"). Confirmed case: triethylamine (4322bdfc).
register_motif(Motif("N-ethyl (Et3N)", lambda one_bond: _et3n_template(one_bond),
                     (OneBondSite("C1", "13C", 2, (120, 160)), OneBondSite("C2", "13C", 3, (115, 150))),
                     "ethyl on a tertiary N with two more N-CH2 groups (triethylamine type)"))

# Secondary amine with ethyl and methyl on N, fast N-H exchange (the NH decoupled): the carbons couple to the other
# side's protons through N (3J(C,N,C,H)). Confirmed case: N-ethylmethylamine (e66a4b08).
register_motif(Motif("CH3CH2-N-CH3", _chain("CH3CH2-N-CH3", [("C1", "C", 3), ("C2", "C", 2), ("N1", "N", 0),
                                                            ("C3", "C", 3)],
                                            [("C1", "C2"), ("C2", "N1"), ("N1", "C3")]),
                     (OneBondSite("C1", "13C", 3, (115, 140)), OneBondSite("C2", "13C", 2, (120, 150)),
                      OneBondSite("C3", "13C", 3, (120, 145))),
                     "N-ethyl-N-methyl amine, fast N-H exchange (N-ethylmethylamine type)"))


@dataclass
class MotifProposal:
    motif: str
    one_bond: Dict[str, float]
    model: HypothesisModel
    chi2: float
    score: float

    def describe(self) -> str:
        j = ", ".join(f"{k}={v:.1f}" for k, v in self.one_bond.items())
        return f"motif {self.motif} [{j}]"


def _site_candidates(site: OneBondSite, inventory: Inventory, per_site: int) -> List[float]:
    lines = sorted(inventory.lines, key=lambda l: -l.snr)
    out = []
    for line in lines:
        for ratio, _ in _unit_lines(site.nucleus, site.n_h):
            j = line.frequency_hz / ratio
            if site.j_range_hz[0] <= j <= site.j_range_hz[1] and all(abs(j - x) > 1.0 for x in out):
                out.append(j)
        if len(out) >= per_site:
            break
    return out


def _quick_chi2(model: HypothesisModel, observed, stick: Yardstick, rate_per_s: float) -> float:
    # A motif is one molecule: its isotopologue ratios are fixed by the labelling, so a minor isotopologue
    # (e.g. 15N) cannot take over a band of another nucleus in the quick score.
    settings = model.settings(RefineSettings(), fixed_ratios=True, fix_unspecified=True)
    param = settings.parameterize(model.interpretation)
    values = param.values()
    for key in values:
        if ".log_rate" in key:
            values[key] = float(np.log(rate_per_s))
    forward = MixtureForward(param, observed, gain_model="complex", background=1,
                             amplitude_map=settings.amplitude_map)
    pred = forward.predict(values=values)
    return stick.chi2(observed.values[observed.selected], pred.model)


def _refine_one_bond(motif, best: "MotifProposal", observed, inventory, stick, labeling, min_ratio, rate_per_s,
                     kind, span_hz, step_hz, passes) -> "MotifProposal":
    current = best
    offsets = np.arange(-span_hz, span_hz + step_hz / 2, step_hz)
    for _ in range(passes):
        improved = False
        for site in motif.one_bond:
            center = current.one_bond[site.site]
            for d in offsets:
                if d == 0:
                    continue
                one_bond = dict(current.one_bond, **{site.site: float(center + d)})
                lo, hi = site.j_range_hz
                if not lo <= one_bond[site.site] <= hi:
                    continue
                try:
                    model = build_model(motif.fragment(one_bond), ranges=inventory.ranges, labeling=labeling,
                                        min_ratio=min_ratio, name=f"motif {motif.name}")
                    chi2 = _quick_chi2(model, observed, stick, rate_per_s)
                except ValueError:
                    continue
                k = free_parameter_count(model, model.settings(RefineSettings()))
                score = criterion(chi2, k, stick.n, kind)
                if score < current.score:
                    current = MotifProposal(motif.name, one_bond, model, chi2, score)
                    improved = True
        if not improved:
            break
    return current


def scan_motifs(observed, inventory: Inventory, stick: Yardstick, motifs: Optional[Sequence[str]] = None,
                labeling: Optional[Labeling] = None, min_ratio: float = 0.0, per_site: int = 6,
                max_combinations: int = 300, keep_per_motif: int = 1, rate_per_s: float = 2.0,
                kind: str = "bic", refine_hz: float = 4.0, refine_step_hz: float = 0.5,
                refine_passes: int = 2, refine_top: int = 8) -> List[MotifProposal]:
    """Best one-bond assignments of every motif, scored on the yardstick (linear solve only).

    The band positions give 1J to a few Hz only (triethylamine and N-ethylmethylamine: the right motif lost
    against CH-CH3 at band-position values). The best combination of each motif is therefore refined by
    coordinate search: each site's 1J on a grid of +-`refine_hz` (step `refine_step_hz`), `refine_passes`
    passes, symmetric sites through the motif's own representatives, for the `refine_top` leading motifs
    (0 disables)."""
    out: List[MotifProposal] = []
    for name in (motifs or list(MOTIFS)):
        motif = MOTIFS[name]
        choices = [_site_candidates(s, inventory, per_site) for s in motif.one_bond]
        if any(not c for c in choices):
            continue
        scored = []
        for combo in list(product(*choices))[:max_combinations]:
            one_bond = {s.site: float(j) for s, j in zip(motif.one_bond, combo)}
            try:
                model = build_model(motif.fragment(one_bond), ranges=inventory.ranges, labeling=labeling,
                                    min_ratio=min_ratio, name=f"motif {name}")
                chi2 = _quick_chi2(model, observed, stick, rate_per_s)
            except ValueError:
                continue
            k = free_parameter_count(model, model.settings(RefineSettings()))
            scored.append(MotifProposal(name, one_bond, model, chi2, criterion(chi2, k, stick.n, kind)))
        scored.sort(key=lambda p: p.score)
        out.extend(scored[:keep_per_motif])
    out.sort(key=lambda p: p.score)
    if refine_hz > 0:
        # only the leading motifs (those the search screens) are worth the local 1J search
        for i, p in enumerate(out[:refine_top]):
            out[i] = _refine_one_bond(MOTIFS[p.motif], p, observed, inventory, stick, labeling, min_ratio,
                                      rate_per_s, kind, refine_hz, refine_step_hz, refine_passes)
        out.sort(key=lambda p: p.score)
    for p in out:
        p.model.name = p.describe()
    return out
