"""Synthetic recognition benchmark for stage 2 (and optionally stage 3).

Each case is a truth fragment with its own couplings (different from the
generic motif values); its spectrum is simulated at natural abundance with
noise, run through `propose_hypotheses`, and the rank of the first proposal
with the same topology is reported for the group-based list and the motif
list. Add cases with `register_case`; run `python -m zulf_hypothesis.benchmark`.
"""
from __future__ import annotations

import time
from typing import Callable, Dict, List, Optional

import numpy as np

from zulf_core.physics import compute_transitions
from zulf_core.render.acquisition import Acquisition
from zulf_core.render.renderer import Renderer
from zulf_core.solver import ObservedSpectrum

from .builder import build_model
from .fragment import Fragment, ProtonGroup, Site, pair
from .knowledge import match_template

CASES: Dict[str, Callable[[], Fragment]] = {}


def register_case(name: str):
    def wrap(fn):
        CASES[name] = fn
        return fn
    return wrap


def same_topology(a: Optional[Fragment], b: Optional[Fragment]) -> bool:
    if a is None or b is None:
        return False
    sizes = lambda f: sorted((f.site(p.site).element, p.size) for p in f.protons)
    return sizes(a) == sizes(b) and match_template(a, b) is not None and match_template(b, a) is not None


def _fragments_of(model) -> List[Fragment]:
    return [m.fragment for m in model.parts] if model.parts else [model.fragment]


def observe_truth(fragment: Fragment, rate: float = 1.0, snr: float = 60.0, seed: int = 0,
                  ranges=((20.0, 400.0),)) -> ObservedSpectrum:
    acq = Acquisition(1000.0, 6000, start_sample=40, sg_window=101, sg_order=2)
    model = build_model(fragment)
    renderer = Renderer(acq)
    fid = sum(renderer.synthesize(compute_transitions(c.system), rate, gain=r)
              for c, r in zip(model.interpretation.components, model.ratios))
    spectrum_peak = np.abs(np.fft.rfft(fid[40:])).max()
    noise = spectrum_peak / snr / np.sqrt(len(fid) - 40)
    fid = fid + np.random.default_rng(seed).normal(0, noise, len(fid))
    return ObservedSpectrum.from_fid(fid, acq, list(ranges))


def run_benchmark(cases: Optional[List[str]] = None, snr: float = 60.0, search: bool = False, search_settings=None,
                  verbose: bool = True) -> List[dict]:
    from .pipeline import propose_hypotheses
    rows = []
    for name in (cases or list(CASES)):
        truth = CASES[name]()
        t0 = time.perf_counter()
        obs = observe_truth(truth, snr=snr)
        ps = propose_hypotheses(obs)
        group_rank = next((i + 1 for i, m in enumerate(ps.models)
                           if len(_fragments_of(m)) == 1 and same_topology(_fragments_of(m)[0], truth)), None)
        motif_rank = next((i + 1 for i, p in enumerate(ps.motif_proposals) if same_topology(p.model.fragment, truth)),
                          None)
        row = {"case": name, "bands": len(ps.inventory.bands), "group_rank": group_rank, "motif_rank": motif_rank,
               "top_group": ps.proposals[0].describe() if ps.proposals else None,
               "top_motif": ps.motif_proposals[0].describe() if ps.motif_proposals else None,
               "seconds": round(time.perf_counter() - t0, 1)}
        if search:
            from .search import search_hypotheses
            res = search_hypotheses(obs, ps, search_settings)
            best = res.best
            row["search_best"] = best.key if best else None
            row["search_correct"] = bool(best and any(same_topology(f, truth) for f in _fragments_of(best.model)))
        rows.append(row)
        if verbose:
            print(row, flush=True)
    return rows


# -- cases (couplings chosen near typical values, deliberately not equal to the generic motif values) ------------
S, P = Site, ProtonGroup


@register_case("methyl")
def _methyl():
    return Fragment("methyl", (S("C1", "C"),), (P("H1", 3, "C1"),), {pair("C1", "H1"): 141.0})


@register_case("CH-CH3")
def _ch_ch3():
    return Fragment("CH-CH3", (S("C1", "C"), S("C2", "C")), (P("H1", 1, "C1"), P("H2", 3, "C2")),
                    {pair("C1", "H1"): 147.1, pair("C2", "H2"): 129.2, pair("C1", "H2"): -4.4, pair("C2", "H1"): -4.1,
                     pair("H1", "H2"): 7.0}, (), (("C1", "C2"),))


@register_case("isopropyl")
def _isopropyl():
    swap = {"C2": "C3", "C3": "C2", "H2": "H3", "H3": "H2"}
    return Fragment("isopropyl", (S("C1", "C"), S("C2", "C"), S("C3", "C")),
                    (P("H1", 1, "C1"), P("H2", 3, "C2"), P("H3", 3, "C3")),
                    {pair("C1", "H1"): 140.0, pair("C2", "H2"): 126.0, pair("C1", "H2"): -4.0, pair("C2", "H1"): -3.8,
                     pair("C2", "H3"): 4.6, pair("H1", "H2"): 6.6}, (swap,), (("C1", "C2"), ("C1", "C3")))


@register_case("ethyl")
def _ethyl():
    return Fragment("ethyl", (S("C1", "C"), S("C2", "C")), (P("H1", 2, "C1"), P("H2", 3, "C2")),
                    {pair("C1", "H1"): 141.0, pair("C2", "H2"): 125.5, pair("C1", "H2"): -4.6, pair("C2", "H1"): -2.4,
                     pair("H1", "H2"): 7.1}, (), (("C1", "C2"),))


@register_case("CH2(CH3)2")
def _propane():
    swap = {"C2": "C3", "C3": "C2", "H2": "H3", "H3": "H2"}
    return Fragment("propane", (S("C1", "C"), S("C2", "C"), S("C3", "C")),
                    (P("H1", 2, "C1"), P("H2", 3, "C2"), P("H3", 3, "C3")),
                    {pair("C1", "H1"): 125.0, pair("C2", "H2"): 124.5, pair("C1", "H2"): -4.2, pair("C2", "H1"): -4.4,
                     pair("C2", "H3"): 4.2, pair("H1", "H2"): 7.3}, (swap,), (("C1", "C2"), ("C1", "C3")))


@register_case("CH2-CH2")
def _succinate():
    swap = {"C1": "C2", "C2": "C1", "H1": "H2", "H2": "H1"}
    return Fragment("CH2-CH2", (S("C1", "C"), S("C2", "C")), (P("H1", 2, "C1"), P("H2", 2, "C2")),
                    {pair("C1", "H1"): 130.0, pair("C1", "H2"): -4.0, pair("H1", "H2"): 7.0}, (swap,), (("C1", "C2"),))


@register_case("CH3-C-CH3")
def _acetone():
    swap = {"C1": "C2", "C2": "C1", "H1": "H2", "H2": "H1"}
    return Fragment("acetone", (S("C1", "C"), S("C2", "C"), S("C0", "C")), (P("H1", 3, "C1"), P("H2", 3, "C2")),
                    {pair("C1", "H1"): 127.0, pair("C1", "H2"): 0.6}, (swap,), (("C0", "C1"), ("C0", "C2")))


@register_case("benzene ring")
def _benzene():
    labels = [f"A{i}" for i in range(1, 7)]
    c = {}
    for i, s in enumerate(labels):
        for k, t in enumerate(labels):
            d = min(abs(i - k), 6 - abs(i - k))
            c[pair(s, f"H{t}")] = {0: 158.4, 1: 1.1, 2: 7.6, 3: -1.3}[d]
            if i < k:
                c[pair(f"H{s}", f"H{t}")] = {1: 7.5, 2: 1.4, 3: 0.7}[d]
    rot = {l: labels[(i + 1) % 6] for i, l in enumerate(labels)}
    rot.update({f"H{l}": f"H{labels[(i + 1) % 6]}" for i, l in enumerate(labels)})
    return Fragment("benzene", tuple(S(l, "C") for l in labels), tuple(P(f"H{l}", 1, l) for l in labels), c, (rot,),
                    tuple((labels[i], labels[(i + 1) % 6]) for i in range(6)))


@register_case("pyridine ring")
def _pyridine():
    CH = {2: {2: 178.0, 3: 3.1, 4: 6.8, 5: -1.4, 6: 11.2}, 3: {3: 162.5, 2: 3.1, 4: 0.8, 5: 6.5, 6: -1.0},
          4: {4: 162.0, 3: 0.9, 5: 0.9, 2: 6.8, 6: 6.8}}
    HH = {(2, 3): 4.9, (2, 4): 1.8, (2, 5): 0.9, (2, 6): -0.1, (3, 4): 7.7, (3, 5): 1.4}
    c = {pair(f"A{ci}", f"HA{hi}"): v for ci, row in CH.items() for hi, v in row.items()}
    c.update({pair(f"HA{a}", f"HA{b}"): v for (a, b), v in HH.items()})
    mirror = {"A2": "A6", "A6": "A2", "A3": "A5", "A5": "A3", "HA2": "HA6", "HA6": "HA2", "HA3": "HA5", "HA5": "HA3"}
    sites = (S("A1", "N"),) + tuple(S(f"A{i}", "C") for i in range(2, 7))
    return Fragment("pyridine", sites, tuple(P(f"HA{i}", 1, f"A{i}") for i in range(2, 7)), c, (mirror,),
                    tuple((f"A{i}", f"A{i % 6 + 1}") for i in range(1, 7)))


@register_case("CH3-15NH3")
def _methylammonium():
    return Fragment("methylammonium", (S("C1", "C"), S("N1", "N")), (P("H1", 3, "C1"), P("H2", 3, "N1")),
                    {pair("C1", "H1"): 145.0, pair("N1", "H2"): -75.0, pair("C1", "H2"): -2.0, pair("N1", "H1"): -1.5,
                     pair("H1", "H2"): 6.5}, (), (("C1", "N1"),))


@register_case("(CH3)2CH-NH3+")
def _isopropylammonium():
    swap = {"C2": "C3", "C3": "C2", "H2": "H3", "H3": "H2"}
    return Fragment("isopropylammonium", (S("C1", "C"), S("C2", "C"), S("C3", "C"), S("N1", "N")),
                    (P("H1", 1, "C1"), P("H2", 3, "C2"), P("H3", 3, "C3"), P("HN", 3, "N1")),
                    {pair("C1", "H1"): 143.0, pair("C2", "H2"): 127.5, pair("C1", "H2"): -4.2, pair("C2", "H1"): -4.0,
                     pair("C2", "H3"): 4.4, pair("H1", "H2"): 6.6, pair("N1", "HN"): -74.5, pair("N1", "H1"): -1.4,
                     pair("N1", "H2"): 1.8, pair("C1", "HN"): -1.5, pair("C2", "HN"): 1.2, pair("H1", "HN"): 6.2},
                    (swap,), (("C1", "C2"), ("C1", "C3"), ("C1", "N1")))


@register_case("CH3-15NH2")
def _methylamine():
    return Fragment("methylamine", (S("C1", "C"), S("N1", "N")), (P("H1", 3, "C1"), P("H2", 2, "N1")),
                    {pair("C1", "H1"): 132.0, pair("N1", "H2"): -64.5, pair("C1", "H2"): -2.2, pair("N1", "H1"): -0.8,
                     pair("H1", "H2"): 5.8}, (), (("C1", "N1"),))


@register_case("ethylenediamine (fast exchange)")
def _ethylenediamine():
    swap = {"C1": "C2", "C2": "C1", "H1": "H2", "H2": "H1", "N1": "N2", "N2": "N1"}
    return Fragment("ethylenediamine", (S("C1", "C"), S("C2", "C"), S("N1", "N"), S("N2", "N")),
                    (P("H1", 2, "C1"), P("H2", 2, "C2")),
                    {pair("C1", "H1"): 133.0, pair("C1", "H2"): -3.5, pair("H1", "H2"): 6.0, pair("N1", "H1"): -1.2,
                     pair("N1", "H2"): 1.5}, (swap,), (("C1", "C2"), ("C1", "N1"), ("C2", "N2")))


if __name__ == "__main__":
    run_benchmark()
