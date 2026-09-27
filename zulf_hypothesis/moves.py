"""Extension moves: turn a fragment into larger alternatives, one step at a time.

A move proposes new fragments from a fragment and the findings of its fit
(checks.py). Moves are registered by name; a search loop (planned) applies
them to the best hypotheses, refines the proposals on identical data and
keeps what an information criterion supports. Add a move by subclassing
`ExtensionMove` and calling `register_move`.

Implemented: `AddCoupledProton` (the H7 -> H8 step on lactic acid: one more
weakly coupled proton, for example a slowly exchanging OH) and
`ChangeProtonCount` (CH <-> CH2 <-> CH3 on one site and its symmetry partners;
4322bdfc: a CH-CH3 fit whose CH was the CH2 of triethylamine).
Planned: add a 15N isotopologue site, change an equivalence (CH3 <-> CH(CH3)2),
add a remote proton group, split a group.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import replace
from typing import Dict, List, Optional, Sequence

from .checks import Finding
from .fragment import Fragment, ProtonGroup, Site, pair


class ExtensionMove(ABC):
    name: str = "move"

    def triggered_by(self) -> Sequence[str]:
        """Finding codes that make this move worth trying (empty: always)."""
        return ()

    def applies(self, fragment: Fragment, findings: Sequence[Finding] = ()) -> bool:
        codes = self.triggered_by()
        return not codes or any(f.code in codes for f in findings)

    @abstractmethod
    def propose(self, fragment: Fragment, findings: Sequence[Finding] = ()) -> List[Fragment]:
        """New fragments (each with complete starting couplings)."""


MOVES: Dict[str, ExtensionMove] = {}


def register_move(move: ExtensionMove) -> ExtensionMove:
    MOVES[move.name] = move
    return move


def propose_all(fragment: Fragment, findings: Sequence[Finding] = (), only: Optional[Sequence[str]] = None,
                allowed: Optional[Sequence[str]] = None) -> List[Fragment]:
    """Proposals of every registered move whose trigger findings are present; moves named in `only`
    run regardless of triggers (an explicit request); `allowed` restricts the moves considered."""
    out = []
    for name, move in MOVES.items():
        if (only is not None and name not in only) or (allowed is not None and name not in allowed):
            continue
        if only is not None or move.applies(fragment, findings):
            out.extend(move.propose(fragment, findings))
    return out


class AddCoupledProton(ExtensionMove):
    """Attach one proton (exchangeable by default) to a site and couple it to everything.

    `host`: existing site label, or a new heteroatom site ('X' of `host_element`)
    bonded to `anchor`. Starting couplings: `initial_hz` to every included
    proton group and labelled site; symmetry is kept only when the new proton
    sits on a site every generator fixes.
    """
    name = "add_coupled_proton"

    def __init__(self, host: Optional[str] = None, anchor: Optional[str] = None, host_element: str = "O",
                 initial_hz: float = 1.0, exchangeable: bool = True, label: str = "HX"):
        self.host, self.anchor, self.host_element = host, anchor, host_element
        self.initial_hz, self.exchangeable, self.label = initial_hz, exchangeable, label

    def triggered_by(self) -> Sequence[str]:
        return ("rate_asymmetry", "misfit")

    def _anchors(self, fragment: Fragment) -> List[str]:
        if self.anchor is not None:
            return [self.anchor]
        # Default: sites carrying exactly one proton group of size 1 (a methine), where a substituent sits;
        # ring atoms are skipped (an aromatic C-H has no free valence for another substituent).
        ring = _ring_sites(fragment)
        return [s.label for s in fragment.sites
                if s.label not in ring and sum(p.size for p in fragment.protons if p.site == s.label) == 1]

    def propose(self, fragment: Fragment, findings: Sequence[Finding] = ()) -> List[Fragment]:
        out = []
        hosts = [self.host] if self.host is not None else [None] * len(self._anchors(fragment))
        for host, anchor in zip(hosts, self._anchors(fragment) if self.host is None else [self.anchor]):
            sites = list(fragment.sites)
            bonds = list(fragment.bonds)
            host_label = host
            if host_label is None:
                host_label = "X" if "X" not in fragment.labels else f"X{len(sites)}"
                sites.append(Site(host_label, self.host_element))
                if anchor is not None:
                    bonds.append((anchor, host_label))
            label = self.label if self.label not in fragment.labels else f"{self.label}{len(fragment.protons)}"
            protons = list(fragment.protons) + [ProtonGroup(label, 1, host_label, self.exchangeable)]
            couplings = dict(fragment.couplings)
            for other in [s.label for s in fragment.sites if s.label_isotopes()] + [p.label for p in fragment.protons]:
                couplings[pair(label, other)] = self.initial_hz
            symmetry = tuple(g for g in fragment.symmetry if g.get(host_label, host_label) == host_label
                             and (anchor is None or g.get(anchor, anchor) == anchor))
            out.append(replace(fragment, name=f"{fragment.name} + {label} on {host_label}", sites=tuple(sites),
                               protons=tuple(protons), couplings=couplings, symmetry=symmetry, bonds=tuple(bonds)))
        return out


# Position of the main zero-field line of an isolated X-Hn group in units of 1J(X,H): XH at J, XH2 at 3/2 J,
# XH3 at J (and 2 J). A proton-count change rescales 1J so that this line stays where the parent put it.
MAIN_LINE_FACTOR = {1: 1.0, 2: 1.5, 3: 1.0}


class ChangeProtonCount(ExtensionMove):
    """One more or one fewer proton on a site (between 1 and `max_protons`), applied to the site's whole symmetry
    orbit so the fragment keeps its symmetry. The one-bond coupling is rescaled so the group's main line stays in
    place (`MAIN_LINE_FACTOR`: a CH fitted at 204 Hz becomes a CH2 with 1J 136 Hz, not 204 Hz; e66a4b08); the
    other couplings keep their values, so a warm start carries the parent's small couplings. Exchangeable groups
    and ring sites are left alone. Not a nested model: acceptance still asks for the score improvement."""
    name = "change_proton_count"

    def __init__(self, max_protons: int = 3, steps: Sequence[int] = (1, -1)):
        self.max_protons, self.steps = max_protons, tuple(steps)

    def triggered_by(self) -> Sequence[str]:
        return ("misfit", "abundance", "collapsed_component", "rate_asymmetry")

    def propose(self, fragment: Fragment, findings: Sequence[Finding] = ()) -> List[Fragment]:
        ring = _ring_sites(fragment)
        out, done = [], set()
        for group in fragment.protons:
            if group.exchangeable or group.site in ring or group.label in done:
                continue
            orbit = {group.label} | {g.get(group.label, group.label) for g in fragment.symmetry}
            done |= orbit
            for step in self.steps:
                size = group.size + step
                if not 1 <= size <= self.max_protons:
                    continue
                protons = tuple(replace(p, size=size) if p.label in orbit else p for p in fragment.protons)
                scale = MAIN_LINE_FACTOR.get(group.size, 1.0) / MAIN_LINE_FACTOR.get(size, 1.0)
                couplings = dict(fragment.couplings)
                for p in fragment.protons:
                    key = pair(p.site, p.label)
                    if p.label in orbit and couplings.get(key) is not None:
                        couplings[key] = couplings[key] * scale
                tag = {1: "H", 2: "H2", 3: "H3"}.get(size, f"H{size}")
                out.append(replace(fragment, name=f"{fragment.name} with {group.site}{tag}", protons=protons,
                                   couplings=couplings))
        return out


def _ring_sites(fragment: Fragment) -> set:
    """Sites on a cycle of the bond graph (an edge is in a cycle if its ends stay connected without it)."""
    bonds = [tuple(b) for b in fragment.bonds]
    out = set()
    for a, b in bonds:
        adj = {}
        for x, y in bonds:
            if {x, y} == {a, b}:
                continue
            adj.setdefault(x, set()).add(y)
            adj.setdefault(y, set()).add(x)
        seen, frontier = {a}, [a]
        while frontier:
            frontier = [y for x in frontier for y in adj.get(x, ()) if y not in seen]
            seen.update(frontier)
        if b in seen:
            out.update((a, b))
    return out


class FreeRemoteCouplings(ExtensionMove):
    """Give every coupling the structure leaves unspecified (built as 0 Hz and held there: 4J, 5J, couplings
    through a heteroatom) a small starting value, so the refinement fits it. e66a4b08: freeing them lowered chi2
    from 35769 to 23353 for the same skeleton. Symmetry orbits are kept (one value per orbit)."""
    name = "free_remote_couplings"

    def __init__(self, initial_hz: float = 0.3, max_new: int = 24):
        self.initial_hz, self.max_new = initial_hz, max_new

    def triggered_by(self) -> Sequence[str]:
        return ("misfit",)

    def propose(self, fragment: Fragment, findings: Sequence[Finding] = ()) -> List[Fragment]:
        groups = [p.label for p in fragment.protons]
        orbits = fragment.coupling_orbits()               # pair -> canonical representative of its orbit
        # labelled sites that already couple to a proton (a site without any would add a new isotopologue)
        labelled = [s.label for s in fragment.sites if s.label_isotopes()
                    and any(fragment.coupling(s.label, g, orbits) is not None for g in groups)]
        members: Dict[object, list] = {}
        for k, rep in orbits.items():
            members.setdefault(rep, []).append(k)
        chosen, seen = {}, set()
        for x in labelled + groups:
            for g in groups:
                if x == g or fragment.coupling(x, g, orbits) is not None:
                    continue
                rep = orbits.get(pair(x, g), pair(x, g))
                if rep in seen:
                    continue
                seen.add(rep)
                for k in members.get(rep, [pair(x, g)]):
                    chosen[k] = self.initial_hz
        if not chosen or len(seen) > self.max_new:
            return []
        couplings = dict(fragment.couplings)
        couplings.update(chosen)
        return [replace(fragment, name=f"{fragment.name} + remote J", couplings=couplings)]


class ModelVariantMove(ABC):
    """A move on a built model that keeps the structure and changes how it is fitted (line shape, ...)."""
    name: str = "model_variant"

    def triggered_by(self) -> Sequence[str]:
        return ()

    def applies(self, model, findings: Sequence[Finding] = ()) -> bool:
        codes = self.triggered_by()
        return not codes or any(f.code in codes for f in findings)

    @abstractmethod
    def propose_model(self, model, findings: Sequence[Finding] = ()) -> list:
        """New HypothesisModel objects (different name, same structure)."""


MODEL_MOVES: Dict[str, ModelVariantMove] = {}


def register_model_move(move: ModelVariantMove) -> ModelVariantMove:
    MODEL_MOVES[move.name] = move
    return move


class GaussianLineShape(ModelVariantMove):
    """Voigt lines: one Gaussian width shared by all components (field inhomogeneity), fitted with the rates.
    Triggered by misfit; the width range comes from the data (no instrument prior)."""
    name = "gaussian_line_shape"

    def triggered_by(self) -> Sequence[str]:
        return ("misfit",)

    def propose_model(self, model, findings: Sequence[Finding] = ()) -> list:
        if model.line_shape.get("gaussian"):
            return []
        from dataclasses import replace as dc_replace
        return [dc_replace(model, name=f"{model.name} + Gaussian width",
                           line_shape=dict(model.line_shape, gaussian=True))]


def propose_model_moves(model, findings: Sequence[Finding] = (), only: Optional[Sequence[str]] = None,
                        allowed: Optional[Sequence[str]] = None) -> list:
    out = []
    for name, move in MODEL_MOVES.items():
        if (only is not None and name not in only) or (allowed is not None and name not in allowed):
            continue
        if only is not None or move.applies(model, findings):
            out.extend(move.propose_model(model, findings))
    return out


register_move(AddCoupledProton())
register_move(ChangeProtonCount())
register_move(FreeRemoteCouplings())
register_model_move(GaussianLineShape())
