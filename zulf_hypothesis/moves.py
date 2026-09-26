"""Extension moves: turn a fragment into larger alternatives, one step at a time.

A move proposes new fragments from a fragment and the findings of its fit
(checks.py). Moves are registered by name; a search loop (planned) applies
them to the best hypotheses, refines the proposals on identical data and
keeps what an information criterion supports. Add a move by subclassing
`ExtensionMove` and calling `register_move`.

Implemented: `AddCoupledProton` (the H7 -> H8 step on lactic acid: one more
weakly coupled proton, for example a slowly exchanging OH).
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


def propose_all(fragment: Fragment, findings: Sequence[Finding] = (), only: Optional[Sequence[str]] = None) -> List[Fragment]:
    """Proposals of every registered move whose trigger findings are present; moves named in `only`
    run regardless of triggers (an explicit request)."""
    out = []
    for name, move in MOVES.items():
        if only is not None and name not in only:
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


register_move(AddCoupledProton())
