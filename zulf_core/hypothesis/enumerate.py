"""Fragment enumeration from group candidates.

1. Covers: the smallest sets of group candidates (top-k per band) whose lines
   explain every band. A 2J partner band is explained by its methyl, so it
   is never taken for a separate group.
2. Topologies per cover: all fragments made of isolated parts, and every
   tree over the group sites (Pruefer sequences; groups bonded directly).
3. Equivalence: a methyl leaf may be doubled into two symmetry-equivalent
   methyls on the same neighbour (isopropyl type).
4. Couplings: 1J from the candidate; longer-range couplings from a
   `CouplingPrior` by bond distance (registry `PRIORS`; generic sp3 values by
   default, unspecified couplings are built as 0 Hz and held fixed).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Dict, List, Optional, Sequence, Tuple

from .fragment import Fragment, ProtonGroup, Site, pair
from .groups import GroupCandidate, by_band
from .inventory import Inventory


class CouplingPrior(ABC):
    name = "prior"

    @abstractmethod
    def heteronuclear(self, nucleus: str, bonds: int) -> Optional[float]:
        """Coupling between a labelled nucleus and a proton `bonds` bonds away (2, 3, ...)."""

    @abstractmethod
    def homonuclear(self, bonds: int) -> Optional[float]:
        """Proton-proton coupling across `bonds` bonds (3 = vicinal)."""


class GenericSp3Prior(CouplingPrior):
    """Generic sp3 starting values; None leaves a coupling unspecified (0 Hz, fixed)."""
    name = "generic_sp3"
    HETERO = {("13C", 2): -4.5, ("13C", 3): 4.5, ("15N", 2): -1.0, ("15N", 3): 1.0}
    HOMO = {3: 7.0}

    def heteronuclear(self, nucleus: str, bonds: int) -> Optional[float]:
        return self.HETERO.get((nucleus, bonds))

    def homonuclear(self, bonds: int) -> Optional[float]:
        return self.HOMO.get(bonds)


PRIORS: Dict[str, CouplingPrior] = {}


def register_prior(prior: CouplingPrior) -> CouplingPrior:
    PRIORS[prior.name] = prior
    return prior


register_prior(GenericSp3Prior())

ELEMENT = {"13C": "C", "15N": "N"}


@dataclass
class FragmentProposal:
    fragments: List[Fragment]              # one per disconnected part
    groups: List[GroupCandidate]
    topology: str
    score: float
    sources: List[str] = field(default_factory=list)
    unexplained: List[int] = field(default_factory=list)   # weak bands this proposal leaves open

    @property
    def n_groups(self) -> int:
        return len(self.groups)

    def describe(self) -> str:
        text = f"{self.topology}  [" + ", ".join(f"{g.pattern} J={g.j_hz:.1f}" for g in self.groups) + "]"
        return text + (f"  unexplained bands {self.unexplained}" if self.unexplained else "")


def covers(candidates: Sequence[GroupCandidate], inventory: Inventory, top_per_band: int = 3, max_groups: int = 4,
           optional_snr: float = 8.0, max_covers: int = 12) -> List[Tuple[Tuple[GroupCandidate, ...], List[int]]]:
    """Minimal candidate sets explaining every strong band; weak bands (peak SNR < optional_snr) may stay open.

    Returns (cover, unexplained weak bands), ordered by size, open bands, score.
    """
    strong = {b.index for b in inventory.bands if b.peak.snr >= optional_snr}
    weak = {b.index for b in inventory.bands} - strong
    pool = []
    for band, cands in sorted(by_band(candidates).items()):
        best = {}
        for c in cands:                       # best J per group type within the band
            best.setdefault(c.pattern, c)
        pool.extend(sorted(best.values(), key=lambda c: c.sort_key())[:top_per_band])
    found = []
    for size in range(1, max_groups + 1):
        for combo in combinations(pool, size):
            explained = set().union(*(c.explains for c in combo))
            if not strong <= explained:
                continue
            if any(strong | (weak & explained) <= set().union(*(d.explains for d in combo if d is not c))
                   for c in combo):
                continue                      # a member adds nothing
            found.append((combo, sorted(weak - explained)))
        if len(found) >= max_covers:
            break
    found.sort(key=lambda t: (len(t[0]), len(t[1]), -sum(c.rank_score for c in t[0])))
    return found[:max_covers]


def _trees(n: int) -> List[List[Tuple[int, int]]]:
    """All labelled trees on n nodes (Pruefer sequences)."""
    if n == 1:
        return [[]]
    if n == 2:
        return [[(0, 1)]]
    out = []
    for seq in product(range(n), repeat=n - 2):
        degree = [1] * n
        for x in seq:
            degree[x] += 1
        edges, seq = [], list(seq)
        for x in seq:
            leaf = min(i for i in range(n) if degree[i] == 1)
            edges.append((min(leaf, x), max(leaf, x)))
            degree[leaf] -= 1
            degree[x] -= 1
        u, v = [i for i in range(n) if degree[i] == 1]
        edges.append((u, v))
        out.append(sorted(edges))
    return out


def _distances(n: int, edges: Sequence[Tuple[int, int]]) -> Dict[Tuple[int, int], int]:
    adj = {i: set() for i in range(n)}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    out = {}
    for s in range(n):
        dist, frontier = {s: 0}, [s]
        while frontier:
            nxt = []
            for x in frontier:
                for y in adj[x]:
                    if y not in dist:
                        dist[y] = dist[x] + 1
                        nxt.append(y)
            frontier = nxt
        for t, d in dist.items():
            out[(s, t)] = d
    return out


def build_fragment(groups: Sequence[GroupCandidate], edges: Sequence[Tuple[int, int]], name: str,
                   prior: CouplingPrior, double: Optional[int] = None) -> Fragment:
    """Fragment from groups bonded along `edges`; `double` duplicates that group as a symmetric copy."""
    nodes = list(range(len(groups)))
    edges = list(edges)
    copy_of = {}
    if double is not None:
        neighbour = [b if a == double else a for a, b in edges if double in (a, b)]
        copy = len(nodes)
        nodes.append(copy)
        copy_of[copy] = double
        edges.append((neighbour[0], copy) if neighbour else (double, copy))
    group = lambda i: groups[copy_of.get(i, i)]
    sites, protons = [], []
    for i in nodes:
        g = group(i)
        sites.append(Site(f"{ELEMENT[g.nucleus]}{i + 1}", ELEMENT[g.nucleus]))
        protons.append(ProtonGroup(f"H{i + 1}", g.n_h, f"{ELEMENT[g.nucleus]}{i + 1}"))
    dist = _distances(len(nodes), edges)
    couplings = {}
    for i in nodes:
        for k in nodes:
            d = dist.get((i, k))
            if i == k:
                sign = -1.0 if group(i).nucleus == "15N" else 1.0
                couplings[pair(sites[i].label, protons[i].label)] = sign * group(i).j_hz
            elif d is not None:
                value = prior.heteronuclear(group(i).nucleus, d + 1)
                if value is not None:
                    couplings[pair(sites[i].label, protons[k].label)] = value
                if i < k:
                    value = prior.homonuclear(d + 2)
                    if value is not None:
                        couplings[pair(protons[i].label, protons[k].label)] = value
    symmetry = ()
    if double is not None:
        c = nodes[-1]
        symmetry = ({sites[double].label: sites[c].label, sites[c].label: sites[double].label,
                     protons[double].label: protons[c].label, protons[c].label: protons[double].label},)
    bonds = tuple((sites[a].label, sites[b].label) for a, b in edges)
    return Fragment(name, tuple(sites), tuple(protons), couplings, symmetry, bonds,
                    f"enumerated; couplings beyond 1J from prior '{prior.name}'")


def enumerate_fragments(candidates: Sequence[GroupCandidate], inventory: Inventory, top_per_band: int = 3,
                        max_groups: int = 4, prior: str = "generic_sp3", allow_equivalent: bool = True,
                        max_proposals: int = 30, open_band_penalty: float = 0.1) -> List[FragmentProposal]:
    """Proposals ordered by group count, then score minus `open_band_penalty` per weak band left open."""
    prior_obj = PRIORS[prior]
    out: List[FragmentProposal] = []
    seen = set()
    for combo, unexplained in covers(candidates, inventory, top_per_band, max_groups):
        score = sum(c.rank_score for c in combo) / len(combo)
        sources = sorted({s for c in combo for s in c.sources})
        label = " + ".join(f"{c.pattern}({c.j_hz:.0f})" for c in combo)
        options = []
        if len(combo) > 1:
            options.append(("isolated", None, None))
        for edges in _trees(len(combo)):
            options.append(("tree", edges, None))
            if allow_equivalent:
                for i, c in enumerate(combo):
                    leaf = sum(i in e for e in edges) <= 1
                    if c.n_h == 3 and leaf and len(combo) > 1:
                        options.append(("tree+equivalent", edges, i))
        for kind, edges, double in options:
            if kind == "isolated":
                frags = [build_fragment([c], [], f"{c.pattern}({c.j_hz:.0f})", prior_obj) for c in combo]
                topo = "isolated: " + " | ".join(f.name for f in frags)
            else:
                bond_text = ", ".join(f"{combo[a].pattern}-{combo[b].pattern}" for a, b in edges)
                extra = f"; 2 x {combo[double].pattern}" if double is not None else ""
                topo = f"bonded: {bond_text or combo[0].pattern}{extra}"
                frags = [build_fragment(combo, edges, f"{label} ({topo})", prior_obj, double)]
            key = (tuple(sorted((c.pattern, round(c.j_hz, 1)) for c in combo)), topo)
            if key in seen:
                continue
            seen.add(key)
            out.append(FragmentProposal(frags, list(combo), topo, score, sources, unexplained))
    out.sort(key=lambda p: (p.n_groups, -(p.score - open_band_penalty * len(p.unexplained)),
                            p.topology.startswith("isolated"), "equivalent" in p.topology))
    return out[:max_proposals]
