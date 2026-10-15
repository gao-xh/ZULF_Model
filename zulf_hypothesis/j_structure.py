"""Route A of J -> structure: rank heavy-atom graphs by how well they explain an observed J network.

The observation is what a fit (or a spectrum reading) gives: carbon units (one per 13C site, with the number of
protons on it from the 1J pattern and the number of symmetry-equivalent copies), proton groups attached to a copy of
a unit, and the fitted couplings between units and proton groups. A candidate structure is a connected graph over
the carbon copies and up to `max_unseen` unseen heavy atoms X (N, O, or a carbon without protons: no 1J, so no
13C line of its own), with single, double or triple bonds. Each observed coupling has a bond count in the candidate (C to a proton: carbon-carbon
distance + 1; proton to proton: distance + 2), and a coupling-size likelihood p(J | nuclei, bonds) scores it
(`couplings_likelihood`: approximate literature ranges, broadened by the fit's sigma and with a wide tail for
unusual values). Bond orders do not change bond counts; they enter through the 1J of each carbon, whose size depends
on its hybridization (`one_bond_likelihood`: sp3 near 125-145 Hz, sp2 near 150-175, sp near 250). log L = sum over
the observed couplings + prior (each unseen atom, ring, multiple bond and missing valence costs a fixed amount). Candidates must keep the copies of a unit equivalent (equal Weisfeiler-Lehman colours) and
distinct units distinct.

The ranking is a conditional result: it orders the graphs this enumeration allows under these likelihoods; it does
not determine a structure. Route B (zulf_model.structure.edge_model) learns the edge evidence from generator
pairs and plugs in through `likelihood`; the existing
hypothesis pipeline (zulf_hypothesis.search) fits candidate fragments to the spectrum directly.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# coupling-size models: (lo, hi) uniform ranges with weights, or "long" (a spike at 0 plus a Laplace tail)
CH_RANGES = {1: [((115.0, 175.0), 1.0)], 2: [((-7.5, 4.0), 0.9), ((-12.0, 12.0), 0.1)],
             3: [((0.0, 9.0), 0.8), ((-4.0, 12.0), 0.2)]}
HH_RANGES = {2: [((-17.0, -6.0), 0.8), ((-4.0, 4.0), 0.2)], 3: [((4.0, 9.5), 0.8), ((0.0, 18.0), 0.2)],
             4: "long+meta"}
TAIL = 0.02                    # probability mass of a value outside every rule (wide uniform -20..20 Hz)
LONG_SPIKE, LONG_SCALE = 0.5, 0.4


def _uniform_conv(j, lo, hi, s):
    """Density of U(lo, hi) convolved with N(0, s) at j."""
    a = (j - lo) / (s * math.sqrt(2.0))
    b = (j - hi) / (s * math.sqrt(2.0))
    return 0.5 * (math.erf(a) - math.erf(b)) / (hi - lo)


def _long(j, s):
    spike = math.exp(-0.5 * (j / s) ** 2) / (s * math.sqrt(2 * math.pi))
    lap = math.exp(-abs(j) / LONG_SCALE) / (2 * LONG_SCALE)
    return LONG_SPIKE * spike + (1 - LONG_SPIKE) * lap


MAX_BONDS = 6                  # bond counts above this share the long-range model (route A) / class (route B)


ONE_BOND = {"sp3": [((118.0, 150.0), 0.95), ((105.0, 185.0), 0.05)],
            "sp2": [((145.0, 180.0), 0.95), ((105.0, 235.0), 0.05)],
            "sp": [((240.0, 260.0), 0.95), ((200.0, 270.0), 0.05)]}


def one_bond_likelihood(j: float, hyb: str, sigma: float = 0.1, h=None) -> float:
    """p(1J(C,H) | hybridization of the carbon), up to a factor shared by all hybridizations."""
    s = math.hypot(max(sigma, 0.05), 0.3)
    return (1 - TAIL) * sum(w * _uniform_conv(j, lo, hi, s) for (lo, hi), w in ONE_BOND[hyb]) + TAIL / 300.0


def hybridization(orders: Sequence[int]) -> str:
    """Hybridization of a carbon from the orders of its bonds."""
    if 3 in orders or list(orders).count(2) >= 2:
        return "sp"
    return "sp2" if 2 in orders else "sp3"


def couplings_likelihood(j: float, kind: str, bonds: int, sigma: float = 0.1, h=None) -> float:
    """p(J | kind, bonds), kind "CH" (13C to proton) or "HH". bonds >= 4 (CH) / 5 (HH): long-range model. h (the
    proton counts of the carbons at the two ends) is accepted for the likelihood interface and not used."""
    s = math.hypot(max(sigma, 0.05), 0.3)
    table = CH_RANGES if kind == "CH" else HH_RANGES
    rule = table.get(bonds)
    if rule is None or rule == "long+meta":
        p = _long(j, s)
        if rule == "long+meta":
            p = 0.9 * p + 0.1 * _uniform_conv(j, 1.0, 3.0, s)
    else:
        p = sum(w * _uniform_conv(j, lo, hi, s) for (lo, hi), w in rule)
    return (1 - TAIL) * p + TAIL / 40.0


@dataclass
class JObservation:
    """units {label: (protons, copies)}; protons {group: (unit label, copy)}; couplings {(a, b): J} with a, b unit
    labels (copy 0) or proton groups; sigma {(a, b): Hz} optional; isotopes {unit or group: isotope} optional
    (units default to 13C, proton groups to 1H; e.g. {"HN": "2H"} for an exchanged group)."""
    units: Dict[str, Tuple[int, int]]
    protons: Dict[str, Tuple[str, int]]
    couplings: Dict[Tuple[str, str], float]
    sigma: Dict[Tuple[str, str], float] = field(default_factory=dict)
    isotopes: Dict[str, str] = field(default_factory=dict)

    def isotope(self, label: str) -> str:
        return self.isotopes.get(label, "1H" if label in self.protons else "13C")

    def labelled(self) -> bool:
        """True when some unit or group is not the default 13C / 1H (route A's ranges assume the defaults)."""
        return any(self.isotope(k) != ("1H" if k in self.protons else "13C") for k in self.isotopes)

    @classmethod
    def from_dict(cls, data: dict) -> "JObservation":
        def key(k):
            a, b = k[2:-1].split(",") if k.startswith("J(") else k.split(",")
            return a.strip(), b.strip()
        units = {u: (int(v["h"]), int(v.get("copies", 1))) for u, v in data["units"].items()}
        protons = {g: (v[0], int(v[1])) for g, v in data["protons"].items()}
        return cls(units, protons, {key(k): float(v) for k, v in data["couplings"].items()},
                   {key(k): float(v) for k, v in data.get("sigma", {}).items()}, dict(data.get("isotopes", {})))

    def atoms(self) -> List[Tuple[str, int]]:
        return [(u, c) for u, (_, copies) in self.units.items() for c in range(copies)]


@dataclass
class Candidate:
    edges: Tuple[Tuple[int, int], ...]
    unseen: int
    log_likelihood: float
    log_prior: float
    explanation: List[dict]
    orders: Tuple[int, ...] = ()

    @property
    def score(self) -> float:
        return self.log_likelihood + self.log_prior


def _distances(n: int, edges) -> np.ndarray:
    d = np.full((n, n), 10 ** 6, int)
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    for s in range(n):
        d[s, s] = 0
        frontier = [s]
        while frontier:
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if d[s, v] > d[s, u] + 1:
                        d[s, v] = d[s, u] + 1
                        nxt.append(v)
            frontier = nxt
    return d


def _wl(n, edges, labels, rounds=4, orders=None):
    orders = orders or [1] * len(edges)
    adj = [[] for _ in range(n)]
    for (a, b), o in zip(edges, orders):
        adj[a].append((b, o))
        adj[b].append((a, o))
    col = list(labels)
    for _ in range(rounds):
        col = [hash((col[i], tuple(sorted((col[k], o) for k, o in adj[i])))) for i in range(n)]
    return col


def describe(obs: JObservation, cand: Candidate) -> str:
    """Bond list with unit labels, e.g. C1-C2, C1-C2', C1-X1."""
    atoms = obs.atoms()
    names = [u + "'" * c for u, c in atoms] + [f"X{k + 1}" for k in range(cand.unseen)]
    orders = cand.orders or (1,) * len(cand.edges)
    return ", ".join(f"{names[a]}{BOND_SYMBOL[o]}{names[b]}" for (a, b), o in zip(cand.edges, orders))


BOND_SYMBOL = {1: "-", 2: "=", 3: "#"}


def rank_structures(obs: JObservation, max_unseen: int = 2, unseen_valence: int = 3, ring_cost: float = 2.0,
                    unseen_cost: float = 1.0, missing_valence_cost: float = 3.0, max_rings: int = 1,
                    top: int = 10, likelihood=None, double_cost: float = 1.5, triple_cost: float = 3.0,
                    max_multiple: int = 2, unseen_max_valence: int = 4) -> List[Candidate]:
    """All connected graphs over the carbon copies and 0..max_unseen unseen atoms (degree caps 4 - protons for
    carbons and unseen_valence for X, at most max_rings independent cycles; then every assignment of up to
    max_multiple double or triple bonds within the valences, 4 - protons for carbons and unseen_max_valence for X), with the copies of a unit equivalent; ranked by
    log L + log prior.

    likelihood(j, kind, bonds, sigma, h) -> p(J | bonds) up to a factor that does not depend on bonds; default
    `couplings_likelihood` (route A); route B passes zulf_model.structure.edge_model.LearnedLikelihood. h is
    (protons on the carbon, protons on the proton's carbon) for CH and the sorted pair for HH. The 1J terms use
    likelihood.hybrid(j, hyb, sigma, h) when the likelihood has it, else `one_bond_likelihood`. A likelihood with
    `accepts_nuclei` also gets nuclei= (the isotopes of the pair, in the order of h; route B in mode K uses them)."""
    likelihood = likelihood or couplings_likelihood
    hybrid = getattr(likelihood, "hybrid", one_bond_likelihood)
    with_nuclei = getattr(likelihood, "accepts_nuclei", False)
    hybrid_nuclei = with_nuclei and hasattr(likelihood, "hybrid")
    atoms = obs.atoms()
    nc = len(atoms)
    cap_c = [4 - obs.units[u][0] for u, _ in atoms]
    index = {a: i for i, a in enumerate(atoms)}
    # observed couplings as (atom index, kind, offset): bonds = distance + offset
    terms, one_bond = [], []
    for (a, b), j in obs.couplings.items():
        ends = []
        for g in (a, b):
            if g in obs.protons:
                ends.append(("H", index[obs.protons[g]]))
            else:
                ends.append(("C", index[(g, 0)]))
        kinds = "".join(sorted(e[0] for e in ends))
        if kinds == "CC":
            continue                                     # 13C-13C couplings: not used (double labels)
        offset = 1 if kinds == "CH" else 2
        same_proton_group = kinds == "HH" and a == b
        if same_proton_group:
            continue
        hs = [obs.units[g][0] if g not in obs.protons else obs.units[obs.protons[g][0]][0] for g in (a, b)]
        h = tuple(hs) if ends[0][0] == "C" else (tuple(hs[::-1]) if kinds == "CH" else tuple(sorted(hs)))
        sig = obs.sigma.get((a, b), 0.1)
        extra = {}
        if with_nuclei:
            pair = (a, b) if ends[0][0] == "C" or kinds == "HH" else (b, a)
            extra["nuclei"] = (obs.isotope(pair[0]), obs.isotope(pair[1]))
        table = [0.0] + [math.log(likelihood(j, kinds, b_, sig, h, **extra)) for b_ in range(1, MAX_BONDS + 1)]
        if kinds == "CH" and ends[0][1] == ends[1][1]:
            hyb_extra = extra if hybrid_nuclei else {}
            hyb_table = {hy: math.log(hybrid(j, hy, sig, h, **hyb_extra)) for hy in ONE_BOND}
            one_bond.append((a, b, ends[0][1], j, table[1], hyb_table))
            continue
        terms.append((a, b, ends[0][1], ends[1][1], kinds, offset, j, table))
    best: Dict[tuple, Candidate] = {}
    for nx in range(max_unseen + 1):
        n = nc + nx
        pairs = [(a, b) for a in range(n) for b in range(a + 1, n) if a < nc or b < nc]   # no X-X bonds
        caps = cap_c + [unseen_valence] * nx
        # backtracking over edges with degree caps
        deg = [0] * n
        chosen: List[Tuple[int, int]] = []

        def emit():
            m = len(chosen)
            if m < n - 1 or m - (n - 1) > max_rings:
                return
            d = _distances(n, chosen)
            if (d >= 10 ** 6).any():
                return
            if any(deg[k] == 0 for k in range(nc, n)):
                return
            ll0, expl0 = 0.0, []
            for a, b, ia, ib, kinds, offset, j, table in terms:
                bonds = int(d[ia, ib]) + offset
                lp = table[min(bonds, MAX_BONDS)]
                ll0 += lp
                expl0.append({"coupling": f"J({a},{b})", "J": j, "bonds": bonds, "log_p": lp})
            rings = m - (n - 1)
            edges = tuple(chosen)
            free = [caps[k] - deg[k] for k in range(nc)] + [unseen_max_valence - deg[k] for k in range(nc, n)]
            for orders in _bond_orders(edges, free, nc, max_multiple):
                score_orders(n, nx, edges, orders, ll0, expl0, rings)

        def score_orders(n, nx, edges, orders, ll0, expl0, rings):
            labels = [hash((u, obs.units[u][0])) for u, _ in atoms] + [hash("X")] * nx
            col = _wl(n, edges, labels, orders=orders)
            for u in obs.units:                          # copies equivalent, distinct units distinct
                cs = {col[index[(u, c)]] for c in range(obs.units[u][1])}
                if len(cs) != 1:
                    return
            ucol = [col[index[(u, 0)]] for u in obs.units]
            if len(set(ucol)) != len(ucol):
                return
            bond_orders = [[] for _ in range(n)]
            valence = [0] * n
            for (a, b), o in zip(edges, orders):
                bond_orders[a].append(o)
                bond_orders[b].append(o)
                valence[a] += o
                valence[b] += o
            ll, expl = ll0, list(expl0)
            for a, b, ia, j, lp_bond, hyb_table in one_bond:
                hy = hybridization(bond_orders[ia])
                lp = lp_bond + hyb_table[hy]
                ll += lp
                expl.append({"coupling": f"J({a},{b})", "J": j, "bonds": 1, "hybridization": hy, "log_p": lp})
            missing = sum(cap_c[k] - valence[k] for k in range(nc))
            prior = (-unseen_cost * nx - ring_cost * rings - missing_valence_cost * missing
                     - double_cost * orders.count(2) - triple_cost * orders.count(3))
            canon = canonical(n, list(edges), labels, orders)
            cand = Candidate(edges, nx, ll, prior, expl, orders)
            if canon not in best or cand.score > best[canon].score:
                best[canon] = cand              # isomorphic relabelings: keep the one that explains the data best

        def walk(k):
            if k == len(pairs):
                emit()
                return
            a, b = pairs[k]
            walk(k + 1)
            if deg[a] < caps[a] and deg[b] < caps[b] and len(chosen) < n - 1 + max_rings:
                deg[a] += 1
                deg[b] += 1
                chosen.append((a, b))
                walk(k + 1)
                chosen.pop()
                deg[a] -= 1
                deg[b] -= 1
        walk(0)
    results = sorted(best.values(), key=lambda c: -c.score)
    return results[:top]


def _bond_orders(edges, free, nc, max_multiple):
    """Every order assignment (1, 2 or 3 per edge) with at most max_multiple multiple bonds and the extra valence
    of each atom within free[k]."""
    orders = [1] * len(edges)
    left = list(free)

    def rec(k, used):
        if k == len(edges):
            yield tuple(orders)
            return
        yield from rec(k + 1, used)
        if used >= max_multiple:
            return
        a, b = edges[k]
        for extra in (1, 2):
            if left[a] >= extra and left[b] >= extra:
                orders[k] = 1 + extra
                left[a] -= extra
                left[b] -= extra
                yield from rec(k + 1, used + 1)
                left[a] += extra
                left[b] += extra
                orders[k] = 1

    yield from rec(0, 0)


def canonical(n, edges, labels, orders=None) -> tuple:
    """Isomorphism key of a labelled graph (WL colours of the nodes and of the edges, with bond orders)."""
    orders = list(orders) if orders else [1] * len(edges)
    col = _wl(n, edges, labels, orders=orders)
    return (n, tuple(sorted(col)), tuple(sorted((*sorted((col[a], col[b])), o) for (a, b), o in zip(edges, orders))))


def parse_bonds(obs: JObservation, text: str, with_orders: bool = False):
    """'C1-C2, C1=C2', C1-X1' -> (number of unseen atoms, edges[, orders]) in the indexing of rank_structures
    (bond symbols - = #)."""
    atoms = obs.atoms()
    names = {u + "'" * c: i for i, (u, c) in enumerate(atoms)}
    bonds = []
    for bond in text.split(","):
        bond = bond.strip()
        sym = next(c for c in "-=#" if c in bond)
        a, b = bond.split(sym)
        bonds.append((a.strip(), b.strip(), {"-": 1, "=": 2, "#": 3}[sym]))
    xs = sorted({t for a, b, _ in bonds for t in (a, b) if t.startswith("X")})
    for k, x in enumerate(xs):
        names[x] = len(atoms) + k
    edges = [tuple(sorted((names[a], names[b]))) for a, b, _ in bonds]
    if with_orders:
        return len(xs), edges, [o for *_, o in bonds]
    return len(xs), edges


def same_structure(obs: JObservation, cand: Candidate, text: str) -> bool:
    """True when the candidate and the bond list are the same structure (unit labels and bond orders respected)."""
    nx, edges, orders = parse_bonds(obs, text, with_orders=True)
    atoms = obs.atoms()
    labels = lambda k: [hash((u, obs.units[u][0])) for u, _ in atoms] + [hash("X")] * k
    return (nx == cand.unseen and canonical(len(atoms) + nx, edges, labels(nx), orders)
            == canonical(len(atoms) + cand.unseen, list(cand.edges), labels(cand.unseen), cand.orders))
