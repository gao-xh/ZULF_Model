"""Spin systems, components and interpretations.

A `SpinSystem` is a list of nucleus symbols and a real symmetric J matrix in Hz
with zero diagonal. Magnetic-equivalence groups are either supplied and
verified, or detected. Nothing in this module assumes a particular spin count.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

from .nuclei import get_registry

EQUIVALENCE_ATOL_HZ = 1e-9


def _as_groups(groups) -> Tuple[Tuple[int, ...], ...]:
    return tuple(tuple(int(i) for i in g) for g in groups)


def detect_equivalence_groups(isotopes: Sequence[str], couplings: np.ndarray,
                              atol: float = EQUIVALENCE_ATOL_HZ) -> Tuple[Tuple[int, ...], ...]:
    """Coarsest magnetic-equivalence partition (singletons included).

    Spins in a group share a nucleus and have identical couplings to every spin
    outside the group. Couplings inside a group are unobservable and ignored.
    Starting from the partition by nucleus, classes are split by their coupling
    profile to spins outside the class until stable; no valid partition is
    ever split, so the result is the unique coarsest one.
    """
    n = len(isotopes)
    classes: List[List[int]] = []
    for symbol in dict.fromkeys(isotopes):
        classes.append([i for i in range(n) if isotopes[i] == symbol])
    changed = True
    while changed:
        changed = False
        refined: List[List[int]] = []
        for cls in classes:
            outside = [m for m in range(n) if m not in cls]
            buckets: List[List[int]] = []
            for i in cls:
                for bucket in buckets:
                    if np.allclose(couplings[i, outside], couplings[bucket[0], outside], atol=atol, rtol=0):
                        bucket.append(i)
                        break
                else:
                    buckets.append([i])
            if len(buckets) > 1:
                changed = True
            refined.extend(buckets)
        classes = refined
    return tuple(tuple(sorted(g)) for g in sorted(classes, key=min))


def verify_equivalence_groups(isotopes, couplings, groups, atol=EQUIVALENCE_ATOL_HZ) -> None:
    n = len(isotopes)
    seen = set()
    for group in groups:
        if not group or any(type(i) is not int or not 0 <= i < n or i in seen for i in group):
            raise ValueError("Equivalence groups must be disjoint zero-based spin indices.")
        seen.update(group)
        if len({isotopes[i] for i in group}) != 1:
            raise ValueError("An equivalence group mixes nuclei.")
        outside = [m for m in range(n) if m not in group]
        rows = couplings[np.ix_(list(group), outside)]
        if len(group) > 1 and outside and not np.allclose(rows, rows[0], atol=atol, rtol=0):
            raise ValueError(f"Group {group} is not magnetically equivalent.")
    if seen != set(range(n)):
        raise ValueError("Equivalence groups must cover every spin (use singletons).")


@dataclass(frozen=True)
class SpinSystem:
    isotopes: Tuple[str, ...]
    couplings_hz: np.ndarray
    groups: Optional[Tuple[Tuple[int, ...], ...]] = None

    def __post_init__(self):
        registry = get_registry()
        isotopes = tuple(str(s) for s in self.isotopes)
        for symbol in isotopes:
            registry[symbol]
        j = np.array(self.couplings_hz, dtype=float, copy=True)
        n = len(isotopes)
        if n < 1 or j.shape != (n, n) or not np.isfinite(j).all():
            raise ValueError("Couplings must be a finite n x n matrix for n >= 1 spins.")
        if not np.allclose(j, j.T, atol=1e-12, rtol=0) or np.any(np.diag(j) != 0):
            raise ValueError("Couplings must be symmetric with zero diagonal; no implicit repair.")
        j = (j + j.T) / 2
        j.setflags(write=False)
        object.__setattr__(self, "isotopes", isotopes)
        object.__setattr__(self, "couplings_hz", j)
        if self.groups is None:
            groups = detect_equivalence_groups(isotopes, j)
        else:
            groups = _as_groups(self.groups)
            verify_equivalence_groups(isotopes, j, groups)
        object.__setattr__(self, "groups", tuple(sorted(groups, key=min)))

    @property
    def n_spins(self) -> int:
        return len(self.isotopes)

    def group_of(self) -> np.ndarray:
        index = np.empty(self.n_spins, dtype=int)
        for g, members in enumerate(self.groups):
            index[list(members)] = g
        return index

    def observable_couplings(self) -> np.ndarray:
        """Couplings with entries inside equivalence groups set to zero (unobservable)."""
        g = self.group_of()
        out = self.couplings_hz.copy()
        out[g[:, None] == g[None, :]] = 0.0
        return out

    def group_couplings(self) -> np.ndarray:
        """Coupling between equivalence groups; diagonal is zero (unobservable)."""
        reps = [g[0] for g in self.groups]
        out = self.couplings_hz[np.ix_(reps, reps)].copy()
        np.fill_diagonal(out, 0.0)
        return out

    def permute(self, order: Sequence[int]) -> "SpinSystem":
        """Return the system with spin `k` of the result equal to spin `order[k]`."""
        order = [int(i) for i in order]
        if sorted(order) != list(range(self.n_spins)):
            raise ValueError("order must be a permutation.")
        inverse = np.argsort(order)
        groups = tuple(tuple(sorted(int(inverse[i]) for i in g)) for g in self.groups)
        return SpinSystem(tuple(self.isotopes[i] for i in order),
                          self.couplings_hz[np.ix_(order, order)], groups)

    def composition(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for s in self.isotopes:
            out[s] = out.get(s, 0) + 1
        return out

    def group_signature(self) -> Tuple[Tuple[str, int], ...]:
        """Sorted multiset of (nucleus, group size)."""
        return tuple(sorted((self.isotopes[g[0]], len(g)) for g in self.groups))

    def to_dict(self) -> dict:
        return {"isotopes": list(self.isotopes), "couplings_hz": self.couplings_hz.tolist(),
                "groups": [list(g) for g in self.groups]}

    @classmethod
    def from_dict(cls, data: dict) -> "SpinSystem":
        groups = data.get("groups")
        return cls(tuple(data["isotopes"]), np.asarray(data["couplings_hz"], float),
                   None if groups is None else _as_groups(groups))

    @classmethod
    def from_group_couplings(cls, group_isotopes: Sequence[str], group_sizes: Sequence[int],
                             group_couplings: np.ndarray) -> "SpinSystem":
        """Expand a group-level description into a spin-level system."""
        g = np.asarray(group_couplings, float)
        if g.shape != (len(group_sizes),) * 2:
            raise ValueError("Group coupling matrix shape mismatch.")
        index = np.repeat(np.arange(len(group_sizes)), group_sizes)
        j = g[np.ix_(index, index)].copy()
        j[index[:, None] == index[None, :]] = 0.0
        isotopes = tuple(np.repeat(np.asarray(group_isotopes, dtype=object), group_sizes))
        groups, start = [], 0
        for size in group_sizes:
            groups.append(tuple(range(start, start + size)))
            start += size
        return cls(isotopes, j, tuple(groups))


# ---------------------------------------------------------------------------
# Canonical ordering (training convenience only; never used for equality)
# ---------------------------------------------------------------------------

def _refined_colors(system: SpinSystem, nuclei_order: Sequence[str], rounds: int = 3,
                    decimals: int = 3) -> List[tuple]:
    j = np.round(system.observable_couplings(), decimals)
    order = {s: i for i, s in enumerate(nuclei_order)}
    sizes = {i: len(g) for g in system.groups for i in g}
    colors = [(order.get(s, len(order)), -sizes[i]) for i, s in enumerate(system.isotopes)]
    for _ in range(rounds):
        colors = [(colors[i], tuple(sorted(((colors[k], -abs(j[i, k]), j[i, k])
                                            for k in range(system.n_spins) if k != i), reverse=True)))
                  for i in range(system.n_spins)]
        # Compress to ranks so nested tuples stay small.
        ranks = {c: r for r, c in enumerate(sorted(set(colors)))}
        colors = [(ranks[c],) for c in colors]
    return colors


def canonical_group_order(system: SpinSystem, nuclei_order: Sequence[str]) -> List[Tuple[int, ...]]:
    """Order equivalence groups by nucleus, size and refined coupling colors."""
    j = system.observable_couplings()
    order = {s: i for i, s in enumerate(nuclei_order)}
    colors = _refined_colors(system, nuclei_order)

    def key(group):
        i = group[0]
        profile = tuple(sorted((-abs(round(float(v), 3)) for k, v in enumerate(j[i]) if k not in group)))
        return (order.get(system.isotopes[i], len(order)), -len(group), profile, colors[i], group)

    return sorted(system.groups, key=key)


def canonical_permutation(system: SpinSystem, nuclei_order: Sequence[str]) -> List[int]:
    return [i for group in canonical_group_order(system, nuclei_order) for i in group]


def canonical_sign(system: SpinSystem) -> int:
    """Global sign making the largest-magnitude observable heteronuclear coupling positive.

    Falls back to the largest observable coupling of any kind; returns +1 for a
    system without observable couplings.
    """
    j = system.observable_couplings()
    iso = np.array(system.isotopes, dtype=object)
    hetero = iso[:, None] != iso[None, :]
    for mask in (hetero, np.ones_like(hetero)):
        values = np.where(mask, j, 0.0)
        idx = np.unravel_index(np.argmax(np.abs(values)), values.shape)
        if abs(values[idx]) > 0:
            return 1 if values[idx] > 0 else -1
    return 1


def canonicalize(system: SpinSystem, nuclei_order: Sequence[str], fix_sign: bool = True) -> SpinSystem:
    """Canonical spin order and (optionally) canonical global J sign."""
    if fix_sign and canonical_sign(system) < 0:
        system = SpinSystem(system.isotopes, -system.couplings_hz, system.groups)
    return system.permute(canonical_permutation(system, nuclei_order))


# ---------------------------------------------------------------------------
# Permutation matching
# ---------------------------------------------------------------------------

@dataclass
class MatchResult:
    permutation: Optional[List[int]]
    rms_error_hz: float
    max_abs_error_hz: float
    exact: bool
    nodes: int
    sign: int = 1

    @property
    def matched(self) -> bool:
        return self.permutation is not None


def best_permutation(reference: SpinSystem, other: SpinSystem, node_limit: int = 200000,
                     allow_global_sign: bool = True) -> MatchResult:
    """Find the isotope-preserving relabelling of `other` closest to `reference`.

    Returns `permutation` and `sign` such that `sign * other.permute(permutation)`
    is compared with `reference` on observable couplings. A global sign flip of
    every J leaves a zero-field spectrum unchanged under the default protocol,
    so it is allowed by default. Uses branch and bound on the squared coupling
    error; `exact` is False when the node limit stopped the search.
    """
    if allow_global_sign:
        plus = best_permutation(reference, other, node_limit, False)
        flipped = SpinSystem(other.isotopes, -other.couplings_hz, other.groups)
        minus = best_permutation(reference, flipped, node_limit, False)
        if minus.matched and minus.rms_error_hz < plus.rms_error_hz:
            minus.sign = -1
            return minus
        return plus
    n = reference.n_spins
    if other.n_spins != n or sorted(reference.isotopes) != sorted(other.isotopes):
        return MatchResult(None, math.inf, math.inf, True, 0)
    a, b = reference.observable_couplings(), other.observable_couplings()
    # Heuristic seed: assignment on sorted absolute coupling profiles.
    pa = np.sort(np.abs(a), axis=1)
    pb = np.sort(np.abs(b), axis=1)
    cost = ((pa[:, None, :] - pb[None, :, :]) ** 2).sum(-1)
    iso_a = np.array(reference.isotopes, dtype=object)
    iso_b = np.array(other.isotopes, dtype=object)
    cost = np.where(iso_a[:, None] == iso_b[None, :], cost, 1e30)
    rows, cols = linear_sum_assignment(cost)
    seed = [0] * n
    for r, c in zip(rows, cols):
        seed[r] = int(c)

    def total(perm):
        return float(((a - b[np.ix_(perm, perm)]) ** 2).sum())

    best_perm, best_cost = list(seed), total(seed)
    order = list(np.argsort(-np.abs(a).sum(1), kind="stable"))
    nodes = 0
    exhausted = True
    assigned = [-1] * n
    used = [False] * n

    def search(depth, partial):
        nonlocal best_perm, best_cost, nodes, exhausted
        if nodes >= node_limit:
            exhausted = False
            return
        nodes += 1
        if depth == n:
            if partial < best_cost:
                best_cost, best_perm = partial, list(assigned)
            return
        i = order[depth]
        placed = [order[d] for d in range(depth)]
        for c in range(n):
            if used[c] or iso_b[c] != iso_a[i]:
                continue
            add = 0.0
            for p in placed:
                add += 2 * (a[i, p] - b[c, assigned[p]]) ** 2
            if partial + add >= best_cost - 1e-15:
                continue
            used[c] = True
            assigned[i] = c
            search(depth + 1, partial + add)
            used[c] = False
            assigned[i] = -1

    search(0, 0.0)
    diff = a - b[np.ix_(best_perm, best_perm)]
    pairs = n * (n - 1) / 2
    rms = math.sqrt(float((diff[np.triu_indices(n, 1)] ** 2).sum()) / pairs) if pairs else 0.0
    return MatchResult(best_perm, rms, float(np.abs(diff).max()) if n > 1 else 0.0, exhausted, nodes)


# ---------------------------------------------------------------------------
# Components and interpretations
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Component:
    system: SpinSystem
    contribution: float = 1.0
    label: str = ""
    metadata: dict = field(default_factory=dict, compare=False)

    def __post_init__(self):
        if not math.isfinite(self.contribution) or self.contribution < 0:
            raise ValueError("Component contribution must be finite and nonnegative.")

    def to_dict(self) -> dict:
        return {"system": self.system.to_dict(), "contribution": self.contribution,
                "label": self.label, "metadata": self.metadata}

    @classmethod
    def from_dict(cls, data: dict) -> "Component":
        return cls(SpinSystem.from_dict(data["system"]), float(data.get("contribution", 1.0)),
                   data.get("label", ""), dict(data.get("metadata", {})))


@dataclass(frozen=True)
class Interpretation:
    components: Tuple[Component, ...]
    score: Optional[float] = None
    metadata: dict = field(default_factory=dict, compare=False)

    def __post_init__(self):
        object.__setattr__(self, "components", tuple(self.components))
        if not self.components:
            raise ValueError("An interpretation needs at least one component.")

    def __len__(self) -> int:
        return len(self.components)

    def to_dict(self) -> dict:
        return {"components": [c.to_dict() for c in self.components], "score": self.score,
                "metadata": self.metadata}

    @classmethod
    def from_dict(cls, data: dict) -> "Interpretation":
        return cls(tuple(Component.from_dict(c) for c in data["components"]), data.get("score"),
                   dict(data.get("metadata", {})))

    def key(self, nuclei_order: Sequence[str], decimals: int = 1) -> tuple:
        """Hashable approximate identity for deduplicating candidate lists."""
        parts = []
        for c in self.components:
            s = canonicalize(c.system, nuclei_order)
            obs = s.observable_couplings()
            parts.append((s.isotopes, tuple(np.round(obs[np.triu_indices(s.n_spins, 1)], decimals) + 0.0)))
        return tuple(sorted(parts))


def all_isotope_preserving_permutations(isotopes: Sequence[str]):
    """Enumerate permutations that preserve the isotope at every position (small n only)."""
    classes: Dict[str, List[int]] = {}
    for i, s in enumerate(isotopes):
        classes.setdefault(s, []).append(i)
    keys = list(classes)
    for choice in itertools.product(*(itertools.permutations(classes[k]) for k in keys)):
        perm = [0] * len(isotopes)
        for k, images in zip(keys, choice):
            for src, dst in zip(classes[k], images):
                perm[src] = dst
        yield perm
