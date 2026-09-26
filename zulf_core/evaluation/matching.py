"""Operational comparison of interpretations (open question Q4, current definition).

* composition match: same nuclei multiset and same magnetic-equivalence group
  signature (nucleus, group size);
* structure match: a one-to-one pairing of components that all match in
  composition (same component count);
* J match: after structure match, every paired component agrees on observable
  couplings within `tolerance_hz` (plus `relative` x |J|), up to spin
  relabelling and one global sign per connected block of the true coupling
  network (the relative sign of disconnected blocks is unobservable; D25).

These are recorded conventions, not claims of physical identity.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

from ..spinsystem import Interpretation, MatchResult, SpinSystem, best_permutation


def composition_match(a: SpinSystem, b: SpinSystem) -> bool:
    return sorted(a.isotopes) == sorted(b.isotopes) and a.group_signature() == b.group_signature()


def coupling_blocks(system: SpinSystem, threshold_hz: float = 1e-9) -> List[List[int]]:
    """Connected components of spins linked by observable couplings above `threshold_hz`."""
    j = np.abs(system.observable_couplings()) > threshold_hz
    n = system.n_spins
    seen, blocks = set(), []
    for start in range(n):
        if start in seen:
            continue
        stack, block = [start], []
        seen.add(start)
        while stack:
            i = stack.pop()
            block.append(i)
            for k in np.flatnonzero(j[i]):
                if int(k) not in seen:
                    seen.add(int(k))
                    stack.append(int(k))
        blocks.append(sorted(block))
    return blocks


def aligned_couplings(truth: SpinSystem, other: SpinSystem) -> Tuple[Optional[np.ndarray], MatchResult]:
    """Observable couplings of `other` relabelled (and sign-fixed per truth block) onto `truth`.

    One global sign is free for a connected truth network. When the truth
    network splits into disconnected blocks, each block's sign is chosen
    independently, and the relabelling is also searched on |J| so that
    opposite block signs cannot mislead it.
    """
    t = truth.observable_couplings()
    blocks = coupling_blocks(truth)
    options = [best_permutation(truth, other)]
    if len(blocks) > 1:
        options.append(best_permutation(SpinSystem(truth.isotopes, np.abs(truth.couplings_hz), truth.groups),
                                        SpinSystem(other.isotopes, np.abs(other.couplings_hz), other.groups),
                                        allow_global_sign=False))
    best, best_match, best_err = None, options[0], math.inf
    for match in options:
        if not match.matched:
            continue
        aligned = match.sign * other.permute(match.permutation).observable_couplings()
        if len(blocks) > 1:
            for block in blocks:
                idx = np.ix_(block, block)
                if np.abs(-aligned[idx] - t[idx]).max(initial=0.0) < np.abs(aligned[idx] - t[idx]).max(initial=0.0):
                    aligned[idx] = -aligned[idx]
        err = float(np.abs(aligned - t).max()) if truth.n_spins > 1 else 0.0
        if err < best_err:
            best, best_match, best_err = aligned, match, err
    return best, best_match


def coupling_error(truth: SpinSystem, other: SpinSystem, relative: float = 0.0) -> Tuple[float, MatchResult]:
    """Max over couplings of |dJ| - relative * |J_true| after best relabelling (Hz)."""
    aligned, match = aligned_couplings(truth, other)
    if aligned is None:
        return math.inf, match
    t = truth.observable_couplings()
    excess = np.abs(aligned - t) - relative * np.abs(t)
    return float(excess.max()) if truth.n_spins > 1 else 0.0, match


@dataclass
class InterpretationMatch:
    structure: bool
    j_match: bool
    pairs: List[Tuple[int, int]] = field(default_factory=list)
    max_j_error_hz: float = math.inf
    component_errors_hz: List[float] = field(default_factory=list)


def match_interpretations(truth: Interpretation, candidate: Interpretation, tolerance_hz: float = 1.0,
                          relative: float = 0.0) -> InterpretationMatch:
    t, c = truth.components, candidate.components
    if len(t) != len(c):
        return InterpretationMatch(False, False)
    n = len(t)
    cost = np.full((n, n), 1e12)
    errors = {}
    for i in range(n):
        for j in range(n):
            if composition_match(t[i].system, c[j].system):
                err, _ = coupling_error(t[i].system, c[j].system, relative)
                errors[(i, j)] = err
                cost[i, j] = err
    rows, cols = linear_sum_assignment(cost)
    if any(cost[r, k] >= 1e12 for r, k in zip(rows, cols)):
        return InterpretationMatch(False, False)
    comp_errors = [errors[(r, k)] for r, k in zip(rows, cols)]
    worst = max(comp_errors)
    return InterpretationMatch(True, worst <= tolerance_hz, list(zip(rows.tolist(), cols.tolist())), worst, comp_errors)


def pair_category(system: SpinSystem, i: int, j: int, strong_hz: float = 50.0) -> str:
    """Coupling category by nucleus pair and magnitude (no bond information required)."""
    pair = "-".join(sorted((system.isotopes[i], system.isotopes[j])))
    return f"{pair}:{'strong' if abs(system.couplings_hz[i, j]) >= strong_hz else 'weak'}"


def coupling_errors_by_category(truth: SpinSystem, other: SpinSystem) -> Dict[str, List[float]]:
    aligned, match = aligned_couplings(truth, other)
    out: Dict[str, List[float]] = {}
    if aligned is None:
        return out
    t = truth.observable_couplings()
    groups = truth.group_of()
    for i in range(truth.n_spins):
        for j in range(i + 1, truth.n_spins):
            if groups[i] == groups[j]:
                continue
            out.setdefault(pair_category(truth, i, j), []).append(float(abs(aligned[i, j] - t[i, j])))
    return out
