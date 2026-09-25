"""Angular momentum operators and collective-spin sector bookkeeping."""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from typing import Dict, List, Sequence, Tuple

import numpy as np


def angular_momentum(spin) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Sx, Sy, Sz) for spin quantum number `spin`, basis m = S, S-1, ..., -S."""
    spin = Fraction(spin)
    dim = int(2 * spin + 1)
    m = np.array([float(spin) - k for k in range(dim)])
    plus = np.zeros((dim, dim), dtype=complex)
    s = float(spin)
    for col in range(1, dim):
        plus[col - 1, col] = np.sqrt(s * (s + 1) - m[col] * (m[col] + 1))
    minus = plus.conj().T
    return (plus + minus) / 2, (plus - minus) / (2j), np.diag(m).astype(complex)


@lru_cache(maxsize=256)
def collective_multiplicities(count: int, spin: Fraction) -> Tuple[Tuple[Fraction, int], ...]:
    """Total-spin sectors of `count` equivalent spins of quantum number `spin`.

    Returns (S, multiplicity) for every S with nonzero multiplicity. The
    multiplicity is the number of independent copies of the irreducible
    representation S, obtained from the distribution of total M.
    """
    spin = Fraction(spin)
    single = [spin - k for k in range(int(2 * spin + 1))]
    counts: Dict[Fraction, int] = {Fraction(0): 1}
    for _ in range(count):
        new: Dict[Fraction, int] = {}
        for total, c in counts.items():
            for m in single:
                new[total + m] = new.get(total + m, 0) + c
        counts = new
    top = count * spin
    out = []
    s = top
    while s >= 0:
        mult = counts.get(s, 0) - counts.get(s + 1, 0)
        if mult > 0:
            out.append((s, mult))
        s -= 1
    total_dim = sum(int(2 * s + 1) * m for s, m in out)
    expected = int(2 * spin + 1) ** count
    if total_dim != expected:
        raise AssertionError("Collective sector dimensions do not add up.")
    return tuple(out)


@lru_cache(maxsize=512)
def product_operators(spins: Tuple[Fraction, ...]):
    """Embedded (Sx, Sy, Sz) for each site of a product space with given site spins.

    Returns a list over sites of 3-tuples of complex matrices, and the list of
    real pair operators S_i . S_j for i < j keyed by (i, j).
    """
    dims = [int(2 * s + 1) for s in spins]
    site_ops: List[Tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for site, s in enumerate(spins):
        local = angular_momentum(s)
        embedded = []
        for op in local:
            full = np.ones((1, 1), dtype=complex)
            for k, d in enumerate(dims):
                full = np.kron(full, op if k == site else np.eye(d))
            embedded.append(full)
        site_ops.append(tuple(embedded))
    pairs = {}
    for i in range(len(spins)):
        for j in range(i + 1, len(spins)):
            value = sum(a @ b for a, b in zip(site_ops[i], site_ops[j]))
            pairs[(i, j)] = np.ascontiguousarray(value.real)
    return site_ops, pairs


def full_space_spins(spins: Sequence[Fraction]) -> Tuple[Fraction, ...]:
    return tuple(Fraction(s) for s in spins)
