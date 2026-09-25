"""Zero-field transition lists (frequency, complex amplitude) for spin systems.

The detected real signal is `s(t) = dc + sum_k Re(A_k exp(2 pi i f_k t))` with
positive `f_k` in Hz. Amplitudes are per molecule when the protocol normalizes
by the full Hilbert dimension. Transition lists do not depend on linewidth and
can be cached and re-rendered.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.linalg import eigh, expm

from ..nuclei import get_registry
from ..spinsystem import SpinSystem
from .operators import collective_multiplicities, product_operators
from .protocol import SUDDEN_DROP, Protocol

MERGE_TOLERANCE_HZ = 1e-7
NUMERICAL_ZERO_RELATIVE = 1e-12


@dataclass
class TransitionList:
    frequencies_hz: np.ndarray
    amplitudes: np.ndarray
    dc: complex = 0.0
    families: Optional[np.ndarray] = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.frequencies_hz = np.asarray(self.frequencies_hz, float)
        self.amplitudes = np.asarray(self.amplitudes, complex)
        if self.frequencies_hz.shape != self.amplitudes.shape or self.frequencies_hz.ndim != 1:
            raise ValueError("Frequencies and amplitudes must be matching 1D arrays.")
        if np.any(self.frequencies_hz <= 0):
            raise ValueError("Transition frequencies must be positive; DC is stored separately.")
        if self.families is None:
            self.families = np.zeros(len(self.frequencies_hz), dtype=int)
        self.families = np.asarray(self.families, int)

    def __len__(self) -> int:
        return len(self.frequencies_hz)

    def signal(self, times_s: np.ndarray) -> np.ndarray:
        """Noiseless undamped real signal at the given times (test reference)."""
        t = np.asarray(times_s, float)
        phase = np.exp(2j * np.pi * t[:, None] * self.frequencies_hz[None, :])
        return (phase @ self.amplitudes).real + np.real(self.dc)

    def scaled(self, factor: complex) -> "TransitionList":
        return TransitionList(self.frequencies_hz.copy(), self.amplitudes * factor, self.dc * factor,
                              self.families.copy(), dict(self.metadata))

    def with_families(self, families: np.ndarray) -> "TransitionList":
        return TransitionList(self.frequencies_hz.copy(), self.amplitudes.copy(), self.dc,
                              np.asarray(families, int), dict(self.metadata))

    def split_families(self, edges_hz: Sequence[float]) -> "TransitionList":
        """Label transitions by frequency interval; a transition on an edge goes up."""
        return self.with_families(np.searchsorted(np.asarray(edges_hz, float), self.frequencies_hz, side="right"))

    def within(self, f_min: float, f_max: float) -> "TransitionList":
        keep = (self.frequencies_hz >= f_min) & (self.frequencies_hz <= f_max)
        return TransitionList(self.frequencies_hz[keep], self.amplitudes[keep], self.dc,
                              self.families[keep], dict(self.metadata))

    @staticmethod
    def concatenate(parts: Sequence["TransitionList"]) -> "TransitionList":
        if not parts:
            return TransitionList(np.zeros(0), np.zeros(0, complex))
        return TransitionList(np.concatenate([p.frequencies_hz for p in parts]),
                              np.concatenate([p.amplitudes for p in parts]),
                              complex(sum(p.dc for p in parts)),
                              np.concatenate([p.families for p in parts]))

    def to_dict(self) -> dict:
        return {"frequencies_hz": self.frequencies_hz.tolist(),
                "amplitudes_real": self.amplitudes.real.tolist(),
                "amplitudes_imag": self.amplitudes.imag.tolist(),
                "dc": [float(np.real(self.dc)), float(np.imag(self.dc))],
                "families": self.families.tolist(), "metadata": self.metadata}

    @classmethod
    def from_dict(cls, data: dict) -> "TransitionList":
        return cls(np.asarray(data["frequencies_hz"]),
                   np.asarray(data["amplitudes_real"]) + 1j * np.asarray(data["amplitudes_imag"]),
                   complex(*data.get("dc", [0.0, 0.0])), np.asarray(data.get("families"), int)
                   if data.get("families") is not None else None, dict(data.get("metadata", {})))


def merge_transitions(frequencies: np.ndarray, amplitudes: np.ndarray,
                      tolerance_hz: float = MERGE_TOLERANCE_HZ,
                      relative_zero: float = NUMERICAL_ZERO_RELATIVE):
    """Merge coincident frequencies by summing amplitudes; drop numerical zeros.

    The merged frequency is the amplitude-magnitude weighted mean of the cluster.
    `relative_zero` removes only numerically vanishing amplitudes (relative to
    the largest), not physically weak lines.
    """
    f = np.asarray(frequencies, float)
    a = np.asarray(amplitudes, complex)
    if not len(f):
        return f, a
    order = np.argsort(f, kind="stable")
    f, a = f[order], a[order]
    breaks = np.flatnonzero(np.diff(f) > tolerance_hz) + 1
    starts = np.r_[0, breaks]
    summed = np.add.reduceat(a, starts)
    weights = np.abs(a)
    wsum = np.add.reduceat(weights, starts)
    fsum = np.add.reduceat(weights * f, starts)
    plain = np.add.reduceat(f, starts) / np.diff(np.r_[starts, len(f)])
    centers = np.where(wsum > 0, fsum / np.where(wsum > 0, wsum, 1), plain)
    scale = np.abs(summed).max() if len(summed) else 0.0
    keep = np.abs(summed) > relative_zero * scale if scale > 0 else np.zeros(len(summed), bool)
    return centers[keep], summed[keep]


def _node_structure(system: SpinSystem):
    """Collective nodes: one per equivalence group, with nucleus and member count."""
    registry = get_registry()
    nodes = []
    for group in system.groups:
        symbol = system.isotopes[group[0]]
        nodes.append((symbol, len(group), registry.spin(symbol)))
    return nodes


def _full_dimension(system: SpinSystem) -> int:
    registry = get_registry()
    dim = 1
    for s in system.isotopes:
        dim *= registry[s].multiplicity
    return dim


def zeeman(site_ops, symbols: Sequence[str], field_ut) -> np.ndarray:
    """Zeeman Hamiltonian in Hz, -sum_n gamma_n B . I_n, for sites of the given nuclei."""
    registry = get_registry()
    h = 0
    for sym, ops in zip(symbols, site_ops):
        g = registry.gamma(sym)
        for b, op in zip(field_ut, ops):
            if b != 0.0:
                h = h - g * b * op
    return h


def _block_amplitudes(h: np.ndarray, rho: np.ndarray, det: np.ndarray, weight: float,
                      tolerance_hz: float, real: bool):
    if real:
        e, v = eigh(h, driver="evr")
    else:
        e, v = eigh((h + h.conj().T) / 2, driver="evr")
    r = v.conj().T @ rho @ v
    d = v.conj().T @ det @ v
    product = r * d.T  # rho_ab D_ba
    gap = e[None, :] - e[:, None]  # E_b - E_a
    upper = gap > tolerance_hz
    degenerate = np.abs(gap) <= tolerance_hz
    dc = weight * product[degenerate].sum()
    return gap[upper], 2 * weight * product[upper], dc


def compute_transitions(system: SpinSystem, protocol: Protocol = SUDDEN_DROP, method: str = "sectors",
                        tolerance_hz: float = MERGE_TOLERANCE_HZ,
                        relative_zero: float = NUMERICAL_ZERO_RELATIVE, m_blocking: bool = True) -> TransitionList:
    """Transition list of one spin system under a protocol.

    `method="sectors"` combines each magnetic-equivalence group into collective
    spins and diagonalizes each product of total-spin sectors separately;
    `method="full"` diagonalizes the full product space (reference, small n).
    With longitudinal preparation and detection and no transverse field or
    pulse, each sector is further split by total M (`m_blocking`), which is
    exact because coherences between different M carry no amplitude.
    """
    if method == "sectors":
        nodes = _node_structure(system)
        gj = system.group_couplings()
    elif method == "full":
        registry = get_registry()
        nodes = [(s, 1, registry.spin(s)) for s in system.isotopes]
        gj = system.couplings_hz.copy()
    else:
        raise ValueError("method must be 'sectors' or 'full'.")
    full_dim = _full_dimension(system)
    norm = 1.0 / full_dim if protocol.normalize_by_dimension else 1.0
    prep = [protocol.preparation_weight(sym) for sym, _, _ in nodes]
    detw = [protocol.detection_weight(sym) for sym, _, _ in nodes]
    choices = [collective_multiplicities(count, spin) for _, count, spin in nodes]
    freqs: List[np.ndarray] = []
    amps: List[np.ndarray] = []
    dc = 0.0 + 0.0j
    pair_index = [(i, j) for i in range(len(nodes)) for j in range(i + 1, len(nodes)) if gj[i, j] != 0]
    m_blocks = m_blocking and not protocol.pulses and protocol.field_ut[0] == 0.0 and protocol.field_ut[1] == 0.0
    for combo in itertools.product(*choices):
        spins = tuple(Fraction(s) for s, _ in combo)
        multiplicity = int(np.prod([m for _, m in combo]))
        site_ops, pairs = product_operators(spins)
        dim = site_ops[0][0].shape[0]
        h = np.zeros((dim, dim))
        for i, j in pair_index:
            h += gj[i, j] * pairs[(i, j)]
        if protocol.has_field:
            h = h + zeeman(site_ops, [sym for sym, _, _ in nodes], protocol.field_ut)
        rho = sum(p * ops[2] for p, ops in zip(prep, site_ops))
        det = sum(d * ops[2] for d, ops in zip(detw, site_ops))
        real = protocol.is_real
        for pulse in protocol.pulses:
            axis = 0 if pulse.axis == "x" else 1
            generator = sum(ops[axis] for (sym, _, _), ops in zip(nodes, site_ops) if sym in pulse.nuclei)
            if isinstance(generator, int):
                continue
            u = expm(-1j * pulse.angle_rad * generator)
            rho = u @ rho @ u.conj().T
            real = False
        if real:
            rho, det = rho.real, det.real
            h_use = np.ascontiguousarray(np.real(h))
        else:
            h_use = np.asarray(h, complex)
        if m_blocks:
            # Longitudinal preparation and detection are diagonal in the product basis and the
            # Hamiltonian conserves total M (zero field or a field along z), so coherences between
            # different M carry no amplitude: diagonalize each M block separately (exact).
            total_m = np.real(np.diag(sum(ops[2] for ops in site_ops)))
            keys = np.round(2 * total_m).astype(int)
            for value in np.unique(keys):
                idx = np.flatnonzero(keys == value)
                block = np.ix_(idx, idx)
                f, a, block_dc = _block_amplitudes(h_use[block], rho[block], det[block], multiplicity * norm,
                                                   tolerance_hz, real)
                freqs.append(f)
                amps.append(a.astype(complex))
                dc += block_dc
            continue
        f, a, block_dc = _block_amplitudes(h_use, rho, det, multiplicity * norm, tolerance_hz, real)
        freqs.append(f)
        amps.append(a.astype(complex))
        dc += block_dc
    f, a = merge_transitions(np.concatenate(freqs), np.concatenate(amps), tolerance_hz, relative_zero)
    return TransitionList(f, a, complex(dc), metadata={"method": method, "protocol": protocol.name,
                                                        "n_spins": system.n_spins})


def reference_signal(system: SpinSystem, times_s: np.ndarray, protocol: Protocol = SUDDEN_DROP) -> np.ndarray:
    """Brute-force full-space propagation with matrix exponentials (tests only)."""
    registry = get_registry()
    spins = tuple(registry.spin(s) for s in system.isotopes)
    site_ops, pairs = product_operators(spins)
    dim = site_ops[0][0].shape[0]
    h = np.zeros((dim, dim), complex)
    j = system.couplings_hz
    for (i, k), op in pairs.items():
        h += j[i, k] * op
    if protocol.has_field:
        h = h + zeeman(site_ops, list(system.isotopes), protocol.field_ut)
    rho = sum(protocol.preparation_weight(s) * ops[2] for s, ops in zip(system.isotopes, site_ops))
    det = sum(protocol.detection_weight(s) * ops[2] for s, ops in zip(system.isotopes, site_ops))
    for pulse in protocol.pulses:
        axis = 0 if pulse.axis == "x" else 1
        generator = sum(ops[axis] for s, ops in zip(system.isotopes, site_ops) if s in pulse.nuclei)
        if isinstance(generator, int):
            continue
        u = expm(-1j * pulse.angle_rad * generator)
        rho = u @ rho @ u.conj().T
    norm = 1.0 / dim if protocol.normalize_by_dimension else 1.0
    out = np.empty(len(times_s))
    for k, t in enumerate(times_s):
        u = expm(-2j * np.pi * h * t)
        out[k] = np.real(np.trace(u @ rho @ u.conj().T @ det)) * norm
    return out


class TransitionCache:
    """Small exact-key cache of transition lists keyed by system content and protocol."""

    def __init__(self, max_entries: int = 4096):
        self.max_entries = max_entries
        self._store: Dict[tuple, TransitionList] = {}
        self.hits = 0
        self.misses = 0

    @staticmethod
    def key(system: SpinSystem, protocol: Protocol) -> tuple:
        return (system.isotopes, system.couplings_hz.tobytes(), system.groups, repr(protocol.to_dict()))

    def get(self, system: SpinSystem, protocol: Protocol = SUDDEN_DROP) -> TransitionList:
        k = self.key(system, protocol)
        if k in self._store:
            self.hits += 1
            return self._store[k]
        self.misses += 1
        value = compute_transitions(system, protocol)
        if len(self._store) >= self.max_entries:
            self._store.pop(next(iter(self._store)))
        self._store[k] = value
        return value
