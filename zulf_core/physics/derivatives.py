"""Derivatives of zero-field transition lists with respect to group couplings.

The detected signal s(t) = w Tr(U rho U^+ D), U = exp(-2 pi i H t), is an
entire function of H, so its derivative along dH/dJ = K exists everywhere,
including at degenerate levels. In the eigenbasis of H (K' = V^+ K V):

    ds = sum_ab exp(2 pi i (E_b - E_a) t) { ([r, C]_ab) d_ba + r_ab ([d, C])_ba
                                            + 2 pi i t [r, K'_D]_ab d_ba },

with r = V^+ rho V, d = V^+ D V, [x, C] = x C - C x, C_mn = K'_mn / (E_n - E_m)
for levels in different degenerate clusters (0 inside a cluster) and K'_D the
part of K' inside clusters (Daleckii-Krein). Rotations inside a cluster leave
the signal unchanged, so the result does not depend on the arbitrary basis
eigh returns for degenerate levels. For nondegenerate levels the bracket terms
are the eigenvector derivatives and the t-term is A df.

Per transition this gives a pair (dA, T) such that, for the list convention
s(t) = sum_k Re(A_k exp(2 pi i f_k t)),

    ds/dJ = sum_k Re((dA_k + 2 pi i t T_k) exp(2 pi i f_k t)).

Transitions are merged by frequency exactly as in `compute_transitions`
(amplitudes, dA and T summed). DC terms are not modeled and are dropped.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Sequence, Tuple

import numpy as np
from scipy.linalg import eigh, expm

from ..spinsystem import SpinSystem
from .operators import collective_multiplicities, product_operators
from .protocol import SUDDEN_DROP, Protocol
from .transitions import MERGE_TOLERANCE_HZ, _full_dimension, _node_structure, zeeman


@dataclass
class TransitionDerivatives:
    frequencies_hz: np.ndarray        # (T,)
    amplitudes: np.ndarray            # (T,) complex
    d_amplitudes: np.ndarray          # (P, T) complex
    t_weights: np.ndarray             # (P, T) complex; A df in the nondegenerate case
    pairs: Tuple[Tuple[int, int], ...]

    def __len__(self) -> int:
        return len(self.frequencies_hz)

    def signal_derivative(self, times_s: np.ndarray) -> np.ndarray:
        """ds/dJ_p at the given times, shape (P, len(times)) (test reference, undamped)."""
        t = np.asarray(times_s, float)
        phase = np.exp(2j * np.pi * np.outer(self.frequencies_hz, t))
        return ((self.d_amplitudes[:, :, None] + 2j * np.pi * t[None, None, :] * self.t_weights[:, :, None])
                * phase[None]).sum(axis=1).real


def _clusters(e: np.ndarray, tolerance_hz: float) -> np.ndarray:
    """Cluster label per (sorted) eigenvalue; neighbours within the tolerance share a label."""
    return np.r_[0, np.cumsum(np.diff(e) > tolerance_hz)] if len(e) else np.zeros(0, int)


def _block(h, rho, det, ks, weight, tolerance_hz, real):
    if real:
        e, v = eigh(h, driver="evr")
    else:
        e, v = eigh((h + h.conj().T) / 2, driver="evr")
    vh = v.conj().T
    r = vh @ rho @ v
    d = vh @ det @ v
    label = _clusters(e, tolerance_hz)
    same = label[:, None] == label[None, :]
    gap_nm = e[None, :] - e[:, None]          # E_n - E_m at [m, n]
    inverse = np.where(same, 0.0, 1.0 / np.where(same, 1.0, gap_nm))
    gap = gap_nm                               # E_b - E_a at [a, b]
    upper = gap > tolerance_hz
    product = r * d.T
    freqs = gap[upper]
    amps = 2 * weight * product[upper]
    d_amps, t_weights = [], []
    for k in ks:
        kp = vh @ k @ v
        c = kp * inverse
        kd = np.where(same, kp, 0.0)
        dr = r @ c - c @ r
        dd = d @ c - c @ d
        rot = dr * d.T + r * dd.T
        tw = (r @ kd - kd @ r) * d.T
        d_amps.append(2 * weight * rot[upper])
        t_weights.append(2 * weight * tw[upper])
    shape = (len(ks), len(freqs))
    return (freqs, amps.astype(complex), np.asarray(d_amps, complex).reshape(shape),
            np.asarray(t_weights, complex).reshape(shape))


def _merge(f, a, da, tw, tolerance_hz):
    order = np.argsort(f, kind="stable")
    f, a, da, tw = f[order], a[order], da[:, order], tw[:, order]
    if not len(f):
        return f, a, da, tw
    starts = np.r_[0, np.flatnonzero(np.diff(f) > tolerance_hz) + 1]
    weights = np.abs(a)
    wsum = np.add.reduceat(weights, starts)
    plain = np.add.reduceat(f, starts) / np.diff(np.r_[starts, len(f)])
    centers = np.where(wsum > 0, np.add.reduceat(weights * f, starts) / np.where(wsum > 0, wsum, 1), plain)
    return (centers, np.add.reduceat(a, starts), np.add.reduceat(da, starts, axis=1),
            np.add.reduceat(tw, starts, axis=1))


def transition_derivatives(system: SpinSystem, pairs: Sequence[Tuple[int, int]],
                           protocol: Protocol = SUDDEN_DROP, tolerance_hz: float = MERGE_TOLERANCE_HZ,
                           relative_zero: float = 1e-12) -> TransitionDerivatives:
    """Transitions and their derivatives with respect to the group couplings J_ab, (a, b) in `pairs`.

    Groups are the system's magnetic-equivalence groups (collective spins), as
    in `compute_transitions(method="sectors")`. A pair may have zero coupling.
    """
    pairs = tuple((min(a, b), max(a, b)) for a, b in pairs)
    nodes = _node_structure(system)
    if any(a == b or b >= len(nodes) for a, b in pairs):
        raise ValueError("Coupling pairs must join two distinct equivalence groups.")
    gj = system.group_couplings()
    norm = 1.0 / _full_dimension(system) if protocol.normalize_by_dimension else 1.0
    prep = [protocol.preparation_weight(sym) for sym, _, _ in nodes]
    detw = [protocol.detection_weight(sym) for sym, _, _ in nodes]
    choices = [collective_multiplicities(count, spin) for _, count, spin in nodes]
    nonzero = [(i, j) for i in range(len(nodes)) for j in range(i + 1, len(nodes)) if gj[i, j] != 0]
    m_blocks = not protocol.pulses and protocol.field_ut[0] == 0.0 and protocol.field_ut[1] == 0.0
    fs: List[np.ndarray] = []
    amps: List[np.ndarray] = []
    das: List[np.ndarray] = []
    tws: List[np.ndarray] = []
    for combo in itertools.product(*choices):
        spins = tuple(Fraction(s) for s, _ in combo)
        weight = int(np.prod([m for _, m in combo])) * norm
        site_ops, pair_ops = product_operators(spins)
        dim = site_ops[0][0].shape[0]
        h = np.zeros((dim, dim))
        for i, j in nonzero:
            h = h + gj[i, j] * pair_ops[(i, j)]
        if protocol.has_field:
            h = h + zeeman(site_ops, [sym for sym, _, _ in nodes], protocol.field_ut)
        rho = sum(p * ops[2] for p, ops in zip(prep, site_ops))
        det = sum(dw * ops[2] for dw, ops in zip(detw, site_ops))
        ks = [pair_ops[p] for p in pairs]
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
            rho, det, h = rho.real, det.real, np.ascontiguousarray(np.real(h))
        else:
            h = np.asarray(h, complex)
        if m_blocks:
            keys = np.round(2 * np.real(np.diag(sum(ops[2] for ops in site_ops)))).astype(int)
            blocks = [np.flatnonzero(keys == value) for value in np.unique(keys)]
        else:
            blocks = [np.arange(dim)]
        for idx in blocks:
            sel = np.ix_(idx, idx)
            f, a, da, tw = _block(h[sel], rho[sel], det[sel], [k[sel] for k in ks], weight, tolerance_hz, real)
            fs.append(f); amps.append(a); das.append(da); tws.append(tw)
    f = np.concatenate(fs)
    a = np.concatenate(amps)
    da = np.concatenate(das, axis=1) if das else np.zeros((len(pairs), 0), complex)
    tw = np.concatenate(tws, axis=1) if tws else np.zeros((len(pairs), 0), complex)
    f, a, da, tw = _merge(f, a, da, tw, tolerance_hz)
    scale = max(np.abs(a).max(initial=0.0), np.abs(da).max(initial=0.0), np.abs(tw).max(initial=0.0))
    keep = (np.abs(a) > relative_zero * scale) | np.any(np.abs(da) > relative_zero * scale, axis=0) | \
        np.any(np.abs(tw) > relative_zero * scale, axis=0) if scale > 0 else np.zeros(len(f), bool)
    return TransitionDerivatives(f[keep], a[keep], da[:, keep], tw[:, keep], pairs)
