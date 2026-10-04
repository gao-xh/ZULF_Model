"""Line table: every transition of a spin system with the couplings that place it.

The zero-field Hamiltonian is linear and homogeneous in the couplings, H = sum_k J_k H_k (H_k = 2 pi I_a . I_b for
the group pair k). Scaling every coupling by s scales every eigenvalue and so every line frequency by s: each
frequency is homogeneous of degree 1 in the couplings, and Euler's theorem gives the exact decomposition

    f = sum_k J_k df/dJ_k,

where df/dJ_k = <a|H_k|a> - <b|H_k|b> (Hellmann-Feynman, `transition_derivatives`). The term J_k df/dJ_k is the
contribution of coupling k to the position of the line (in Hz); df/dJ_k also predicts how the line moves,
f(J + dJ) ~ f(J) + sum_k df/dJ_k dJ_k, to first order. For a line of a weakly coupled heteronuclear group the
coefficient of its one-bond coupling is a simple number (1 for X-H, 3/2 for X-H2, 1 or 2 for X-H3) and the small
couplings enter with their own coefficients.

Second order (`second_order=True`): the coefficients df/dJ change with J because the eigenstates do, so f is not
linear in J. d2f/dJ_k dJ_l (cross terms included) is taken by central differences of the analytic df/dJ. They are
large where levels nearly coincide (second-order perturbation: 1/(E_a - E_n)); there a Taylor step is only valid for
changes much smaller than the gap and the spectrum must be recomputed. Two indicators per line: `nearest_line_hz`
(the closest other line of the same system) and `trust_step_hz`, the coupling change for which the largest
second-order term reaches `tolerance_hz`. Homogeneity also gives sum_l J_l d2f/dJ_k dJ_l = 0 (checked in the tests).

Exactly degenerate transitions are merged by `transition_derivatives`; their df/dJ is the amplitude-weighted mean.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..spinsystem import SpinSystem
from .derivatives import transition_derivatives
from .protocol import SUDDEN_DROP, Protocol


@dataclass
class Line:
    frequency_hz: float
    amplitude: complex
    relative_amplitude: float
    df_dj: Dict[str, float]                 # coupling name -> df/dJ (dimensionless)
    contribution_hz: Dict[str, float] = field(default_factory=dict)   # name -> J df/dJ; sums to frequency_hz
    hessian: Optional[Dict[Tuple[str, str], float]] = None             # (name, name) -> d2f/dJ dJ (1/Hz)
    nearest_line_hz: Optional[float] = None                            # distance to the closest other line
    trust_step_hz: Optional[float] = None                              # largest 2nd-order term = tolerance_hz

    def formula(self, digits: int = 3, minimum: float = 1e-3) -> str:
        """f = sum c_k J_k with the coefficients c_k = df/dJ_k at this point (terms below `minimum` dropped)."""
        terms = [f"{c:+.{digits}f} {n}" for n, c in sorted(self.df_dj.items(), key=lambda kv: -abs(kv[1]))
                 if abs(c) >= minimum]
        return f"{self.frequency_hz:.3f} Hz = " + " ".join(terms)

    def near_degenerate(self, threshold_hz: float = 0.1) -> bool:
        """Another line of the same system lies within threshold_hz: its eigenstates mix strongly and the Taylor
        expansion (any order) holds only for coupling changes well below that spacing; recompute instead."""
        return self.nearest_line_hz is not None and self.nearest_line_hz < threshold_hz

    def shifted(self, dj: Dict[str, float], second_order: bool = False) -> float:
        """Frequency after changing couplings by dj (Hz): first order, plus the second-order and cross terms
        0.5 sum_kl d2f/dJk dJl dJk dJl when second_order (needs the Hessian)."""
        f = self.frequency_hz + sum(self.df_dj.get(n, 0.0) * v for n, v in dj.items())
        if second_order:
            if self.hessian is None:
                raise ValueError("No Hessian: build the table with second_order=True.")
            f += 0.5 * sum(self.hessian.get((a, b), 0.0) * va * vb for a, va in dj.items() for b, vb in dj.items())
        return f


def _with_group_change(system: SpinSystem, pairs: Sequence[Tuple[int, int]], dj: float) -> SpinSystem:
    j = np.array(system.couplings_hz, float)
    groups = system.groups
    for a, b in pairs:
        for i in groups[a]:
            for k in groups[b]:
                j[i, k] += dj
                j[k, i] += dj
    return SpinSystem(system.isotopes, j, groups)


def line_table(system: SpinSystem, names: Optional[Dict[Tuple[int, int], str]] = None,
               protocol: Protocol = SUDDEN_DROP, band_hz: Optional[Tuple[float, float]] = None,
               min_relative_amplitude: float = 0.0, second_order: bool = False, step_hz: float = 1e-3,
               tolerance_hz: float = 0.02) -> List[Line]:
    """Every transition of `system` (optionally only inside band_hz and above min_relative_amplitude of the largest)
    with df/dJ and the exact contribution J df/dJ of every group coupling (zero ones included). `names` maps a group
    pair (a, b), a < b, to a label (default "J(a,b)"); pairs with the same label are summed (tied couplings)."""
    gj = system.group_couplings()
    g = gj.shape[0]
    # every group pair, zero couplings included: a coupling held at 0 still moves lines (its df/dJ is reported, its
    # contribution is 0); pairs sharing a label (symmetry-tied couplings) are summed into it
    pairs = [(a, b) for a in range(g) for b in range(a + 1, g)]
    if not pairs:
        return []
    label = {p: (names or {}).get(p, f"J({p[0]},{p[1]})") for p in pairs}
    d = transition_derivatives(system, pairs, protocol)
    amp = np.abs(d.amplitudes)
    top = float(amp.max()) if len(amp) else 1.0
    out = []
    for i in range(len(d)):
        f = float(d.frequencies_hz[i])
        if band_hz is not None and not (band_hz[0] <= f <= band_hz[1]):
            continue
        if amp[i] < min_relative_amplitude * top or amp[i] == 0.0:
            continue
        dfdj, contrib = {}, {}
        for k, p in enumerate(pairs):
            v = float((d.t_weights[k, i] / d.amplitudes[i]).real)
            dfdj[label[p]] = dfdj.get(label[p], 0.0) + v
            contrib[label[p]] = contrib.get(label[p], 0.0) + float(gj[p]) * v
        out.append(Line(f, complex(d.amplitudes[i]), float(amp[i] / top), dfdj, contrib))
    out.sort(key=lambda line: line.frequency_hz)
    every = np.sort(np.asarray(d.frequencies_hz, float))
    for line in out:
        others = every[np.abs(every - line.frequency_hz) > 1e-9]
        line.nearest_line_hz = float(np.min(np.abs(others - line.frequency_hz))) if len(others) else None
    if second_order and out:
        groups_of = {}
        for p in pairs:
            groups_of.setdefault(label[p], []).append(p)
        keys = list(groups_of)
        cols = {}
        for key in keys:
            plus = line_table(_with_group_change(system, groups_of[key], step_hz), names, protocol)
            minus = line_table(_with_group_change(system, groups_of[key], -step_hz), names, protocol)
            cols[key] = (plus, minus)
        for line in out:
            hess = {}
            for kb in keys:
                plus, minus = cols[kb]
                lp = min(plus, key=lambda q: (abs(q.frequency_hz - line.frequency_hz),
                                              abs(q.relative_amplitude - line.relative_amplitude)))
                lm = min(minus, key=lambda q: (abs(q.frequency_hz - line.frequency_hz),
                                               abs(q.relative_amplitude - line.relative_amplitude)))
                for ka in keys:
                    hess[(ka, kb)] = (lp.df_dj.get(ka, 0.0) - lm.df_dj.get(ka, 0.0)) / (2 * step_hz)
            for ka in keys:                      # symmetrize (exact Hessians are symmetric)
                for kb in keys:
                    if keys.index(kb) > keys.index(ka):
                        v = 0.5 * (hess[(ka, kb)] + hess[(kb, ka)])
                        hess[(ka, kb)] = hess[(kb, ka)] = v
            line.hessian = hess
            top = max(abs(v) for v in hess.values())
            line.trust_step_hz = float(np.sqrt(2 * tolerance_hz / top)) if top > 0 else float("inf")
    return out


def predict_lines(lines: Sequence[Line], dj: Dict[str, float]) -> np.ndarray:
    """First-order line positions after the coupling changes dj (Hz), from a line table."""
    return np.array([line.shifted(dj) for line in lines])
