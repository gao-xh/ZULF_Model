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

    def formula(self, digits: int = 3, minimum: float = 1e-3) -> str:
        """f = sum c_k J_k with the coefficients c_k = df/dJ_k at this point (terms below `minimum` dropped)."""
        terms = [f"{c:+.{digits}f} {n}" for n, c in sorted(self.df_dj.items(), key=lambda kv: -abs(kv[1]))
                 if abs(c) >= minimum]
        return f"{self.frequency_hz:.3f} Hz = " + " ".join(terms)

    def shifted(self, dj: Dict[str, float]) -> float:
        """First-order frequency after changing couplings by dj (Hz)."""
        return self.frequency_hz + sum(self.df_dj.get(n, 0.0) * v for n, v in dj.items())


def line_table(system: SpinSystem, names: Optional[Dict[Tuple[int, int], str]] = None,
               protocol: Protocol = SUDDEN_DROP, band_hz: Optional[Tuple[float, float]] = None,
               min_relative_amplitude: float = 0.0) -> List[Line]:
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
    return out


def predict_lines(lines: Sequence[Line], dj: Dict[str, float]) -> np.ndarray:
    """First-order line positions after the coupling changes dj (Hz), from a line table."""
    return np.array([line.shifted(dj) for line in lines])
