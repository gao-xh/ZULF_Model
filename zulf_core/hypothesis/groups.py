"""Group candidates: which X-Hn group, with which 1J, could produce each band.

An isolated XHn group has lines at fixed multiples of J (XH: J; XH2: 3J/2;
XH3: J and 2J, 0.8 : 1), computed once with `compute_transitions` and
scaled. For every band peak and every pattern line, the hypothesis "this peak
is that line" fixes J; the other lines of the pattern are then checked in the
data: observed (support), unobservable (outside the fitted ranges or on an
instrument line: neutral) or missing (penalised). Candidates record which
bands their lines explain, so a 2J partner band is not taken for a separate
group. Patterns are a registry (`register_pattern`) for new nuclei or groups.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..physics import compute_transitions
from ..spinsystem import SpinSystem
from .inventory import Inventory


@dataclass(frozen=True)
class GroupPattern:
    name: str               # e.g. "13CH3"
    nucleus: str
    n_h: int
    j_range_hz: Tuple[float, float]
    prior: float = 1.0      # search weight (15N: lower abundance, N-H usually exchange-decoupled)

    @property
    def unit_lines(self) -> Tuple[Tuple[float, float], ...]:
        return _unit_lines(self.nucleus, self.n_h)


@lru_cache(maxsize=None)
def _unit_lines(nucleus: str, n_h: int) -> Tuple[Tuple[float, float], ...]:
    """(frequency / J, relative amplitude) of an isolated X-Hn group."""
    s = SpinSystem.from_group_couplings([nucleus, "1H"], [1, n_h], np.array([[0, 100.0], [100.0, 0]]))
    tl = compute_transitions(s)
    a = np.abs(tl.amplitudes)
    keep = a >= 0.05 * a.max()
    return tuple((float(f / 100.0), float(x / a.max())) for f, x in zip(tl.frequencies_hz[keep], a[keep]))


PATTERNS: Dict[str, GroupPattern] = {}


def register_pattern(pattern: GroupPattern) -> GroupPattern:
    PATTERNS[pattern.name] = pattern
    return pattern


for _n, _name in ((1, "13CH"), (2, "13CH2"), (3, "13CH3")):
    register_pattern(GroupPattern(_name, "13C", _n, (115.0, 230.0)))
for _n, _name in ((1, "15NH"), (2, "15NH2"), (3, "15NH3")):
    register_pattern(GroupPattern(_name, "15N", _n, (60.0, 100.0), prior=0.6))


@dataclass
class GroupCandidate:
    pattern: str
    nucleus: str
    n_h: int
    j_hz: float
    anchor_band: int
    lines: List[Tuple[float, float, str, float]]   # (frequency, relative amplitude, status, support sigma)
    explains: Tuple[int, ...]                        # bands containing a predicted line
    score: float
    prior: float = 1.0
    hint_support: float = 0.0
    sources: List[str] = field(default_factory=lambda: ["data"])
    notes: List[str] = field(default_factory=list)

    @property
    def rank_score(self) -> float:
        return self.score * self.prior + self.hint_support

    @property
    def support(self) -> float:
        return float(sum(sp for _, _, st, sp in self.lines if st != "unobservable"))

    def sort_key(self):
        return (-self.rank_score, -len(self.explains), -self.support, self.n_h)

    def to_dict(self) -> dict:
        return {"pattern": self.pattern, "j_hz": round(self.j_hz, 3), "anchor_band": self.anchor_band,
                "explains": list(self.explains), "score": round(self.score, 3), "prior": self.prior,
                "hint_support": round(self.hint_support, 3), "sources": self.sources,
                "lines": [(round(f, 2), round(a, 2), st, round(sp, 1)) for f, a, st, sp in self.lines],
                "notes": self.notes}


def group_candidates(inventory: Inventory, patterns: Optional[Sequence[str]] = None, tol_hz: float = 1.5,
                     partner_threshold: float = 2.0, min_score: float = 0.3,
                     explain_fraction: float = 0.4, multiplet_hz: float = 10.0) -> List[GroupCandidate]:
    """Candidates for every band, best first within each band.

    score = (sum of relative amplitudes of observed lines, weighted by
    min(support / need, 1)) / (sum over observable lines), where `need` is
    `partner_threshold` for the anchor and, for partner lines, also half the
    height predicted from the anchor (isolated-group ratios; coupled groups
    can have weaker partners, so this only lowers, never removes).
    A predicted line explains another band only when its height, scaled from
    the anchor line, reaches `explain_fraction` of that band's peak (a weak
    partner cannot account for a strong band). Small couplings split a group's
    lines into multiplets up to about `multiplet_hz` wide, so a line explains
    every band within that distance.
    """
    chosen = [PATTERNS[p] for p in (patterns or PATTERNS)]
    out: List[GroupCandidate] = []
    for band in inventory.bands:
        anchors = sorted({round(l.frequency_hz, 3) for l in band.lines}, key=lambda f: -band_line_snr(band, f))[:3]
        seen = set()
        for pattern in chosen:
            for f_anchor in anchors:
                for ratio, anchor_amp in pattern.unit_lines:
                    j = f_anchor / ratio
                    scale = band_line_snr(band, f_anchor) / anchor_amp
                    if not (pattern.j_range_hz[0] <= j <= pattern.j_range_hz[1]):
                        continue
                    key = (pattern.name, round(j, 1))
                    if key in seen:
                        continue
                    seen.add(key)
                    lines, got, possible, explains = [], 0.0, 0.0, set()
                    for r, amp in pattern.unit_lines:
                        f = r * j
                        if not inventory.observable(f):
                            lines.append((f, amp, "unobservable", 0.0))
                            continue
                        sup = inventory.support(f, tol_hz)
                        # A partner must reach half its height predicted from the anchor (or the threshold).
                        need = partner_threshold if abs(f - f_anchor) < 1e-9 else max(partner_threshold,
                                                                                    0.5 * scale * amp)
                        weight = min(sup / need, 1.0)
                        possible += amp
                        got += amp * weight
                        status = "observed" if weight >= 1.0 else ("weak" if weight > 0 else "missing")
                        lines.append((f, amp, status, sup))
                        if weight > 0:
                            for b in inventory.bands:
                                if b.contains(f, multiplet_hz) and (b.index == band.index or
                                                                    scale * amp >= explain_fraction * b.peak.snr):
                                    explains.add(b.index)
                    if possible <= 0:
                        continue
                    score = got / possible
                    if score < min_score:
                        continue
                    out.append(GroupCandidate(pattern.name, pattern.nucleus, pattern.n_h, float(j), band.index,
                                              lines, tuple(sorted(explains)), float(score), pattern.prior))
    out.sort(key=lambda c: (c.anchor_band,) + c.sort_key())
    return out


def band_line_snr(band, f: float) -> float:
    return max((l.snr for l in band.lines if abs(l.frequency_hz - f) < 1e-3), default=0.0)


def by_band(candidates: Sequence[GroupCandidate]) -> Dict[int, List[GroupCandidate]]:
    out: Dict[int, List[GroupCandidate]] = {}
    for c in candidates:
        out.setdefault(c.anchor_band, []).append(c)
    for v in out.values():
        v.sort(key=lambda c: c.sort_key())
    return out
