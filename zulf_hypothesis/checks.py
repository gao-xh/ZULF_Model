"""Physical checks on refined hypotheses, as a registry of small rules.

Each check reads a refinement summary (`RefinementResult.summary()` or the
saved JSON of one), the built `HypothesisModel` and, optionally, peer
results of the same hypothesis (other starts or variants), and returns
findings. Add a rule with `@register_check`; `run_checks` applies all of
them. The rules encode skills/zulf-hypothesis-refinement/references/checks.md.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Mapping, Optional, Sequence

import numpy as np

from .builder import HypothesisModel


@dataclass
class Finding:
    code: str
    severity: str          # "reject" | "warn" | "info"
    message: str
    data: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"code": self.code, "severity": self.severity, "message": self.message, "data": self.data}


@dataclass
class CheckContext:
    model: HypothesisModel
    result: Mapping                     # refinement summary
    peers: Sequence[Mapping] = ()       # other results of the same hypothesis
    options: Mapping = field(default_factory=dict)

    @property
    def amplitudes(self) -> np.ndarray:
        return np.array([abs(complex(*g)) if isinstance(g, (list, tuple)) else abs(g) for g in self.result["gains"]])

    def rates(self) -> np.ndarray:
        p = self.result["parameters"]
        return np.array([math.exp(p[f"c{c}.log_rate0"]) for c in range(len(self.model.component_labels))
                         if f"c{c}.log_rate0" in p])


CHECKS: Dict[str, Callable[[CheckContext], List[Finding]]] = {}


def register_check(code: str):
    def wrap(fn):
        CHECKS[code] = fn
        return fn
    return wrap


def run_checks(model: HypothesisModel, result: Mapping, peers: Sequence[Mapping] = (),
               only: Optional[Sequence[str]] = None, **options) -> List[Finding]:
    ctx = CheckContext(model, result, peers, options)
    out: List[Finding] = []
    for code, fn in CHECKS.items():
        if only is None or code in only:
            out.extend(fn(ctx))
    return out


@register_check("bounds")
def _bounds(ctx: CheckContext) -> List[Finding]:
    hits = list(ctx.result.get("boundary_hits") or [])
    return [Finding("bounds", "warn", f"parameters at bounds: {hits}", {"parameters": hits})] if hits else []


@register_check("abundance")
def _abundance(ctx: CheckContext) -> List[Finding]:
    """Free amplitudes against the natural-abundance ratios of the built model."""
    if not ctx.model.abundance_known:
        return []
    amp = ctx.amplitudes
    expected = np.asarray(ctx.model.ratios, float)
    if len(amp) != len(expected) or len(amp) < 2 or amp.max() <= 0:
        return []
    tol = float(ctx.options.get("abundance_factor", 1.6))
    observed = amp / amp[0]
    out = []
    for c, (o, e) in enumerate(zip(observed, expected)):
        if c == 0:
            continue
        if o < e / 20:
            out.append(Finding("abundance", "warn", f"{ctx.model.component_labels[c]} collapsed: {o:.3g} vs {e:.3g}",
                               {"component": c, "observed": o, "expected": e}))
        elif not (e / tol <= o <= e * tol):
            out.append(Finding("abundance", "warn", f"{ctx.model.component_labels[c]} amplitude {o:.3g} vs natural "
                               f"{e:.3g} (outside x{tol})", {"component": c, "observed": o, "expected": e}))
    return out


@register_check("background_component")
def _background_component(ctx: CheckContext) -> List[Finding]:
    """A component much stronger than natural abundance allows and much broader than its partners."""
    if not ctx.model.abundance_known:
        return []
    amp, rates = ctx.amplitudes, ctx.rates()
    expected = np.asarray(ctx.model.ratios, float)
    if len(amp) < 2 or len(rates) != len(amp):
        return []
    per_site = amp / expected
    out = []
    for c in range(len(amp)):
        scale = np.median(np.delete(per_site, c))
        if per_site[c] > 3 * scale and rates[c] > 2 * np.median(np.delete(rates, c)):
            out.append(Finding("background_component", "warn",
                               f"{ctx.model.component_labels[c]} acts as background (amplitude "
                               f"{per_site[c] / scale:.1f}x the others, rate {rates[c]:.2g} 1/s)",
                               {"component": c}))
    return out


@register_check("rate_asymmetry")
def _rate_asymmetry(ctx: CheckContext) -> List[Finding]:
    """Isotopologues of one molecule should decay alike; a much faster one hides couplings."""
    rates = ctx.rates()
    if len(rates) < 2 or rates.min() <= 0:
        return []
    factor = float(ctx.options.get("rate_factor", 2.0))
    if rates.max() / rates.min() < factor:
        return []
    c = int(np.argmax(rates))
    return [Finding("rate_asymmetry", "info", f"{ctx.model.component_labels[c]} decays {rates.max() / rates.min():.1f}x "
                    "faster than the slowest isotopologue: unresolved couplings on it?",
                    {"component": c, "rates": rates.tolist()})]


@register_check("undetermined")
def _undetermined(ctx: CheckContext) -> List[Finding]:
    """Named couplings that differ between comparably good peer results."""
    if not ctx.peers:
        return []
    key = ctx.options.get("residual_key", "data_region_residual")
    best = ctx.result.get(key)
    near = [p for p in ctx.peers if best is None or p.get(key) is None or p[key] <= 1.1 * best + 1e-12]
    if not near:
        return []
    tol = float(ctx.options.get("coupling_tolerance_hz", 0.5))
    mine = ctx.model.named_couplings(ctx.result["parameters"])
    spread = {}
    for k, v in mine.items():
        values = [v] + [ctx.model.named_couplings(p["parameters"]).get(k, v) for p in near]
        if max(values) - min(values) > tol:
            spread[k] = [round(x, 3) for x in values]
    if not spread:
        return []
    return [Finding("undetermined", "info", f"couplings not determined across comparable fits: {sorted(spread)}",
                    {"spread": spread})]


@register_check("misfit")
def _misfit(ctx: CheckContext) -> List[Finding]:
    """Reduced chi2 on the common yardstick (set by the search as summary['reduced_chi2']) well above 1:
    the data hold structure the hypothesis does not explain."""
    value = ctx.result.get("reduced_chi2")
    limit = float(ctx.options.get("misfit_limit", 3.0))
    if value is None or value <= limit:
        return []
    return [Finding("misfit", "info", f"reduced chi2 {value:.1f} on the data cores (> {limit:g}): unexplained structure",
                    {"reduced_chi2": value})]
