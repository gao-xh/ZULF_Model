"""Processing parameters chosen per dataset, with the reason for every choice."""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from zulf_core.render.acquisition import Acquisition

from .diagnostics import RawDiagnostics

# Defaults that worked on every confirmed sample (configs/confirmed_samples.json "processing").
DEFAULTS = {"start_sample": 200, "stop_sample": 4200, "sg_window": 201, "sg_order": 2, "remove_mean": True,
            "apodization_rate_per_s": 0.6, "zero_fill": 4,
            "ranges": [[62.0, 119.0], [121.0, 179.0], [181.0, 293.0], [295.0, 320.0]],
            "instrument_lines_hz": [60, 120, 180, 240, 294.0, 300, 322.9]}


@dataclass
class ProcessingPlan:
    sampling_rate_hz: float
    points: int
    start_sample: int
    stop_sample: int
    sg_window: int
    sg_order: int
    remove_mean: bool
    apodization_rate_per_s: float
    zero_fill: int
    ranges: List[Tuple[float, float]]
    instrument_lines_hz: List[float]
    reasons: Dict[str, str] = field(default_factory=dict)

    def acquisition(self) -> Acquisition:
        """The zulf_core acquisition operator of this plan (the only processing implementation)."""
        return Acquisition(self.sampling_rate_hz, self.points, start_sample=self.start_sample,
                           stop_sample=self.stop_sample, sg_window=self.sg_window, sg_order=self.sg_order,
                           remove_mean=self.remove_mean, apodization_rate_per_s=self.apodization_rate_per_s)

    def to_dict(self) -> dict:
        out = asdict(self)
        out["ranges"] = [list(r) for r in self.ranges]
        return out

    @classmethod
    def from_dict(cls, d: dict) -> "ProcessingPlan":
        d = dict(d)
        d["ranges"] = [tuple(r) for r in d["ranges"]]
        return cls(**d)


def plan_for_dataset(points: int, sampling_rate_hz: float, diagnostics: Optional[RawDiagnostics] = None,
                     defaults: Optional[dict] = None, **overrides) -> ProcessingPlan:
    """Processing parameters for one dataset.

    Crop start: after the switching ringing ends (diagnostics.ringing_end_s) and never before the default; crop
    stop, SG window and order, mean removal, apodization, zero filling and ranges: the defaults (they worked on
    every confirmed sample) unless overridden. Every choice records its reason; the next optimisation step is to
    choose the window, SG and apodization from the data (SNR curve, baseline rates) as well."""
    d = dict(DEFAULTS, **(defaults or {}))
    fs = float(sampling_rate_hz)
    reasons = {}
    start = int(d["start_sample"])
    reasons["start_sample"] = f"default {start}"
    if diagnostics is not None and diagnostics.ringing_end_s:
        after = int(math.ceil(diagnostics.ringing_end_s * fs))
        if after > start:
            start = after
            reasons["start_sample"] = f"after the switching ringing ends ({diagnostics.ringing_end_s * 1e3:.2f} ms)"
        else:
            reasons["start_sample"] += f" (ringing ends earlier, {diagnostics.ringing_end_s * 1e3:.2f} ms)"
    stop = min(int(d["stop_sample"]) + (start - int(d["start_sample"])), points)
    reasons["stop_sample"] = f"default window {(int(d['stop_sample']) - int(d['start_sample'])) / fs:.3f} s"
    for key in ("sg_window", "sg_order", "remove_mean", "apodization_rate_per_s", "zero_fill", "ranges",
                "instrument_lines_hz"):
        reasons.setdefault(key, "default")
    plan = ProcessingPlan(fs, int(points), start, stop, int(d["sg_window"]), int(d["sg_order"]),
                          bool(d["remove_mean"]), float(d["apodization_rate_per_s"]), int(d["zero_fill"]),
                          [tuple(r) for r in d["ranges"]], [float(x) for x in d["instrument_lines_hz"]], reasons)
    for key, value in overrides.items():
        setattr(plan, key, value)
        plan.reasons[key] = "set by the caller"
    return plan
