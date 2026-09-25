"""Throughput measurement for transition lists and rendered samples (Phase 0 table)."""
from __future__ import annotations

import platform
import time

import numpy as np

from .generator import build_default_sampler
from .physics import compute_transitions
from .render import ProcessingConfig, SampleRenderer
from .spec import ProblemSpec
from .timing import Timer


def measure(samples: int = 50, points: int = 16384, spec: ProblemSpec = None) -> dict:
    spec = spec or ProblemSpec()
    sampler = build_default_sampler(spec)
    drawn = list(sampler.generate(samples, seed=0))
    systems = [c.system for s in drawn for c in s.interpretation.components]
    start = time.perf_counter()
    counts = [len(compute_transitions(s)) for s in systems]
    per_system = (time.perf_counter() - start) / max(1, len(systems))
    timer = Timer()
    renderer = SampleRenderer(spec, ProcessingConfig(mode="randomized", points=points), timer=timer)
    rng = np.random.default_rng(0)
    start = time.perf_counter()
    for s in drawn:
        renderer.render(s.interpretation, rng)
    per_sample = (time.perf_counter() - start) / max(1, len(drawn))
    return {"machine": platform.processor() or platform.machine(), "python": platform.python_version(),
            "spin_counts": list(spec.spin_counts), "systems": len(systems),
            "transitions_median": float(np.median(counts)), "transitions_max": int(max(counts)),
            "transition_list_ms": 1000 * per_system, "rendered_sample_ms": 1000 * per_sample,
            "record_points": points, "sampler_acceptance": sampler.acceptance_rate, "timing": timer.report()}
