"""Classify why a candidate list failed on a development sample.

Categories:
* `success`: some refined candidate matches the truth within tolerance;
* `structure_miss`: no candidate has the right components and group structure;
* `j_outside_basin`: the structure was proposed but refinement did not reach the truth;
* `unidentifiable`: a refined candidate differs from the truth in J yet fits
  held-out data as well as the truth does (observationally equivalent).
Only development data may be classified for feedback; frozen test failures are
never fed back into training.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence

from ..evaluation.matching import match_interpretations
from ..generator.sampler import Sample
from ..spinsystem import Interpretation

CATEGORIES = ("success", "structure_miss", "j_outside_basin", "unidentifiable")


@dataclass
class FailureRecord:
    sample: Sample
    category: str
    candidates: List[Interpretation] = field(default_factory=list)
    details: dict = field(default_factory=dict)


def classify_failure(sample: Sample, candidates: Sequence[Interpretation],
                     refined: Optional[Sequence] = None, tolerance_hz: float = 0.1,
                     truth_validation_residual: Optional[float] = None,
                     equivalence_ratio: float = 1.05) -> FailureRecord:
    truth = sample.interpretation
    pre = [match_interpretations(truth, c, tolerance_hz) for c in candidates]
    if refined:
        post = [match_interpretations(truth, r.interpretation, tolerance_hz) for r in refined]
        if any(m.j_match for m in post):
            return FailureRecord(sample, "success", list(candidates))
        if truth_validation_residual is not None:
            for r, m in zip(refined, post):
                if (m.structure and r.validation_residual is not None
                        and r.validation_residual <= equivalence_ratio * truth_validation_residual):
                    return FailureRecord(sample, "unidentifiable", list(candidates),
                                         {"validation_residual": r.validation_residual,
                                          "truth_validation_residual": truth_validation_residual,
                                          "max_j_error_hz": m.max_j_error_hz})
    elif any(m.j_match for m in pre):
        return FailureRecord(sample, "success", list(candidates))
    if not any(m.structure for m in pre):
        return FailureRecord(sample, "structure_miss", list(candidates))
    errors = [m.max_j_error_hz for m in pre if m.structure]
    return FailureRecord(sample, "j_outside_basin", list(candidates), {"best_pre_refine_j_error_hz": min(errors)})
