"""Candidate-quality metrics: coverage@k and permutation-matched J errors."""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Sequence

import numpy as np

from ..evaluation.matching import coupling_errors_by_category, match_interpretations
from ..spinsystem import Interpretation


class CoverageMetrics:
    def __init__(self, ks: Sequence[int] = (1, 3, 10), tolerance_hz: float = 1.0):
        self.ks = tuple(sorted(ks))
        self.tolerance_hz = tolerance_hz
        self.count = 0
        self.structure_hits = defaultdict(int)
        self.j_hits = defaultdict(int)
        self.component_count_hits = 0
        self.category_errors: Dict[str, List[float]] = defaultdict(list)
        self.empty = 0

    def update(self, truth: Interpretation, candidates: Sequence[Interpretation]) -> None:
        self.count += 1
        if not candidates:
            self.empty += 1
            return
        matches = [match_interpretations(truth, c, self.tolerance_hz) for c in candidates]
        for k in self.ks:
            top = matches[:k]
            self.structure_hits[k] += any(m.structure for m in top)
            self.j_hits[k] += any(m.j_match for m in top)
        self.component_count_hits += len(candidates[0]) == len(truth)
        best = next((i for i, m in enumerate(matches) if m.structure), None)
        if best is not None:
            for ti, ci in matches[best].pairs:
                errors = coupling_errors_by_category(truth.components[ti].system,
                                                     candidates[best].components[ci].system)
                for key, values in errors.items():
                    self.category_errors[key].extend(values)

    def summary(self) -> dict:
        n = max(1, self.count)
        out = {"count": self.count, "empty_candidate_lists": self.empty,
               "component_count_accuracy_top1": self.component_count_hits / n,
               "tolerance_hz": self.tolerance_hz}
        for k in self.ks:
            out[f"structure_coverage@{k}"] = self.structure_hits[k] / n
            out[f"j_coverage@{k}"] = self.j_hits[k] / n
        out["j_error_median_hz"] = {key: float(np.median(v)) for key, v in sorted(self.category_errors.items()) if v}
        out["j_error_p90_hz"] = {key: float(np.percentile(v, 90)) for key, v in sorted(self.category_errors.items()) if v}
        return out
