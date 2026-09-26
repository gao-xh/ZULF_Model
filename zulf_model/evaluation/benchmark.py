"""Benchmark proposers + solver under a common protocol.

For each test sample the ground truth is rendered into a training observation
and independent held-out observations (fresh noise), every proposer returns k
candidates, each candidate is refined, and structure / J coverage are
recorded before and after refinement, together with evaluations and time.
"""
from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..generator.sampler import Sample
from zulf_core.physics.protocol import SUDDEN_DROP, Protocol
from zulf_core.render.acquisition import Acquisition
from ..render.perturb import PerturbationConfig, mixture_fid, sample_render_params
from zulf_core.render.renderer import Renderer
from zulf_core.physics.transitions import TransitionCache
from zulf_core.solver.observed import ObservedSpectrum
from zulf_core.solver.refine import RefineSettings, refine_candidates
from zulf_core.timing import Timer
from zulf_core.evaluation.matching import match_interpretations
from .proposers import CandidateProposer


@dataclass(frozen=True)
class BenchmarkConfig:
    k: int = 5
    held_out: int = 2
    ranges: tuple = ((1.0, 400.0),)
    tolerance_hz: float = 0.1
    refine: bool = True
    seed: int = 0


def observe_sample(sample: Sample, acquisition: Acquisition, perturbation: PerturbationConfig, ranges,
                   rng: np.random.Generator, held_out: int, protocol: Protocol = SUDDEN_DROP,
                   cache: Optional[TransitionCache] = None):
    """One training and several held-out observations sharing signal parameters, independent noise."""
    cache = cache or TransitionCache()
    renderer = Renderer(acquisition)
    transitions = [cache.get(c.system, protocol) for c in sample.interpretation.components]
    contributions = [c.contribution for c in sample.interpretation.components]
    params = sample_render_params(rng, transitions, perturbation)
    clean = mixture_fid(renderer, transitions, params, contributions)
    peak = float(np.abs(np.fft.rfft(clean)).max() / acquisition.n) or 1.0
    observations = []
    for k in range(1 + held_out):
        noise_std = peak / (params.snr or np.inf) * np.sqrt(acquisition.n)
        noisy = clean + np.random.default_rng([params.noise_seed, k]).normal(0, noise_std, acquisition.points)
        observations.append(ObservedSpectrum.from_fid(noisy, acquisition, ranges, label=f"obs{k}"))
    return observations[0], observations[1:], params


def run_benchmark(samples: Sequence[Sample], proposers: Sequence[CandidateProposer], acquisition: Acquisition,
                  perturbation: PerturbationConfig, config: BenchmarkConfig = BenchmarkConfig(),
                  settings: RefineSettings = RefineSettings(starts=1), output: Optional[Path] = None) -> dict:
    timer = Timer()
    rng = np.random.default_rng(config.seed)
    cache = TransitionCache()
    rows: List[dict] = []
    for index, sample in enumerate(samples):
        train, held, params = observe_sample(sample, acquisition, perturbation, config.ranges, rng,
                                             config.held_out, cache=cache)
        for proposer in proposers:
            start = time.perf_counter()
            with timer.section(f"propose.{proposer.name}"):
                candidates = proposer.propose(train, config.k)
            propose_s = time.perf_counter() - start
            pre = [match_interpretations(sample.interpretation, c, config.tolerance_hz) for c in candidates]
            row = {"sample": index, "family_id": sample.family_id, "proposer": proposer.name,
                   "candidates": len(candidates), "propose_s": propose_s,
                   "structure_hit_pre": any(m.structure for m in pre),
                   "j_hit_pre": any(m.j_match for m in pre)}
            if config.refine and candidates:
                start = time.perf_counter()
                with timer.section(f"refine.{proposer.name}"):
                    results = refine_candidates(candidates, train, held, settings, timer=timer)
                post = [match_interpretations(sample.interpretation, r.interpretation, config.tolerance_hz)
                        for r in results]
                row.update(refine_s=time.perf_counter() - start,
                           evaluations=int(sum(r.evaluations for r in results)),
                           j_hit_post=any(m.j_match for m in post),
                           top1_j_hit_post=bool(post and post[0].j_match),
                           best_validation_residual=results[0].validation_residual if results else None)
            rows.append(row)
    summary: Dict[str, dict] = {}
    for proposer in proposers:
        mine = [r for r in rows if r["proposer"] == proposer.name]
        n = max(1, len(mine))
        summary[proposer.name] = {key: sum(bool(r.get(key)) for r in mine) / n
                                  for key in ("structure_hit_pre", "j_hit_pre", "j_hit_post", "top1_j_hit_post")}
        summary[proposer.name]["mean_evaluations"] = float(np.mean([r.get("evaluations", 0) for r in mine]))
        summary[proposer.name]["mean_seconds"] = float(np.mean([r["propose_s"] + r.get("refine_s", 0) for r in mine]))
    report = {"config": asdict(config), "summary": summary, "rows": rows, "timing": timer.report()}
    if output is not None:
        Path(output).write_text(json.dumps(report, indent=2, default=float), encoding="utf-8")
    return report
