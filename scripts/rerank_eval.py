"""Evaluate candidate decoding and solver reranking on simulated observations.

For validation samples of a run configuration, render a noisy FID (random phase
and delay, as in training), phase it automatically, propose k candidates with a
trained checkpoint under each J decoding, then refine every candidate with a
short solver run and rerank by residual. Reports structure and J coverage
before and after reranking.

Usage:
    python scripts/rerank_eval.py RUN_CONFIG MODEL.pt OUTPUT.json [--samples 30] [--k 10]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from zulf_model.evaluation import ModelProposer
from zulf_model.evaluation.benchmark import observe_sample
from zulf_core.evaluation.matching import match_interpretations
from zulf_model.models import load_model
from zulf_model.render import SpectrumGrid
from zulf_core.solver import ParameterPolicy, RefineSettings, refine
from zulf_model.training import TrainingSetup


def coverage(truth, candidates, tolerance_hz):
    matches = [match_interpretations(truth, c, tolerance_hz) for c in candidates]
    return {"structure@1": bool(matches and matches[0].structure), "structure@k": any(m.structure for m in matches),
            "j@1": bool(matches and matches[0].j_match), "j@k": any(m.j_match for m in matches)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_config")
    parser.add_argument("model")
    parser.add_argument("output")
    parser.add_argument("--samples", type=int, default=30)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--tolerance-hz", type=float, default=1.0)
    parser.add_argument("--refine-seconds", type=float, default=10.0)
    parser.add_argument("--no-rerank", action="store_true")
    args = parser.parse_args()
    setup = TrainingSetup(args.run_config)
    model = load_model(args.model).eval()
    acquisition = setup.processing.base()
    grid = SpectrumGrid.from_spec(model.spec.grid, acquisition)
    ranges = [(model.spec.grid.f_min_hz, model.spec.grid.f_max_hz)]
    rng = np.random.default_rng(args.seed)
    settings = RefineSettings(starts=1, continuation_rates_per_s=(3.0, 1.0, 0.0), max_nfev=80,
                              max_seconds=args.refine_seconds,
                              policy=ParameterPolicy(coupling_margin_hz=3.0, rate_bounds_per_s=(0.05, 20.0)))
    rows = []
    for index, sample in enumerate(setup.sampler.generate(args.samples, seed=args.seed, split="val")):
        observed, _, params = observe_sample(sample, acquisition, setup.perturbation, ranges, rng, held_out=0)
        truth = sample.interpretation
        row = {"sample": index, "components": len(truth.components),
               "spins": [c.system.n_spins for c in truth.components]}
        decoded = {}
        for decoding in ("offset", "expectation"):
            proposer = ModelProposer(model, grid, model.spec.grid.channels, phasing={"auto": True},
                                     j_decoding=decoding)
            candidates = proposer.propose(observed, args.k)
            decoded[decoding] = candidates
            row[decoding] = coverage(truth, candidates, args.tolerance_hz)
        if not args.no_rerank and decoded["expectation"]:
            start = time.perf_counter()
            scored = []
            for candidate in decoded["expectation"]:
                try:
                    result = refine(candidate, observed, settings)
                except ValueError:
                    continue
                scored.append((result.relative_residual, result.interpretation))
            scored.sort(key=lambda item: item[0])
            reranked = [interp for _, interp in scored]
            row["reranked"] = coverage(truth, reranked, args.tolerance_hz)
            row["reranked_tight"] = coverage(truth, reranked, 0.1)
            row["rerank_seconds"] = time.perf_counter() - start
        rows.append(row)
        print(json.dumps(row), flush=True)
    summary = {}
    for key in ("offset", "expectation", "reranked", "reranked_tight"):
        mine = [r[key] for r in rows if key in r]
        if mine:
            summary[key] = {metric: float(np.mean([m[metric] for m in mine])) for metric in mine[0]}
    summary["samples"] = len(rows)
    Path(args.output).write_text(json.dumps({"summary": summary, "rows": rows, "args": vars(args)}, indent=2),
                                 encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
