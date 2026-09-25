"""Solver recovery on simulated systems (plan Phase 3).

For systems drawn from the generator, render a noisy observation, perturb
every observable coupling of the truth by a normal draw of a given scale,
refine, and count recoveries (max coupling error after relabelling within the
tolerance). Strategies: local refinement with matched linewidth continuation,
and global pattern search followed by refinement.

Usage:
    python scripts/solver_recovery.py OUTPUT.json [--systems 8] [--trials 3] [--scales 0.5,2,5]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from zulf_model.evaluation.benchmark import observe_sample
from zulf_model.evaluation.identifiability import basin_of_attraction
from zulf_model.generator import build_default_sampler
from zulf_model.render import Acquisition, PerturbationConfig
from zulf_model.solver import ParameterPolicy, RefineSettings
from zulf_model.spec import ProblemSpec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--systems", type=int, default=8)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--scales", default="0.5,2,5")
    parser.add_argument("--spin-counts", default="4,5,6")
    parser.add_argument("--tolerance-hz", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--strategies", default="local,search,anneal")
    parser.add_argument("--budget-s", type=float, default=900.0,
                        help="Wall-clock budget per refinement; budget exhaustion is recorded in the flags.")
    args = parser.parse_args()
    spec = ProblemSpec(spin_counts=tuple(int(v) for v in args.spin_counts.split(",")), max_components=1)
    sampler = build_default_sampler(spec, {"mixture": {"extra_molecule_probability": 0.0}})
    acquisition = Acquisition(1000.0, 8192, start_sample=50, sg_window=201, sg_order=2)
    perturbation = PerturbationConfig(snr_range=(100.0, 300.0), gaussian_probability=0.0, drift_probability=0.0,
                                      interference_probability=0.0, family_split_probability=0.0,
                                      rate_range_per_s=(0.5, 2.0))
    scales = [float(v) for v in args.scales.split(",")]
    policy = ParameterPolicy(coupling_margin_hz=5.0, rate_bounds_per_s=(0.1, 10.0))
    strategies = {
        "local": RefineSettings(starts=1, policy=policy, max_seconds=args.budget_s, max_evaluations=20000),
        "search": RefineSettings(policy=policy, max_seconds=2 * args.budget_s, max_evaluations=20000,
                                 search=dict(popsize=10, maxiter=40, max_seconds=90, solutions=3)),
        "anneal": RefineSettings(policy=policy, max_seconds=2 * args.budget_s, max_evaluations=20000,
                                 search=dict(method="dual_annealing", annealing_maxiter=200, max_seconds=90,
                                             solutions=3)),
    }
    chosen = [s for s in args.strategies.split(",") if s in strategies]
    rng = np.random.default_rng(args.seed)
    rows = []
    for index, sample in enumerate(sampler.generate(args.systems, seed=args.seed, split="val")):
        if len(sample.interpretation.components) != 1:
            continue
        ranges = [(1.0, 400.0)]
        observed, _, _ = observe_sample(sample, acquisition, perturbation, ranges, rng, held_out=0)
        system = sample.interpretation.components[0].system
        for name in chosen:
            start = time.perf_counter()
            result = basin_of_attraction(sample.interpretation, observed, scales, trials=args.trials,
                                         tolerance_hz=args.tolerance_hz, settings=strategies[name], seed=index)
            for row in result:
                rows.append({"system": index, "isotopes": list(system.isotopes),
                             "groups": [len(g) for g in system.groups], "strategy": name,
                             "seconds": time.perf_counter() - start, **row})
            print(json.dumps(rows[-len(result):]), flush=True)
    summary = {}
    for name in chosen:
        for scale in scales:
            mine = [r for r in rows if r["strategy"] == name and r["scale_hz"] == scale]
            if mine:
                summary[f"{name}@{scale}"] = {"success_rate": float(np.mean([r["success_rate"] for r in mine])),
                                              "median_evaluations": float(np.median([r["median_evaluations"]
                                                                                     for r in mine])),
                                              "systems": len(mine)}
    Path(args.output).write_text(json.dumps({"summary": summary, "rows": rows, "args": vars(args)}, indent=2),
                                 encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
