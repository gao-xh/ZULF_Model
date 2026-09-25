"""Evaluate several checkpoints on several validation distributions (same samples for every model).

Usage:
    python scripts/cross_eval.py OUTPUT.json --models name=path ... --data name=run_config ... [--count 512]

Each data entry renders `count` validation samples (split "val", fixed seed)
with its run configuration; every model proposes k candidates for every sample
and coverage metrics are computed with the same matching rules as training.
"""
import argparse
import json
import time
from pathlib import Path

import torch

from zulf_model.models import load_model
from zulf_model.training import CoverageMetrics, TrainingSetup
from zulf_model.training.data import FixedDataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--models", nargs="+", required=True)
    parser.add_argument("--data", nargs="+", required=True)
    parser.add_argument("--count", type=int, default=512)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--batch", type=int, default=32)
    parser.add_argument("--tolerance-hz", type=float, default=1.0)
    args = parser.parse_args()
    models = {name: load_model(path).eval() for name, path in (m.split("=", 1) for m in args.models)}
    report = {"args": vars(args), "results": {}}
    for data_name, config in (d.split("=", 1) for d in args.data):
        setup = TrainingSetup(config)
        samples = list(setup.sampler.generate(args.count, seed=args.seed, split="val"))
        dataset = FixedDataset(samples, setup.renderer, setup.codec, seed=args.seed)
        collate = setup.collator()
        batches = [collate([dataset[i] for i in range(s, min(s + args.batch, len(dataset)))])
                   for s in range(0, len(dataset), args.batch)]
        for model_name, model in models.items():
            start = time.perf_counter()
            metrics = CoverageMetrics(ks=(1, 3, args.k), tolerance_hz=args.tolerance_hz)
            with torch.no_grad():
                for batch in batches:
                    proposals = model.propose(batch["features"], batch["frequency_hz"], k=args.k)
                    for truth, candidates in zip(batch["interpretations"], proposals):
                        metrics.update(truth, candidates)
            summary = metrics.summary()
            summary["seconds"] = time.perf_counter() - start
            report["results"][f"{model_name}@{data_name}"] = summary
            print(model_name, data_name, json.dumps({k: round(v, 3) for k, v in summary.items()
                                                     if isinstance(v, float)}), flush=True)
    Path(args.output).write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
