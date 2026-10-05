"""Benchmark J -> structure routes A (rule likelihood) and B (learned bond-count likelihood).

    python scripts/j_structure_benchmark.py [--samples 100] [--heavy-atoms 2,5] [--train-graphs 4000]
        [--model runs/models/j_edges.json] [--out runs/processed/j_structure_benchmark.json]

Route B is trained on generator pairs (seed --train-seed, configs/couplings_v1.json) unless --model names an
existing model file (it is written there after training otherwise). The synthetic test set uses a different seed
and the same coupling rules, so it favours route B; the real J networks (configs/j_networks/*.json) are the fair
comparison. Prints top-1/top-3 rates and the rank of the truth per real network.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_hypothesis.j_structure import JObservation, rank_structures, same_structure  # noqa: E402
from zulf_model.generator.couplings import CouplingRules  # noqa: E402
from zulf_model.generator.graphs import GraphConfig, random_graph  # noqa: E402
from zulf_model.structure.edge_model import (EdgeModel, LearnedLikelihood, accuracy, generate_rows,  # noqa: E402
                                             train_edge_model)
from zulf_model.structure.observations import observation_from_graph  # noqa: E402


def truth_rank(obs, truth, likelihood, top):
    ranked = rank_structures(obs, top=top, likelihood=likelihood)
    k = next((i for i, c in enumerate(ranked) if same_structure(obs, c, truth)), None)
    margin = ranked[0].score - ranked[1].score if k == 0 and len(ranked) > 1 else None
    return k, margin


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--heavy-atoms", default="2,5")
    ap.add_argument("--max-atoms", type=int, default=5, help="skip observations with more carbon copies")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--train-graphs", type=int, default=4000)
    ap.add_argument("--train-seed", type=int, default=1)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--model", default=None)
    ap.add_argument("--rules", default=str(ROOT / "configs" / "couplings_v1.json"))
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    rules = CouplingRules.load(args.rules)
    lo, hi = (int(v) for v in args.heavy_atoms.split(","))
    t0 = time.time()
    if args.model and Path(args.model).exists():
        model = EdgeModel.from_dict(json.loads(Path(args.model).read_text()))
        print(f"route B model: {args.model}")
    else:
        rows, hyb = generate_rows(args.train_graphs, seed=args.train_seed, rules=rules, with_hybrid=True)
        model = train_edge_model(rows, epochs=args.epochs, seed=args.train_seed, hybrid_rows=hyb)
        held = generate_rows(500, seed=args.train_seed + 1000, rules=rules)
        print(f"route B trained on {len(rows)} couplings from {args.train_graphs} graphs in {time.time() - t0:.0f} s;"
              f" held-out bond-count accuracy {accuracy(model, held):.3f}")
        if args.model:
            Path(args.model).parent.mkdir(parents=True, exist_ok=True)
            Path(args.model).write_text(json.dumps(model.to_dict()))
    routes = {"A": None, "B": LearnedLikelihood(model)}
    rng = np.random.default_rng(args.seed)
    cfg = GraphConfig(heavy_atoms=(lo, hi))
    ranks = {r: [] for r in routes}
    n = 0
    t1 = time.time()
    while n < args.samples:
        out = observation_from_graph(rng, random_graph(rng, cfg), rules)
        if out is None:
            continue
        obs = JObservation.from_dict(out[0])
        if len(obs.atoms()) > args.max_atoms:
            continue
        n += 1
        for r, lik in routes.items():
            ranks[r].append(truth_rank(obs, out[1], lik, args.top)[0])
    report = {"synthetic": {}, "real": {}, "samples": n, "heavy_atoms": [lo, hi]}
    for r, ks in ranks.items():
        top1 = float(np.mean([k == 0 for k in ks]))
        top3 = float(np.mean([k is not None and k < 3 for k in ks]))
        report["synthetic"][r] = {"top1": top1, "top3": top3}
        print(f"synthetic ({n}, heavy atoms {lo}-{hi}) route {r}: top1 {top1:.2f}  top3 {top3:.2f}")
    print(f"  ({time.time() - t1:.0f} s)")
    for path in sorted((ROOT / "configs" / "j_networks").glob("*.json")):
        data = json.loads(path.read_text())
        obs = JObservation.from_dict(data)
        row = {}
        for r, lik in routes.items():
            k, margin = truth_rank(obs, data["truth"], lik, args.top)
            row[r] = {"rank": None if k is None else k + 1, "margin": margin}
        report["real"][path.stem] = row
        print(f"{path.stem:20s} " + "  ".join(
            f"{r}: rank {v['rank']}" + (f" (margin {v['margin']:.2f})" if v["margin"] is not None else "")
            for r, v in row.items()))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
