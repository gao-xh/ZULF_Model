"""Route B (learned J -> structure likelihood) trained on J versus on the reduced coupling K (D51, D52).

    python scripts/j_k_learning_comparison.py [--seeds 1,2,3] [--train-graphs 4000] [--epochs 20]
        [--held-graphs 500] [--samples 100] [--out runs/processed/j_vs_k]

For every training seed both models are trained on the same generator rows (configs/couplings_v1.json, 13C / 1H)
and tested on the same data:
- bond-count accuracy on held-out protonated rows (seed + 1000);
- bond-count accuracy on the same held-out rows fully deuterated: J(13C,2H) = J(13C,1H) g_2H / g_1H and
  J(2H,2H) = J(1H,1H) (g_2H / g_1H)^2, the element pair (kind) unchanged; a transfer test, nothing deuterated is
  trained (the primary isotope effect on K is neglected, so this isolates the gamma scaling);
- rank of the true structure on synthetic observations (top-1, top-3; seed --test-seed) and on the real amine
  J networks (configs/j_networks).
Mean and sample standard deviation over the seeds. Writes OUT/summary.json and prints a table.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from j_structure_benchmark import truth_rank  # noqa: E402
from zulf_core.nuclei import get_registry  # noqa: E402
from zulf_hypothesis.j_structure import JObservation  # noqa: E402
from zulf_model.generator.couplings import CouplingRules  # noqa: E402
from zulf_model.generator.graphs import GraphConfig, random_graph  # noqa: E402
from zulf_model.structure.edge_model import LearnedLikelihood, accuracy, generate_rows, train_edge_model  # noqa: E402
from zulf_model.structure.observations import observation_from_graph  # noqa: E402

DEUTERATED = {"CH": ("13C", "2H"), "HH": ("2H", "2H")}


def deuterate(rows):
    r = get_registry().gamma("2H") / get_registry().gamma("1H")
    return [(j * (r if kind == "CH" else r * r), kind, h, b) for j, kind, h, b in rows]


def synthetic_cases(rules, n, seed, heavy=(2, 5), max_atoms=5):
    rng = np.random.default_rng(seed)
    cfg = GraphConfig(heavy_atoms=heavy)
    cases = []
    while len(cases) < n:
        out = observation_from_graph(rng, random_graph(rng, cfg), rules)
        if out is None:
            continue
        obs = JObservation.from_dict(out[0])
        if len(obs.atoms()) <= max_atoms:
            cases.append((obs, out[1]))
    return cases


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", default="1,2,3")
    ap.add_argument("--train-graphs", type=int, default=4000)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--held-graphs", type=int, default=500)
    ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--test-seed", type=int, default=7)
    ap.add_argument("--rules", default=str(ROOT / "configs" / "couplings_v1.json"))
    ap.add_argument("--out", default="runs/processed/j_vs_k")
    args = ap.parse_args()
    rules = CouplingRules.load(args.rules)
    seeds = [int(v) for v in args.seeds.split(",")]
    cases = synthetic_cases(rules, args.samples, args.test_seed)
    real = [(p.stem, json.loads(p.read_text())) for p in sorted((ROOT / "configs" / "j_networks").glob("*.json"))]
    runs = []
    for seed in seeds:
        t0 = time.time()
        rows, hyb = generate_rows(args.train_graphs, seed=seed, rules=rules, with_hybrid=True)
        held = generate_rows(args.held_graphs, seed=seed + 1000, rules=rules)
        held_d = deuterate(held)
        for mode in ("J", "K"):
            model = train_edge_model(rows, epochs=args.epochs, seed=seed, hybrid_rows=hyb, mode=mode)
            lik = LearnedLikelihood(model)
            ranks = [truth_rank(obs, truth, lik, 10)[0] for obs, truth in cases]
            row = {"seed": seed, "mode": mode, "train_rows": len(rows),
                   "accuracy_1H": accuracy(model, held), "accuracy_2H": accuracy(model, held_d, DEUTERATED),
                   "top1": float(np.mean([k == 0 for k in ranks])),
                   "top3": float(np.mean([k is not None and k < 3 for k in ranks])),
                   "real_ranks": {}}
            for name, data in real:
                k = truth_rank(JObservation.from_dict(data), data["truth"], lik, 10)[0]
                row["real_ranks"][name] = None if k is None else k + 1
            runs.append(row)
            print(f"seed {seed} {mode}: acc 1H {row['accuracy_1H']:.3f}  acc 2H {row['accuracy_2H']:.3f}  "
                  f"top1 {row['top1']:.2f}  top3 {row['top3']:.2f}  real {row['real_ranks']}  "
                  f"({time.time() - t0:.0f} s)", flush=True)
    summary = {"args": vars(args), "runs": runs, "mean": {}, "std": {}}
    for mode in ("J", "K"):
        sel = [r for r in runs if r["mode"] == mode]
        for key in ("accuracy_1H", "accuracy_2H", "top1", "top3"):
            vals = np.array([r[key] for r in sel])
            summary["mean"].setdefault(mode, {})[key] = float(vals.mean())
            summary["std"].setdefault(mode, {})[key] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
    print("\nmean +- std over seeds", seeds)
    for mode in ("J", "K"):
        m, s = summary["mean"][mode], summary["std"][mode]
        print(f"  {mode}: " + "  ".join(f"{k} {m[k]:.3f} +- {s[k]:.3f}" for k in m))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(out / "summary.json")


if __name__ == "__main__":
    main()
