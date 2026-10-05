"""Rank candidate structures from an observed J network (route A of J -> structure; zulf_hypothesis.j_structure).

    python scripts/j_structure.py configs/j_networks/isopropylamine.json [--max-unseen 2] [--top 8] [--explain 1]
        [--model runs/models/j_edges.json]     (route B: the learned likelihood of scripts/j_structure_benchmark.py)

The JSON holds "units" ({label: {"h": protons, "copies": n}}), "protons" ({group: [unit, copy]}), "couplings"
({"J(a,b)": Hz}), optional "sigma" and "truth" (a bond list as printed, for checking).
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from zulf_hypothesis.j_structure import JObservation, describe, rank_structures, same_structure  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("observation", nargs="+")
    ap.add_argument("--max-unseen", type=int, default=2)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--explain", type=int, default=1, help="print the per-coupling terms of the best N")
    ap.add_argument("--model", default=None, help="route B model JSON (default: route A rule likelihood)")
    args = ap.parse_args()
    likelihood = None
    if args.model:
        from zulf_model.structure.edge_model import EdgeModel, LearnedLikelihood
        likelihood = LearnedLikelihood(EdgeModel.from_dict(json.load(open(args.model))))
    for path in args.observation:
        data = json.load(open(path))
        obs = JObservation.from_dict(data)
        ranked = rank_structures(obs, max_unseen=args.max_unseen, top=args.top, likelihood=likelihood)
        truth = data.get("truth")
        print(f"== {Path(path).stem}" + (f" (truth: {truth})" if truth else ""))
        for k, c in enumerate(ranked):
            mark = "  <- truth" if truth and same_structure(obs, c, truth) else ""
            print(f"  {k + 1:2d}. score {c.score:8.2f} (logL {c.log_likelihood:7.2f}, prior {c.log_prior:5.1f})  "
                  f"{describe(obs, c)}{mark}")
        for c in ranked[:args.explain]:
            for e in c.explanation:
                print(f"      {e['coupling']:12s} {e['J']:+9.3f} Hz  {e['bonds']} bonds  log p {e['log_p']:6.2f}")


if __name__ == "__main__":
    main()
