"""Rank candidate structures from an observed J network: route A and route B in both input modes, J and K (D56).

    python scripts/j_structure.py configs/j_networks/isopropylamine.json [--max-unseen 2] [--top 8] [--explain 1]
        [--model-dir runs/models] [--model FILE ...] [--routes A,J,K] [--out ranking.json]

Every run ranks with route A (rule likelihood, zulf_hypothesis.j_structure) and with route B (the learned
likelihood) once per input mode: J (coupling in Hz) and K (reduced coupling, isotope-independent). The route B
models are <model-dir>/j_edges_J.json and j_edges_K.json, written by scripts/j_structure_benchmark.py; --model
adds or replaces a model file (its saved mode decides the route). A missing model is reported, not trained here.

The JSON holds "units" ({label: {"h": protons, "copies": n}}), "protons" ({group: [unit, copy]}), "couplings"
({"J(a,b)": Hz}), optional "sigma", "isotopes" ({unit or group: isotope}, default 13C / 1H; route A assumes the
defaults, route B in mode K uses them) and "truth" (a bond list as printed, for checking).
"""
import argparse
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from zulf_hypothesis.j_structure import JObservation, describe, rank_structures, same_structure  # noqa: E402

DEFAULT_MODEL_DIR = Path(os.environ.get("ZULF_MODEL_WORKSPACE", ROOT / "runs")) / "models"


def route_likelihoods(model_dir, model_files=(), routes=("A", "J", "K")):
    """{route name: likelihood} in the order of `routes`: "A" (None, the rule likelihood) and "B-J" / "B-K" for every
    mode with a model (model_files override the files in model_dir); missing modes are listed in the second value."""
    from zulf_model.structure.edge_model import EdgeModel, LearnedLikelihood, load_route_b
    modes = [r for r in routes if r != "A"]
    found = load_route_b(model_dir, modes) if modes else {}
    for path in model_files:
        lik = LearnedLikelihood(EdgeModel.from_dict(json.loads(Path(path).read_text())))
        found[lik.mode] = lik
    out = {}
    for r in routes:
        if r == "A":
            out["A"] = None
        elif r in found:
            out[f"B-{r}"] = found[r]
    return out, [m for m in modes if m not in found]


def rank_all(obs, likelihoods, max_unseen=2, top=8):
    return {name: rank_structures(obs, max_unseen=max_unseen, top=top, likelihood=lik)
            for name, lik in likelihoods.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("observation", nargs="+")
    ap.add_argument("--max-unseen", type=int, default=2)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--explain", type=int, default=1, help="print the per-coupling terms of the best N")
    ap.add_argument("--model-dir", default=str(DEFAULT_MODEL_DIR), help="folder of j_edges_J.json / j_edges_K.json")
    ap.add_argument("--model", action="append", default=[], help="route B model file (repeatable; its mode decides)")
    ap.add_argument("--routes", default="A,J,K", help="A (rules), J and K (route B input modes)")
    ap.add_argument("--out", default=None, help="write every ranking as JSON")
    args = ap.parse_args()
    likelihoods, missing = route_likelihoods(args.model_dir, args.model, [r.strip() for r in args.routes.split(",")])
    if missing:
        print(f"no route B model for mode {', '.join(missing)} in {args.model_dir}; train with "
              f"python scripts/j_structure_benchmark.py --model-dir {args.model_dir}")
    report = {}
    for path in args.observation:
        data = json.load(open(path))
        obs = JObservation.from_dict(data)
        truth = data.get("truth")
        ranked = rank_all(obs, likelihoods, args.max_unseen, args.top)
        print(f"== {Path(path).stem}" + (f" (truth: {truth})" if truth else ""))
        if obs.labelled():
            print(f"   isotopes {obs.isotopes}: route A assumes 13C / 1H, route B-K is the isotope-independent one")
        entry = {}
        for name, cands in ranked.items():
            print(f"  -- route {name}")
            rows = []
            for k, c in enumerate(cands):
                hit = bool(truth and same_structure(obs, c, truth))
                print(f"  {k + 1:2d}. score {c.score:8.2f} (logL {c.log_likelihood:7.2f}, prior {c.log_prior:5.1f})  "
                      f"{describe(obs, c)}" + ("  <- truth" if hit else ""))
                rows.append({"structure": describe(obs, c), "score": c.score, "log_likelihood": c.log_likelihood,
                             "log_prior": c.log_prior, "truth": hit})
            for c in cands[:args.explain]:
                for e in c.explanation:
                    print(f"      {e['coupling']:12s} {e['J']:+9.3f} Hz  {e['bonds']} bonds  log p {e['log_p']:6.2f}")
            entry[name] = rows
        tops = {name: (describe(obs, cands[0]) if cands else None) for name, cands in ranked.items()}
        agree = len(set(tops.values())) == 1
        print("  routes agree on the best structure" if agree else
              "  routes differ: " + "; ".join(f"{n}: {t}" for n, t in tops.items()))
        report[Path(path).stem] = {"truth": truth, "routes": entry, "best": tops, "agree": agree}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
