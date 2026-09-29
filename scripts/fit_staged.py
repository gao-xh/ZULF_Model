"""Staged known-structure fit of a processed real spectrum: some couplings held first, then everything free.

    python scripts/fit_staged.py --freq F.npy --values V.npy --id NAME \\
        --structure '{"motif": "pyridine ring", "one_bond": {...}}' --couplings '{"J(HA2,HA3)": 4.9, ...}' \\
        [--hold-first J(HA] [--starts 12] [--spread 1.0] [--variant ratios] [--range 140,200] [--out DIR]

Stage 1 holds every coupling whose key starts with one of the `--hold-first` prefixes at its starting value (e.g.
the H-H couplings at accepted values) and fits the rest; stage 2 frees all couplings from the stage-1 optimum, plus
`--starts - 1` starts perturbed around it (normal, `--spread` Hz on couplings). On pyridine (x = 1.00) this found a
lower minimum than a direct all-free fit from the same accepted values (ANALYSIS_LOG). Outputs OUT/fit.json (stage
results, every stage-2 start, linearised uncertainties), OUT/uncertainty.md and OUT/fit.png (data, model,
residual). The result is a conditional numerical fit: check the start agreement and the residual.
"""
import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from zulf_core.solver import ObservedSpectrum, refine                  # noqa: E402
from zulf_hypothesis.builder import build_model                        # noqa: E402
from zulf_hypothesis.fit import default_fit_base, protonated           # noqa: E402
from zulf_hypothesis.scoring import free_parameter_count               # noqa: E402
from zulf_hypothesis.search import Evaluated, _settings_for            # noqa: E402
from zulf_hypothesis.uncertainty import coupling_uncertainties, start_agreement, uncertainty_markdown  # noqa: E402


def figure(obs, prediction, title, path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    f, y, p = obs.frequencies_hz, obs.values.real, np.asarray(prediction).real
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(f, y, color="#222222", lw=0.8, label="experiment")
    ax.plot(f, p, color="#d1495b", lw=0.9, label="simulation (fit)")
    ax.plot(f, y - p - 0.25 * np.abs(y).max(), color="#888888", lw=0.7, label="residual (offset)")
    ax.set_xlabel("frequency (Hz)")
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required=True)
    ap.add_argument("--values", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--protonate", action="store_true")
    ap.add_argument("--hold-first", default="J(HA", help="comma-separated key prefixes held in stage 1")
    ap.add_argument("--starts", type=int, default=12, help="stage-2 starts (the stage-1 optimum plus perturbed)")
    ap.add_argument("--spread", type=float, default=1.0)
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    f, v = np.load(args.freq).astype(float), np.load(args.values).astype(float)
    lo, hi = (float(x) for x in args.range.split(","))
    obs = ObservedSpectrum.from_spectrum(f, v, [(lo, hi)], record=None, real_only=True, label=args.id)
    spec = json.loads(args.structure)
    spec.setdefault("compound", args.id)
    fragment = reg.structure_for(spec)
    if args.protonate:
        fragment = protonated(fragment)
    fragment = override_couplings(fragment, json.loads(args.couplings))
    model = build_model(fragment, ranges=[(lo, hi)])
    base = dataclasses.replace(default_fit_base(), starts=1, max_evaluations=10 ** 6,
                               max_seconds=600.0 * max(args.starts, 1))
    settings = _settings_for(model, args.variant, base)
    out = Path(args.out or f"runs/processed/{args.id}_staged")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    prefixes = tuple(p for p in args.hold_first.split(",") if p)
    held = {n for key, names in model.coupling_names.items() if key.startswith(prefixes) for n in names}
    p1 = settings.parameterize(model.interpretation)
    p1.fix(*[n for n in p1.order if n in held])
    r1 = refine(model.interpretation, obs, settings, parameterization=p1)
    p2 = settings.parameterize(model.interpretation)
    for n, value in r1.parameters.items():
        if n in p2.parameters:
            p2.set(n, value)
    x0 = p2.vector()
    lower, upper = p2.bounds()
    coupling = np.array([p2.parameters[n].kind == "coupling" for n in p2.free_names])
    rng = np.random.default_rng(args.seed)
    starts = [x0] + [np.clip(x0 + np.where(coupling, rng.normal(0, args.spread, len(x0)), 0.0), lower, upper)
                     for _ in range(max(args.starts - 1, 0))]
    r2 = refine(model.interpretation, obs, settings, parameterization=p2, initial_points=starts)
    e = Evaluated(model.name, args.variant, "staged", model, r2.summary(), r2.prediction,
                  k=free_parameter_count(model, settings))
    e.couplings = model.named_couplings(r2.parameters)
    row = {"id": args.id, "model": model.name, "variant": args.variant, "seconds": round(time.time() - t0),
           "hold_first": list(prefixes),
           "stage1": {"data_region_residual": r1.data_region_residual,
                      "couplings": {k: round(x, 3) for k, x in model.named_couplings(r1.parameters).items()}},
           "stage2": {"data_region_residual": r2.data_region_residual, "relative_residual": r2.relative_residual,
                      "flags": r2.flags, "boundary_hits": r2.boundary_hits,
                      "couplings": {k: round(x, 3) for k, x in e.couplings.items()},
                      "rates_per_s": {n: float(np.exp(x)) for n, x in r2.parameters.items() if "log_rate" in n},
                      "delay_ms": 1e3 * r2.parameters.get("phase_delay", 0.0)},
           "k": e.k}
    try:
        u = coupling_uncertainties(e, obs, base)
        row["uncertainty"] = u
        agreement = start_agreement(e, base)
        row["start_agreement"] = agreement
        (out / "uncertainty.md").write_text(
            uncertainty_markdown(u) + f"\n\nStage-2 starts: {agreement.get('starts')}, within 1 % of the best "
            f"score: {agreement.get('near_best')}\n")
    except Exception as exc:
        row["uncertainty_error"] = f"{type(exc).__name__}: {exc}"
    json.dump(row, open(out / "fit.json", "w"), indent=1, default=str)
    np.save(out / "prediction.npy", np.asarray(r2.prediction))
    figure(obs, r2.prediction, f"{args.id}: {model.name} [{args.variant}], staged fit; relative residual on the "
           f"data cores {r2.data_region_residual:.3f}, k {e.k}", out / "fit.png")
    print(json.dumps({"id": args.id, "residual": round(r2.data_region_residual, 4), "k": e.k,
                      "seconds": row["seconds"], "near_best": row.get("start_agreement", {}).get("near_best"),
                      "starts": row.get("start_agreement", {}).get("starts")}), flush=True)


if __name__ == "__main__":
    main()
