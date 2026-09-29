"""Staged known-structure fit of a processed real spectrum: some couplings held first, then everything free.

    python scripts/fit_staged.py --freq F.npy --values V.npy --id NAME \\
        --structure '{"motif": "pyridine ring", "one_bond": {...}}' --couplings '{"J(HA2,HA3)": 4.9, ...}' \\
        [--hold-first J(HA] [--starts 12] [--spread 1.0] [--variant ratios] [--range 140,200] [--out DIR]

Chemical exchange (`--exchange HA1`): the listed proton groups (labels, comma-separated) exchange with the
solvent at a fitted rate k (1/s, log-parameterized, one value tied across the isotopologues; start
`--kex-start`, bounds `--kex-bounds`); the model is then a Liouville-space one (zulf_core.physics.exchange) and
its derivatives are numerical. `--device numpy|cpu|cuda` selects the dense linear-algebra backend of the
exchange model (cpu/cuda: PyTorch; cuda needs a GPU); the environment variable ZULF_LINALG_DEVICE does the same.

Optional Gaussian priors (`--prior-sigma-hh`, `--prior-sigma-ch`, `--prior-weight`) keep the small couplings
(|J| < 50 Hz) near their starting values (e.g. literature) with a weight, while every coupling stays free
(RefineSettings.priors). Stage 1 holds every coupling whose key starts with one of the `--hold-first` prefixes at its starting value (e.g.
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


def group_index(model, label: str, component: int) -> int:
    """Group index of a proton-group label in one component, from the coupling parameter names."""
    sets = []
    for key, names in model.coupling_names.items():
        members = key[2:-1].split(",")
        if label not in members:
            continue
        for n in names:
            if n.startswith(f"c{component}.J"):
                a, b = n.split(".J")[1].split("-")
                sets.append({int(a), int(b)})
    common = set.intersection(*sets) if len(sets) >= 2 else set()
    if len(common) != 1:
        raise ValueError(f"Cannot locate group {label} in component {component}.")
    return common.pop()


def add_exchange(param, model, labels, start, bounds, fixed=False):
    """One exchange-rate parameter per label, tied across the components that contain the label."""
    names = []
    for label in labels:
        leader = None
        for c in range(len(model.component_labels)):
            try:
                g = group_index(model, label, c)
            except ValueError:
                continue
            name = param.add_exchange(c, g, start, bounds, free=not fixed, name=f"c{c}.log_kex_{label}")
            if leader is None:
                leader = name
                names.append(name)
            else:
                param.tie(leader, name)
    return names


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
    ap.add_argument("--prior-sigma-hh", type=float, default=0.0,
                    help="Gaussian prior width (Hz) of every H-H coupling around its starting value (0: none)")
    ap.add_argument("--prior-sigma-ch", type=float, default=0.0,
                    help="prior width (Hz) of every other small coupling (|J| < 50 Hz) around its starting value")
    ap.add_argument("--prior-weight", type=float, default=1.0,
                    help="prior weight: 1 = a one-sigma deviation costs one data point one noise sigma off")
    ap.add_argument("--exchange", default="", help="proton group labels exchanging with the solvent, e.g. HA1")
    ap.add_argument("--kex-start", type=float, default=10.0, help="starting exchange rate (1/s)")
    ap.add_argument("--kex-bounds", default="0.01,1e5", help="exchange rate bounds (1/s)")
    ap.add_argument("--kex-fixed", action="store_true", help="hold the exchange rate at --kex-start")
    ap.add_argument("--delay-bounds", default="", help="lo,hi in s for the fitted delay (default -0.01,0.01; write --delay-bounds=-0.03,0.03)")
    ap.add_argument("--device", default="", help="exchange linear algebra: numpy (default), cpu or cuda (torch)")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out or f"runs/processed/{args.id}_staged"), "fit_staged")
    if args.device:
        from zulf_core.physics import exchange as exchange_backend
        exchange_backend.set_backend("numpy") if args.device == "numpy" else \
            exchange_backend.set_backend("torch", args.device)
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
    delay_bounds = tuple(float(v) for v in args.delay_bounds.split(",")) if args.delay_bounds else None
    base = dataclasses.replace(default_fit_base(delay_bounds), starts=1, max_evaluations=10 ** 6,
                               max_seconds=600.0 * max(args.starts, 1))
    settings = _settings_for(model, args.variant, base)
    # priors around the starting couplings (e.g. literature values): all couplings stay free
    start_values = settings.parameterize(model.interpretation).values()
    priors = []
    for key, names in model.coupling_names.items():
        value = start_values[names[0]]
        if abs(value) >= 50.0:
            continue
        sigma = args.prior_sigma_hh if key.startswith("J(H") else args.prior_sigma_ch
        if sigma > 0:
            priors += [(n, value, sigma) for n in names]
    settings = dataclasses.replace(settings, priors=tuple(priors), prior_weight=args.prior_weight)
    out = Path(args.out or f"runs/processed/{args.id}_staged")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    prefixes = tuple(p for p in args.hold_first.split(",") if p)
    held = {n for key, names in model.coupling_names.items() if key.startswith(prefixes) for n in names}
    exch_labels = [x for x in args.exchange.split(",") if x]
    kex_bounds = tuple(float(v) for v in args.kex_bounds.split(","))
    p1 = settings.parameterize(model.interpretation)
    add_exchange(p1, model, exch_labels, args.kex_start, kex_bounds, args.kex_fixed)
    p1.fix(*[n for n in p1.order if n in held])
    r1 = refine(model.interpretation, obs, settings, parameterization=p1)
    p2 = settings.parameterize(model.interpretation)
    kex_names = add_exchange(p2, model, exch_labels, args.kex_start, kex_bounds, args.kex_fixed)
    for n, value in r1.parameters.items():
        if n in p2.parameters:
            p2.set(n, value)
    x0 = p2.vector()
    lower, upper = p2.bounds()
    coupling = np.array([p2.parameters[n].kind == "coupling" for n in p2.free_names])
    rng = np.random.default_rng(args.seed)
    starts = [x0] + [np.clip(x0 + np.where(coupling, rng.normal(0, args.spread, len(x0)), 0.0), lower, upper)
                     for _ in range(max(args.starts - 1, 0))]
    if kex_names and not args.kex_fixed:
        # exchange-rate starts spread over decades (the rate is not known in advance)
        col = p2.free_names.index(kex_names[0])
        for k in (1.0, 100.0, 3000.0):
            x = x0.copy()
            x[col] = np.clip(np.log(k), lower[col], upper[col])
            starts.append(x)
    r2 = refine(model.interpretation, obs, settings, parameterization=p2, initial_points=starts)
    e = Evaluated(model.name, args.variant, "staged", model, r2.summary(), r2.prediction,
                  k=free_parameter_count(model, settings))
    e.couplings = model.named_couplings(r2.parameters)
    row = {"id": args.id, "model": model.name, "variant": args.variant, "seconds": round(time.time() - t0),
           "hold_first": list(prefixes),
           "priors": {"sigma_hh": args.prior_sigma_hh, "sigma_ch": args.prior_sigma_ch, "weight": args.prior_weight,
                      "centres": {k: round(start_values[n[0]], 3) for k, n in model.coupling_names.items()
                                  if abs(start_values[n[0]]) < 50.0}},
           "stage1": {"data_region_residual": r1.data_region_residual,
                      "couplings": {k: round(x, 3) for k, x in model.named_couplings(r1.parameters).items()}},
           "stage2": {"data_region_residual": r2.data_region_residual, "relative_residual": r2.relative_residual,
                      "flags": r2.flags, "boundary_hits": r2.boundary_hits,
                      "couplings": {k: round(x, 3) for k, x in e.couplings.items()},
                      "rates_per_s": {n: float(np.exp(x)) for n, x in r2.parameters.items() if "log_rate" in n},
                      "exchange_rates_per_s": {n.split("log_kex_")[1]: float(np.exp(r2.parameters[n]))
                                               for n in kex_names},
                      "delay_ms": 1e3 * r2.parameters.get("phase_delay", 0.0)},
           "k": e.k}
    try:
        if kex_names:
            raise NotImplementedError("linearised uncertainties do not include the exchange model yet")
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
    run_log.finish({"residual_data_cores": r2.data_region_residual, "k": e.k, "model": model.name,
                    "stage1_residual": r1.data_region_residual, "flags": r2.flags,
                    "boundary_hits": r2.boundary_hits, "exchange_rates_per_s": row["stage2"]["exchange_rates_per_s"],
                    "couplings": row["stage2"]["couplings"]}, figures=["fit.png"])
    print(json.dumps({"id": args.id, "residual": round(r2.data_region_residual, 4), "k": e.k,
                      "seconds": row["seconds"], "near_best": row.get("start_agreement", {}).get("near_best"),
                      "starts": row.get("start_agreement", {}).get("starts")}), flush=True)


if __name__ == "__main__":
    main()
