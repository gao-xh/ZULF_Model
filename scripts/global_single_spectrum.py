"""Global search of the couplings of one processed spectrum (differential evolution), then a local fit.

    python scripts/global_single_spectrum.py --freq F.npy --values V.npy --id NAME --structure '{...}'
        --couplings '{...}' [--from-joint JOINT/fit.json] [--half-width 5] [--popsize 10] [--maxiter 60]
        [--workers 4] [--prior-sigma-hh 0.5 --prior-sigma-ch 1.5 --prior-weight 10]
        [--signal-threshold 2.5 --signal-taper 4] [--out DIR]

The coupling bounds are the starting value +- half_width (Hz). Decay rates and the delay are held at the values of
--from-joint (this spectrum's id) or the defaults during the global stage; the gains, shared phase and
background are solved inside every evaluation (variable projection). The best member is then refined by least
squares with every parameter free. Motivation (ANALYSIS_LOG): on synthetic spectra like the Blake pyridine series
the right minimum is about 1 Hz wide in the coupling levels, so random local starts rarely reach it.
"""
import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from zulf_core.physics.protocol import SUDDEN_DROP                      # noqa: E402
from zulf_core.solver import ObservedSpectrum                           # noqa: E402
from zulf_core.solver.forward import MixtureForward                     # noqa: E402
from zulf_core.solver.refine import _signal_kwargs                      # noqa: E402
from zulf_hypothesis.builder import build_model                         # noqa: E402
from zulf_hypothesis.fit import default_fit_base                        # noqa: E402
from zulf_hypothesis.search import _settings_for                        # noqa: E402

_TASK = None


def _objective(xc):
    fw, x_full, idx, prior_idx, prior_mean, prior_scale = _TASK
    x = x_full.copy()
    x[idx] = xc
    r = fw.predict(x).residual
    extra = (x[prior_idx] - prior_mean) * prior_scale
    return float(r @ r + extra @ extra)


def main():
    global _TASK
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required=True)
    ap.add_argument("--values", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--from-joint", default="")
    ap.add_argument("--half-width", type=float, default=5.0)
    ap.add_argument("--popsize", type=int, default=10)
    ap.add_argument("--maxiter", type=int, default=60)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--prior-sigma-hh", type=float, default=0.0)
    ap.add_argument("--prior-sigma-ch", type=float, default=0.0)
    ap.add_argument("--prior-weight", type=float, default=1.0)
    ap.add_argument("--signal-threshold", type=float, default=0.0)
    ap.add_argument("--signal-taper", type=float, default=0.0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    lo, hi = (float(v) for v in args.range.split(","))
    obs = ObservedSpectrum.from_spectrum(np.load(args.freq).astype(float), np.load(args.values).astype(float),
                                         [(lo, hi)], record=None, real_only=True, label=args.id)
    spec = json.loads(args.structure)
    spec.setdefault("compound", args.id)
    model = build_model(override_couplings(reg.structure_for(spec), json.loads(args.couplings)), ranges=[(lo, hi)])
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    settings = _settings_for(model, args.variant, base)
    param = settings.parameterize(model.interpretation)
    fw = MixtureForward(param, obs, SUDDEN_DROP, settings.gain_model, settings.background_order,
                        settings.band_weighting, **_signal_kwargs(settings))
    free = param.free_names
    col = {n: i for i, n in enumerate(free)}
    x0 = param.vector()
    if args.from_joint:
        prev = json.load(open(args.from_joint))
        for n, v in prev.get("spectrum_parameters", {}).get(args.id, {}).items():
            if n in col:
                x0[col[n]] = v
    lower, upper = param.bounds()
    idx = np.array([col[n] for n in free if param.parameters[n].kind == "coupling"])
    centre = x0[idx]
    bounds = [(max(c - args.half_width, lower[i]), min(c + args.half_width, upper[i])) for c, i in zip(centre, idx)]
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = param.ties.get(n, n)
            if leader in col and leader not in key_of:
                key_of[leader] = key
    pi, pm, ps = [], [], []
    for i in idx:
        n = free[i]
        if abs(x0[i]) < 50.0:
            sd = args.prior_sigma_hh if key_of.get(n, n).startswith("J(H") else args.prior_sigma_ch
            if sd > 0:
                pi.append(i)
                pm.append(x0[i])
                ps.append(np.sqrt(args.prior_weight) / sd / fw.norm)
    _TASK = (fw, x0, idx, np.array(pi, int), np.array(pm), np.array(ps))
    t0 = time.time()
    if args.workers > 1:
        import multiprocessing
        pool = multiprocessing.get_context("fork").Pool(args.workers)
        de = differential_evolution(_objective, bounds, popsize=args.popsize, maxiter=args.maxiter, seed=args.seed,
                                    polish=False, workers=pool.map, updating="deferred", tol=1e-8)
        pool.close()
    else:
        de = differential_evolution(_objective, bounds, popsize=args.popsize, maxiter=args.maxiter, seed=args.seed,
                                    polish=False, tol=1e-8)
    t_de = time.time() - t0
    x = x0.copy()
    x[idx] = de.x
    prior_idx, prior_mean, prior_scale = np.array(pi, int), np.array(pm), np.array(ps)

    def residual(xx):
        return np.r_[fw.predict(xx).residual, (xx[prior_idx] - prior_mean) * prior_scale]

    def jacobian(xx):
        rows = np.zeros((len(prior_idx), len(xx)))
        rows[np.arange(len(prior_idx)), prior_idx] = prior_scale
        return np.vstack([fw.jacobian(xx), rows])

    sol = least_squares(residual, np.clip(x, lower + 1e-9, upper - 1e-9), jac=jacobian, bounds=(lower, upper),
                        x_scale="jac", max_nfev=300, ftol=1e-10, xtol=1e-10, gtol=1e-10)
    values = param.values(sol.x)
    result = {"id": args.id, "de_score": float(de.fun), "de_evaluations": int(de.nfev), "de_seconds": round(t_de),
              "score": float(2 * sol.cost), "couplings": {k: float(values[names[0]])
                                                           for k, names in model.coupling_names.items()},
              "spectrum_parameters": {n: float(v) for n, v in zip(free, sol.x) if n not in set(free[i] for i in idx)}}
    out = Path(args.out or f"runs/processed/global_{args.id}")
    out.mkdir(parents=True, exist_ok=True)
    json.dump(result, open(out / "fit.json", "w"), indent=1)
    print(json.dumps({k: result[k] for k in ("id", "de_score", "de_evaluations", "de_seconds", "score")}), flush=True)


if __name__ == "__main__":
    main()
