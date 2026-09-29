"""Known-structure fit of a spectrum that was processed elsewhere (frequency and value arrays, real part).

    python scripts/fit_processed_spectrum.py --freq F.npy --values V.npy --id NAME \\
        --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162, "A4": 162}}' \\
        [--couplings '{"J(HA2,HA3)": 4.9, ...}'] [--range 140,200] [--record 4000,88457]
        [--out runs/processed/NAME] [--workers 4]

The spectrum is compared on its real part. `--record fs,points` gives the acquisition it came from (a plain FFT of
`points` samples at `fs`; e.g. a frequency axis k fs / points): the model is then rendered on that finite record
with analytic derivatives for every line shape (with record=None ideal Lorentzian lines are used and fits with a
Gaussian width fall back to finite differences, several times slower). A residual common phase and a delay are fitted with the
couplings, as for any processed spectrum. Without a FID there is no global pattern search, so the starting
couplings matter: `--couplings` overrides fragment couplings by key (symmetry partners follow). Writes the
report of zulf_hypothesis.write_report (OUT/structure.{json,md,png}), OUT/fit.json and the linearised coupling
uncertainties of the best fit (OUT/uncertainty.md; zulf_hypothesis.uncertainty).
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
from zulf_core.render.acquisition import Acquisition            # noqa: E402
from zulf_core.solver import ObservedSpectrum                     # noqa: E402
from zulf_hypothesis import fit_settings, fit_structure            # noqa: E402
from zulf_hypothesis.fit import protonated                        # noqa: E402
from zulf_hypothesis.fragment import pair                          # noqa: E402
from zulf_hypothesis.uncertainty import coupling_uncertainties, start_agreement, uncertainty_markdown  # noqa: E402


def override_couplings(fragment, couplings):
    """Set couplings by key 'J(a,b)'; every symmetry image of the pair gets the same value."""
    for key, value in couplings.items():
        a, b = key[2:-1].split(",")
        orbit = {pair(a, b)}
        grown = True
        while grown:
            grown = False
            for g in fragment.symmetry:
                for p in list(orbit):
                    x, y = tuple(p)
                    q = pair(g.get(x, x), g.get(y, y))
                    if q not in orbit:
                        orbit.add(q)
                        grown = True
        for p in orbit:
            fragment.couplings[p] = float(value)
    fragment.validate()
    return fragment


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", required=True)
    ap.add_argument("--values", required=True)
    ap.add_argument("--id", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--record", default="", help="fs,points of the acquisition (plain FFT, no processing)")
    ap.add_argument("--exchange", default="auto", choices=["auto", "fast", "slow", "both"],
                    help="exchange regimes (fit_structure); 'fast' fits the given protonation state only")
    ap.add_argument("--protonate", action="store_true", help="add one proton on every N that has room first")
    ap.add_argument("--moves", default="", help="comma-separated refining moves (default: fit_settings' "
                    "REFINING_MOVES); e.g. free_remote_couplings to skip the Gaussian width, whose derivatives "
                    "need --record and are slow for long records")
    ap.add_argument("--starts", type=int, default=0, help="refinement starts per model (default: fit_settings)")
    ap.add_argument("--spread", type=float, default=0.0, help="coupling spread of the perturbed starts, Hz")
    ap.add_argument("--delay-bounds", default="", help="lo,hi in s for the fitted delay (e.g. -1e-5,1e-5 for a "
                    "spectrum already phase-corrected)")
    ap.add_argument("--variants", default="", help="comma-separated fit variants (default: fit_settings)")
    ap.add_argument("--out", default="")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    from run_log import RunLog
    run_log = RunLog(Path(args.out or f"runs/processed/{args.id}"), "fit_processed_spectrum")
    import regression_confirmed as reg
    f, v = np.load(args.freq).astype(float), np.load(args.values).astype(float)
    lo, hi = (float(x) for x in args.range.split(","))
    record = None
    if args.record:
        fs, points = args.record.split(",")
        record = Acquisition.pure(float(fs), int(points))
    obs = ObservedSpectrum.from_spectrum(f, v, [(lo, hi)], record=record, real_only=True, label=args.id)
    spec = json.loads(args.structure)
    spec.setdefault("compound", args.id)
    fragment = reg.structure_for(spec)
    if args.protonate:
        fragment = protonated(fragment)
    fragment = override_couplings(fragment, json.loads(args.couplings))
    out = Path(args.out or f"runs/processed/{args.id}")
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    changes = {"moves": tuple(m for m in args.moves.split(",") if m)} if args.moves else {}
    if args.variants:
        changes["variants"] = tuple(v for v in args.variants.split(",") if v)
    delay = tuple(float(x) for x in args.delay_bounds.split(",")) if args.delay_bounds else None
    settings = fit_settings(workers=args.workers, phase_delay_bounds_s=delay, **changes)
    base = settings.base
    if args.starts:
        # every start gets the time budget one start had (fit_settings: 900 s for 4 starts)
        base = dataclasses.replace(base, starts=args.starts,
                                   max_seconds=max(base.max_seconds, base.max_seconds * args.starts / 4))
    if args.spread:
        base = dataclasses.replace(base, start_spread_hz=args.spread)
    settings = dataclasses.replace(settings, base=base)
    fit = fit_structure(fragment, obs, settings, exchange=args.exchange,
                        report=str(out / "structure"))
    best = fit.best
    row = {"id": args.id, "best": best.key, "reduced_chi2": round(best.chi2 / max(best.n - best.k, 1), 2),
           "k": best.k, "seconds": round(time.time() - t), "couplings": {k: round(float(x), 3) for k, x in best.couplings.items()},
           "findings": [x.code for x in best.findings], "table": fit.table()[:8]}
    try:
        u = coupling_uncertainties(best, obs, settings.base)
        row["uncertainty"] = u
        agreement = start_agreement(best, settings.base)
        row["start_agreement"] = agreement
        (out / "uncertainty.md").write_text(
            uncertainty_markdown(u) + f"\n\nStarts: {agreement.get('starts')}, within 1 % of the best score: "
            f"{agreement.get('near_best')}\n")
    except Exception as exc:                  # errors are a report item, not a reason to lose the fit
        row["uncertainty_error"] = f"{type(exc).__name__}: {exc}"
    json.dump(row, open(out / "fit.json", "w"), indent=1, default=str)
    run_log.finish({k: row.get(k) for k in ("id", "best", "reduced_chi2", "k", "seconds")})
    print(json.dumps({k: row[k] for k in ("id", "best", "reduced_chi2", "k", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
