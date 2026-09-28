"""Known-structure fit of a spectrum that was processed elsewhere (frequency and value arrays, real part).

    python scripts/fit_processed_spectrum.py --freq F.npy --values V.npy --id NAME \\
        --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162, "A4": 162}}' \\
        [--couplings '{"J(HA2,HA3)": 4.9, ...}'] [--range 140,200] [--out runs/processed/NAME] [--workers 4]

The spectrum is compared on its real part with ideal Lorentzian (infinite-record) lines
(ObservedSpectrum.from_spectrum with record=None); a residual common phase and a delay are fitted with the
couplings, as for any processed spectrum. Without a FID there is no global pattern search, so the starting
couplings matter: `--couplings` overrides fragment couplings by key (symmetry partners follow). Writes the
report of zulf_hypothesis.write_report (OUT/structure.{json,md,png}) and OUT/fit.json.
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
from zulf_core.solver import ObservedSpectrum                     # noqa: E402
from zulf_hypothesis import fit_settings, fit_structure            # noqa: E402
from zulf_hypothesis.fragment import pair                          # noqa: E402


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
    ap.add_argument("--out", default="")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    import regression_confirmed as reg
    f, v = np.load(args.freq).astype(float), np.load(args.values).astype(float)
    lo, hi = (float(x) for x in args.range.split(","))
    obs = ObservedSpectrum.from_spectrum(f, v, [(lo, hi)], record=None, real_only=True, label=args.id)
    spec = json.loads(args.structure)
    spec.setdefault("compound", args.id)
    fragment = override_couplings(reg.structure_for(spec), json.loads(args.couplings))
    out = Path(args.out or f"runs/processed/{args.id}")
    out.mkdir(parents=True, exist_ok=True)
    t = time.time()
    fit = fit_structure(fragment, obs, fit_settings(workers=args.workers), report=str(out / "structure"))
    best = fit.best
    row = {"id": args.id, "best": best.key, "reduced_chi2": round(best.chi2 / max(best.n - best.k, 1), 2),
           "k": best.k, "seconds": round(time.time() - t), "couplings": {k: round(float(x), 3) for k, x in best.couplings.items()},
           "findings": [x.code for x in best.findings], "table": fit.table()[:8]}
    json.dump(row, open(out / "fit.json", "w"), indent=1, default=str)
    print(json.dumps({k: row[k] for k in ("id", "best", "reduced_chi2", "k", "seconds")}), flush=True)


if __name__ == "__main__":
    main()
