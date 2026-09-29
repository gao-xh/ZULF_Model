"""Synthetic recovery test for a concentration series: can the joint fit recover known couplings from spectra
like the measured ones?

    python scripts/synthetic_recovery_series.py --series series.json --truth JOINT/fit.json \\
        --structure '{...}' --couplings '{...}' [--signal-threshold 2.5 --signal-taper 4] [--seed 0]
        --out runs/processed/synthetic_series

The truth is a joint fit (couplings at every concentration, and every spectrum's rates and delay). Each
synthetic spectrum is that model rendered on the measured frequency grid (same components, gains and phase as fitted
to the measured spectrum) plus white noise with the measured spectrum's noise level (the robust sigma of the
signal weighting). Writes the synthetic spectra, OUT/series.json for fit_joint_series.py and OUT/truth.json.
If the joint fit of these spectra (same settings, starts from the accepted values) does not return the truth, the
couplings are not determined by spectra of this kind, whatever the measured data contain.
"""
import argparse
import dataclasses
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from fit_joint_series import JointSeries                                 # noqa: E402
from zulf_core.solver import ObservedSpectrum                            # noqa: E402
from zulf_hypothesis.builder import build_model                          # noqa: E402
from zulf_hypothesis.fit import default_fit_base                         # noqa: E402
from zulf_hypothesis.search import _settings_for                         # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--structure", required=True)
    ap.add_argument("--couplings", default="{}")
    ap.add_argument("--variant", default="ratios")
    ap.add_argument("--range", default="140,200")
    ap.add_argument("--signal-threshold", type=float, default=0.0)
    ap.add_argument("--signal-taper", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    import regression_confirmed as reg
    from fit_processed_spectrum import override_couplings
    series = json.load(open(args.series))
    truth = json.load(open(args.truth))
    lo, hi = (float(v) for v in args.range.split(","))
    obs = [ObservedSpectrum.from_spectrum(np.load(e["freq"]).astype(float), np.load(e["values"]).astype(float),
                                          [(lo, hi)], record=None, real_only=True, label=e["id"]) for e in series]
    spec = json.loads(args.structure)
    spec.setdefault("compound", "series")
    model = build_model(override_couplings(reg.structure_for(spec), json.loads(args.couplings)), ranges=[(lo, hi)])
    base = default_fit_base()
    if args.signal_threshold:
        base = dataclasses.replace(base, signal_threshold=args.signal_threshold)
    if args.signal_taper:
        base = dataclasses.replace(base, signal_taper_hz=args.signal_taper)
    settings = _settings_for(model, args.variant, base)
    joint = JointSeries(model, settings, obs, [e["x"] for e in series])
    key_of = {}
    for key, names in model.coupling_names.items():
        for n in names:
            leader = joint.params[0].ties.get(n, n)
            if leader in joint.coupling and leader not in key_of:
                key_of[leader] = key
    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for s, (e, f) in enumerate(zip(series, joint.forwards)):
        x = joint.params[s].vector()
        node = int(np.argmin(np.abs(np.asarray(truth["x"]) - e["x"])))
        for n in joint.coupling:
            x[joint.col[n]] = truth["couplings"][key_of[n]]["J_at_x"][node]
        for n, v in truth.get("spectrum_parameters", {}).get(e["id"], {}).items():
            if n in joint.col:
                x[joint.col[n]] = v
        clean = np.asarray(f.predict(x).model).real
        sigma = float(f.noise_sigma) if getattr(f, "noise_sigma", None) is not None else float(
            1.4826 * np.median(np.abs(f.y.real - np.median(f.y.real))))
        values = clean + rng.normal(0.0, sigma, len(clean))
        np.save(out / f"{e['id']}_frequency.npy", f.f)
        np.save(out / f"{e['id']}_values.npy", values)
        rows.append({"id": e["id"], "x": e["x"], "freq": str(out / f"{e['id']}_frequency.npy"),
                     "values": str(out / f"{e['id']}_values.npy"), "noise_sigma": sigma})
    json.dump(rows, open(out / "series.json", "w"), indent=1)
    json.dump({"source": args.truth, "x": truth["x"],
               "couplings": {k: v["J_at_x"] for k, v in truth["couplings"].items()},
               "spectrum_parameters": truth.get("spectrum_parameters", {})}, open(out / "truth.json", "w"), indent=1)
    print(json.dumps({"spectra": len(rows), "noise_sigma": [round(r["noise_sigma"], 3) for r in rows]}))


if __name__ == "__main__":
    main()
