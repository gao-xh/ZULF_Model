"""Field sensitivity dnu/d|B| of the model lines of a fit and its fitted decay rates, for comparing two fits of the
same compound (e.g. two sample volumes in one fringe field).

    python scripts/field_sensitivity.py --fit runs/processed/A_fit --fit runs/processed/B_fit --label 2mL --label 250uL \
        [--bands 126,147;262,285] [--figure runs/figures/OUT.png] [--json OUT.json]

dnu/d|B|: the transition frequencies of every isotopologue of the fit's structure at the fit's field and at the
field magnitude scaled by 1 + 1e-3 (direction fixed; zulf_core compute_transitions), finite difference over the change
of |B|, at the field of each fit. Decay rates: the fit's rate of each model line (ZULF Studio session lines(), one rate
per family). Lines of the two fits are matched by frequency (0.08 Hz). Expected for field inhomogeneity:
rate = R0 + pi |dnu/dB| dB, so the difference of two fits is pi |dnu/dB| (dB_A - dB_B). Conditional numerical
results of the fitted model; no uncertainties.
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def fit_lines(run, bands, eps=1e-3, min_relative=0.03):
    """[(frequency, relative, rate, slope Hz/nT)] of the model lines of a fit inside the bands."""
    from zulf_core.physics.protocol import Protocol
    from zulf_core.physics.transitions import compute_transitions
    from zulf_studio.session import StudioSession
    run = Path(run)
    status = json.loads((run / "monitor" / "status.json").read_text())
    argv = status["argv"]
    series = argv[argv.index("--series") + 1]
    s = StudioSession(log_file=False)
    s.load_spectrum(series=series)
    s.apply_fit(str(run))
    fit = json.loads((run / "fit.json").read_text())
    sp = next(iter(fit["spectrum_parameters"].values()))
    field = np.array([sp.get("field_transverse_ut", 0.0), 0.0, sp.get("field_z_ut", 0.0)])
    bmag = 1e3 * float(np.linalg.norm(field))                      # nT
    out = []
    lines = [r for r in s.lines(min_relative) if any(a <= r["frequency_hz"] <= b for a, b in bands)]
    _, model = s.model()                                           # the fitted couplings (apply_fit)
    systems = {lab: c.system for lab, c in zip(model.component_labels, model.interpretation.components)}
    for r in lines:
        out.append([r["frequency_hz"], r["relative"], r["decay_per_s"], None, r["component"]])
    # slopes: transitions of each component's spin system at B and B (1 + eps), matched by order of frequency
    for label in {r[4] for r in out}:
        system = systems[label]
        def freqs(b):
            tl = compute_transitions(system, Protocol(field_ut=tuple(b)))
            f = np.asarray(tl.frequencies_hz)
            a = np.abs(np.asarray(tl.amplitudes))
            keep = a > 1e-6 * a.max()
            return f[keep]
        f0, f1 = freqs(field), freqs(field * (1 + eps))
        for row in out:
            if row[4] != label:
                continue
            i = int(np.argmin(np.abs(f0 - row[0])))
            row[3] = (f1[i] - f0[i]) / (eps * bmag)
    return [tuple(r[:4]) for r in out], bmag


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", action="append", required=True)
    ap.add_argument("--label", action="append", required=True)
    ap.add_argument("--bands", default="126,147;262,285")
    ap.add_argument("--figure", default="")
    ap.add_argument("--json", default="")
    args = ap.parse_args()
    bands = [tuple(float(v) for v in q.split(",")) for q in args.bands.split(";")]
    res = [fit_lines(f, bands) for f in args.fit]
    (la, ba), (lb, bb) = res[0], res[1]
    rows = []
    for f, rel, rate, slope in la:
        match = min(lb, key=lambda x: abs(x[0] - f))
        rate_b = match[2] if abs(match[0] - f) < 0.08 else float("nan")
        diff = rate - rate_b
        db = diff / (math.pi * abs(slope)) if abs(slope) > 1e-3 else float("nan")
        rows.append({"frequency_hz": f, "slope_hz_per_nt": slope, "rate_a": rate, "rate_b": rate_b,
                     "difference": diff, "delta_b_nt": db, "relative_a": rel})
    print(f"|B|: {args.label[0]} {ba:.1f} nT, {args.label[1]} {bb:.1f} nT")
    print(f"{'line (Hz)':>10} {'dnu/dB (Hz/nT)':>15} {args.label[0]:>8} {args.label[1]:>8} {'diff':>7} {'dB (nT)':>8}")
    for r in rows:
        print(f"{r['frequency_hz']:10.2f} {r['slope_hz_per_nt']:+15.4f} {r['rate_a']:8.2f} {r['rate_b']:8.2f} "
              f"{r['difference']:+7.2f} {r['delta_b_nt']:8.1f}")
    if args.json:
        Path(args.json).write_text(json.dumps({"labels": args.label, "fits": args.fit, "rows": rows}, indent=1))
    if args.figure:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        for r in rows:
            band = "J band" if r["frequency_hz"] < 200 else "2J band"
            ax.plot(abs(r["slope_hz_per_nt"]), r["difference"], "o" if band == "J band" else "s",
                    color="#d4682e" if band == "J band" else "#1d1d1f", ms=6)
            ax.annotate(f"{r['frequency_hz']:.1f}", (abs(r["slope_hz_per_nt"]), r["difference"]),
                        textcoords="offset points", xytext=(5, 3), fontsize=7)
        ax.plot([], [], "o", color="#d4682e", label="J band")
        ax.plot([], [], "s", color="#1d1d1f", label="2J band")
        ax.axhline(0, color="0.8", lw=0.6)
        x = np.linspace(0, 0.055, 50)
        for db, ls in ((10, ":"), (20, "--")):
            ax.plot(x, math.pi * x * db, ls, color="0.5", lw=0.8, label=f"pi |dnu/dB| x {db} nT")
        ax.set_xlabel("|dnu/dB| (Hz/nT)")
        ax.set_ylabel(f"decay rate {args.label[0]} - {args.label[1]} (1/s)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(args.figure, dpi=150)


if __name__ == "__main__":
    main()
