"""Scan the small couplings of a finished fit against the residual in one frequency window.

    python scripts/coupling_scan.py RUN_FOLDER OUT.json [--window 115,140] [--max-abs 20] [--span 8] [--step 0.25]
        [--pairs 3] [--pair-step 0.5]

Everything else (rates, phase, delay, field) stays at the fit's final vector; the linear gains are solved again at
every point (the fit's own forward model and objective, scripts/fit_monitor.py). For each coupling with |J| below
`--max-abs` Hz: a 1D scan over J0 +- span; the window residual (sum of |data - model|^2 inside `--window`) and the
whole objective at every point. Then 2D grids for the `--pairs` couplings whose 1D scans lower the window residual
most, all pairs among them. Writes the scans and the best points; prints a summary. Conditional numerical results
of the given structure (rates held), not a new fit.
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from fit_monitor import SpectrumCache  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--window", default="115,140")
    ap.add_argument("--max-abs", type=float, default=20.0)
    ap.add_argument("--span", type=float, default=8.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--pairs", type=int, default=3)
    ap.add_argument("--pair-step", type=float, default=0.5)
    args = ap.parse_args()
    run = str(Path(args.run))
    fit = json.load(open(Path(run) / "fit.json"))
    joint = SpectrumCache().get(run).joint
    z0 = np.array(fit["z_final"], float)
    fw = joint.forwards[0]
    f = np.asarray(fw.f)
    lo, hi = (float(v) for v in args.window.split(","))
    win = (f > lo) & (f < hi)
    names = list(joint.coupling)

    def value(z, k):
        return float(z[k * joint.m])

    def setv(z, k, v):
        z[k * joint.m:(k + 1) * joint.m] = v

    def evaluate(z):
        pred = fw.predict(joint.spectrum_vector(z, 0))
        res = np.asarray(pred.model) - np.asarray(fw.y)
        r = joint.residual(z)
        return float(np.sum(np.abs(res[win]) ** 2)), float(r @ r)

    t0 = time.time()
    base_win, base_obj = evaluate(z0)
    small = [k for k, n in enumerate(names) if abs(value(z0, k)) < args.max_abs]
    print(f"{len(small)} small couplings of {len(names)}; window residual {base_win:.4g}, objective {base_obj:.4g} "
          f"({time.time() - t0:.1f} s per point)", flush=True)
    one = {}
    for k in small:
        j0 = value(z0, k)
        grid = np.round(np.arange(j0 - args.span, j0 + args.span + 1e-9, args.step), 4)
        rows = []
        for v in grid:
            z = z0.copy()
            setv(z, k, v)
            rows.append(evaluate(z))
        rows = np.array(rows)
        i = int(np.argmin(rows[:, 0]))
        one[names[k]] = {"start": j0, "grid": grid.tolist(), "window": rows[:, 0].tolist(),
                         "objective": rows[:, 1].tolist(), "best": float(grid[i]),
                         "best_window": float(rows[i, 0]), "best_objective": float(rows[i, 1])}
        print(f"  {names[k]:18s} start {j0:7.2f}  best {grid[i]:7.2f}  window {base_win:.4g} -> {rows[i, 0]:.4g}  "
              f"objective -> {rows[i, 1]:.4g}  ({time.time() - t0:.0f} s)", flush=True)
    ranked = sorted(one, key=lambda n: one[n]["best_window"])[: args.pairs]
    two = {}
    for a, b in itertools.combinations(ranked, 2):
        ka, kb = names.index(a), names.index(b)
        ga = np.round(np.arange(value(z0, ka) - args.span, value(z0, ka) + args.span + 1e-9, args.pair_step), 4)
        gb = np.round(np.arange(value(z0, kb) - args.span, value(z0, kb) + args.span + 1e-9, args.pair_step), 4)
        win_grid = np.zeros((len(ga), len(gb)))
        obj_grid = np.zeros_like(win_grid)
        for i, va in enumerate(ga):
            for j, vb in enumerate(gb):
                z = z0.copy()
                setv(z, ka, va)
                setv(z, kb, vb)
                win_grid[i, j], obj_grid[i, j] = evaluate(z)
        i, j = np.unravel_index(int(np.argmin(win_grid)), win_grid.shape)
        two[f"{a} x {b}"] = {"a": a, "b": b, "grid_a": ga.tolist(), "grid_b": gb.tolist(),
                             "window": win_grid.tolist(), "objective": obj_grid.tolist(),
                             "best": [float(ga[i]), float(gb[j])], "best_window": float(win_grid[i, j]),
                             "best_objective": float(obj_grid[i, j])}
        print(f"  2D {a} x {b}: best {ga[i]:.2f}, {gb[j]:.2f}  window -> {win_grid[i, j]:.4g}  objective -> "
              f"{obj_grid[i, j]:.4g}  ({time.time() - t0:.0f} s)", flush=True)
    json.dump({"run": run, "window_hz": [lo, hi], "base_window": base_win, "base_objective": base_obj,
               "one_d": one, "two_d": two, "settings": vars(args)}, open(args.out, "w"), indent=1)


if __name__ == "__main__":
    main()
