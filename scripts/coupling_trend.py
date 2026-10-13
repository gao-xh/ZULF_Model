"""Couplings of a series of fits against a series variable (e.g. the pyridine mole fraction), one figure.

    python scripts/coupling_trend.py OUT.png --x 2,10,33 --runs RUN1,RUN2,RUN3 [--xlabel "pyridine (mol %)"]
        [--title TEXT] [--csv OUT.csv]

Each fit's couplings (fit.json J_at_x) with their linearised standard errors (J_std_at_x); one-bond couplings
(|J| > 60 Hz), long-range heteronuclear and homonuclear couplings in separate panels. The linearised errors are
known to be 10-30 times too small (docs/DECISIONS.md); they are drawn as given and labelled so. Conditional
numerical results of the given structure, not determined assignments.
"""
import argparse
import csv
import json
import re
from pathlib import Path

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--x", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--xlabel", default="series variable")
    ap.add_argument("--title", default="")
    ap.add_argument("--csv", default="")
    args = ap.parse_args()
    xs = [float(v) for v in args.x.split(",")]
    runs = [Path(r) for r in args.runs.split(",")]
    fits = [json.load(open(r / "fit.json"))["couplings"] for r in runs]
    keys = sorted(set().union(*fits))
    values = {k: np.array([f[k]["J_at_x"][0] if k in f else np.nan for f in fits]) for k in keys}
    errors = {k: np.array([f[k].get("J_std_at_x", [np.nan])[0] if k in f else np.nan for f in fits]) for k in keys}

    def kind(k):
        a, b = re.match(r"J\(([^,]+),([^)]+)\)", k).groups()
        if np.nanmedian(np.abs(values[k])) > 60:
            return 0
        return 2 if a.startswith("H") and b.startswith("H") else 1

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    titles = ["one-bond 1J(C,H)", "long-range J(C,H)", "J(H,H)"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))
    cmap = plt.get_cmap("tab10")
    for p in range(3):
        ax = axes[p]
        group = [k for k in keys if kind(k) == p]
        for i, k in enumerate(group):
            style = "-o" if i < 10 else "--s"
            ax.errorbar(xs, values[k], yerr=errors[k], fmt=style, ms=4, lw=1.2, capsize=2,
                        color=cmap(i % 10), label=k)
        ax.set_title(titles[p], fontsize=11)
        ax.set_xlabel(args.xlabel)
        ax.set_ylabel("J (Hz)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7, frameon=False, ncol=2 if len(group) > 6 else 1)
    fig.suptitle((args.title + "; " if args.title else "") + "error bars: linearised standard errors "
                 "(known to be 10-30x too small)", fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out, dpi=120)
    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["coupling"] + [f"J at {x:g}" for x in xs] + [f"std at {x:g}" for x in xs])
            for k in keys:
                w.writerow([k] + [f"{v:.4f}" for v in values[k]] + [f"{e:.4f}" for e in errors[k]])
    print(args.out)


if __name__ == "__main__":
    main()
