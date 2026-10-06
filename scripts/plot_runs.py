"""Compare finished fit_joint_series runs of one sample: data and model (real part) over the whole fit range,
zoom windows below with their residuals.

    python scripts/plot_runs.py OUT.png RUN_DIR [RUN_DIR ...] [--zooms "lo,hi;lo,hi"] [--title TEXT]

Each run is rebuilt from its RUN_LOG.md (or monitor/status.json) with its own fit.json (run_problem.load_run).
The legend gives the fit.json score and the data-region residual. When the run used model-line passes, the
score is under the final weights; the objective under the original weights is printed as well.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_problem import load_run, segments, spectrum_arrays  # noqa: E402

COLORS = ["#c0469e", "#2a78d6", "#2f8f5b", "#eb6834", "#7b2f8f", "#8a6d1d"]


def summary(run):
    fit = run.fit or {}
    score = (fit.get("scores") or [float("nan")])[0]
    residual = next(iter((fit.get("data_region_residuals") or {"": float("nan")}).values()))
    original = (fit.get("model_line_passes") or {}).get("objective_original_weights")
    return score, residual, original


def plot(out, run_dirs, zooms=(), title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    runs = [load_run(r) for r in run_dirs]
    arrays = [spectrum_arrays(r.prob) for r in runs]
    ncol = max(len(zooms), 1)
    rows = 3 if zooms else 2
    fig = plt.figure(figsize=(4.2 * max(ncol, 3), 2.6 * rows))
    gs = fig.add_gridspec(rows, ncol, height_ratios=([1.2, 1.2, 0.8] if zooms else [1.2, 0.6]), hspace=0.45,
                          wspace=0.25)
    ax0 = fig.add_subplot(gs[0, :])

    def line(ax, x, v, **kw):
        for seg in segments(x):
            ax.plot(x[seg], v[seg], **kw)
            kw.pop("label", None)
    f, y = arrays[0][0], arrays[0][1]
    line(ax0, f, y.real, color="black", lw=0.6, label="data")
    rows_out = []
    for k, (run, (ff, yy, m, _, _)) in enumerate(zip(runs, arrays)):
        score, residual, original = summary(run)
        extra = f", original weights {original:.4f}" if original is not None else ""
        line(ax0, ff, m.real, color=COLORS[k % len(COLORS)], lw=0.6, alpha=0.85,
             label=f"{run.name}: objective {score:.4f}{extra}, residual {residual:.3f}")
        rows_out.append((run.name, score, original, residual))
    ax0.legend(fontsize=7, loc="upper left", frameon=False)
    ax0.set_title(title or ", ".join(r.name for r in runs), fontsize=10, loc="left")
    ax0.set_xlabel("Frequency [Hz]", fontsize=8)
    ax0.tick_params(labelsize=7)
    if not zooms:
        axr = fig.add_subplot(gs[1, :], sharex=ax0)
        for k, (ff, yy, m, _, _) in enumerate(arrays):
            line(axr, ff, (yy - m).real, color=COLORS[k % len(COLORS)], lw=0.6, alpha=0.85)
        axr.set_ylabel("residual", fontsize=7)
        axr.tick_params(labelsize=7)
    for jz, (lo, hi) in enumerate(zooms):
        ax = fig.add_subplot(gs[1, jz])
        axr = fig.add_subplot(gs[2, jz], sharex=ax)
        w = (f >= lo) & (f <= hi)
        line(ax, f[w], y.real[w], color="black", lw=0.8)
        for k, (ff, yy, m, _, _) in enumerate(arrays):
            ww = (ff >= lo) & (ff <= hi)
            line(ax, ff[ww], m.real[ww], color=COLORS[k % len(COLORS)], lw=0.8, alpha=0.85)
            line(axr, ff[ww], (yy - m).real[ww], color=COLORS[k % len(COLORS)], lw=0.7, alpha=0.85)
        ax.set_xlim(lo, hi)
        ax.set_title(f"{lo:g}-{hi:g} Hz", fontsize=8)
        ax.tick_params(labelsize=7)
        axr.tick_params(labelsize=7)
        if jz == 0:
            axr.set_ylabel("residual", fontsize=7)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return rows_out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("runs", nargs="+", help="run output directories")
    ap.add_argument("--zooms", default="", help="lo,hi;lo,hi windows (Hz)")
    ap.add_argument("--title", default="")
    args = ap.parse_args()
    zooms = [tuple(float(v) for v in z.split(",")) for z in args.zooms.split(";") if z.strip()]
    for name, score, original, residual in plot(args.out, args.runs, zooms, args.title):
        extra = f"  (original weights {original:.5f})" if original is not None else ""
        print(f"{name}: objective {score:.5f}{extra}, data residual {residual:.3f}")
    print(args.out)


if __name__ == "__main__":
    main()
