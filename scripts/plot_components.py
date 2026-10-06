"""The components (isotopologues, species) of a finished fit_joint_series run, one panel each, and what every
component has to explain.

    python scripts/plot_components.py RUN_DIR OUT.png [--band lo,hi] [--peaks 0.15]

Panels, real part, one common scale:
a  data and the total model;
b  each component alone, at its fitted amplitude, with its couplings to the own 13C (site from the label);
c  for each component: data minus the other components and the background, overlaid on that component; its
   peaks are marked (local maxima above --peaks x its largest value);
d  residual of the total model.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_problem import load_run, segments, spectrum_arrays  # noqa: E402

COLORS = ["#c0469e", "#2a78d6", "#2f8f5b", "#eb6834", "#7b2f8f", "#8a6d1d"]


def component_peaks(x, v, fraction):
    from scipy.signal import find_peaks
    a = np.abs(v).max() if len(v) else 0.0
    if not a > 0:
        return []
    pk, _ = find_peaks(v, prominence=0.08 * a, height=fraction * a)
    return [float(x[k]) for k in pk]


def plot(run_dir, out, band=None, fraction=0.15):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    run = load_run(run_dir)
    prob = run.prob
    j = prob.joint
    f, y, m, comps, pred = spectrum_arrays(prob)
    labels = list(prob.model.component_labels)
    lo, hi = band if band else (float(f.min()), float(f.max()))
    w = (f >= lo) & (f <= hi)
    background = m - sum(comps)
    gains = np.abs(np.asarray(pred.gains))
    n = len(comps)
    fig, axes = plt.subplots(2 + 2 * n, 1, figsize=(12, 2.3 * (2 + 2 * n) + 1), sharex=True,
                             gridspec_kw={"height_ratios": [1.3] + [1.0] * n + [1.0] * n + [0.6], "hspace": 0.15})

    def line(ax, v, **kw):
        xs, vs = f[w], v[w]
        for seg in segments(xs):
            ax.plot(xs[seg], vs[seg], **kw)
            kw.pop("label", None)
    ymax = float(np.abs(y.real[w]).max()) * 1.08 or 1.0
    score = ((run.fit or {}).get("scores") or [float("nan")])[0]
    line(axes[0], y.real, color="black", lw=0.9, label="data")
    line(axes[0], m.real, color="#7b2f8f", lw=0.9, label="total model")
    axes[0].set_title(f"{run.name}: data and the model of {n} components (objective {score:.4g})", loc="left",
                      fontsize=10)
    key_of = prob.key_of
    out_rows = []
    for c in range(n):
        ax = axes[1 + c]
        color = COLORS[c % len(COLORS)]
        ratio = gains[c] / gains.max() if gains.max() > 0 else float("nan")
        line(ax, comps[c].real, color=color, lw=1.1, label=f"{labels[c]} (amplitude {ratio:.2f} of the largest)")
        prefix = labels[c].split(":")[0] + ":" if ":" in labels[c].split("@")[0] else ""
        site = labels[c].split("@")[-1].split(" ")[0] if "@" in labels[c] else None
        own = []
        for k, name in enumerate(j.coupling):
            key = key_of.get(name, name)
            if site and key.startswith(f"{prefix}J({site},"):
                own.append(f"{key[len(prefix):]} {float(j.coupling_values(prob.z0, k)[0]):.2f}")
        ax.text(0.99, 0.95, "\n".join([", ".join(own[i:i + 3]) for i in range(0, len(own), 3)]),
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=color)
        ax.set_title(f"b{c + 1}  {labels[c]} alone", loc="left", fontsize=9.5)
        ax2 = axes[1 + n + c]
        rest = y - background - sum(comps[k] for k in range(n) if k != c)
        line(ax2, rest.real, color="black", lw=0.9, label="data minus the other components and background")
        line(ax2, comps[c].real, color=color, lw=1.1, label=labels[c])
        peaks = component_peaks(f[w], comps[c].real[w], fraction)
        for i, p in enumerate(peaks):
            ax2.axvline(p, color=color, lw=0.6, ls=":")
            ax2.text(p, 0.97 - 0.08 * (i % 2), f"{p:.1f}", transform=ax2.get_xaxis_transform(), ha="center",
                     va="top", fontsize=7.5, color=color)
        ax2.set_title(f"c{c + 1}  what {labels[c]} has to explain", loc="left", fontsize=9.5)
        out_rows.append((labels[c], ratio, peaks))
    line(axes[-1], (y - m).real, color="#7b2f8f", lw=0.8)
    axes[-1].set_title("d  residual", loc="left", fontsize=9.5)
    for k, ax in enumerate(axes):
        ax.axhline(0, color="#9aa4b2", lw=0.5)
        ax.tick_params(labelsize=7.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if k < len(axes) - 1:
            ax.set_ylim(-0.55 * ymax, ymax)
            ax.legend(fontsize=7.5, frameon=False, loc="upper left")
        else:
            ax.set_ylim(-0.25 * ymax, 0.25 * ymax)
    axes[-1].set_xlabel("Frequency [Hz]")
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--band", default="", help="lo,hi (Hz); default: the whole fit range")
    ap.add_argument("--peaks", type=float, default=0.15, help="marked component peaks: height >= this x its maximum")
    args = ap.parse_args()
    band = tuple(float(v) for v in args.band.split(",")) if args.band else None
    for label, ratio, peaks in plot(args.run, args.out, band, args.peaks):
        print(f"{label}: amplitude {ratio:.3f} of the largest; peaks {', '.join(f'{p:.2f}' for p in peaks)} Hz")
    print(args.out)


if __name__ == "__main__":
    main()
