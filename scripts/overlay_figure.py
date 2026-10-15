"""One overlay panel of a fit for slides: experiment and simulation on one axis, nothing else.

    python scripts/overlay_figure.py <the fit_joint_series options of the fit> --fit RUN/fit.json \\
        --fid DATA/average_fid.npy [--view 100,275] [--display-window 0.3] [--size 8.6,3.4] \\
        [--sim-label Simulation] --figure OUT.png [--formats png,svg,pdf]

The display processing, the model rendering (gains of the fit ranges) and the display baseline are those of
scripts/paper_figure.py; both curves are scaled by the largest experimental value in the range.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import fit_joint_series as fj                                            # noqa: E402
from j_tuner import load_fit                                             # noqa: E402
from paper_figure import display_baseline, display_spectrum, render_model  # noqa: E402

EXPERIMENT, SIMULATION = "#1d1d1f", "#d4682e"


def main():
    import json
    ap = fj.make_parser()
    ap.add_argument("--fit", required=True)
    ap.add_argument("--fid", required=True)
    ap.add_argument("--view", default="100,275", help="frequency range shown (Hz)")
    ap.add_argument("--display-crop", type=float, default=0.1)
    ap.add_argument("--display-window", type=float, default=0.3)
    ap.add_argument("--display-zero-fill", type=int, default=4)
    ap.add_argument("--size", default="8.6,3.4", help="figure width,height (inches)")
    ap.add_argument("--sim-label", default="Simulation")
    ap.add_argument("--sim-width", type=float, default=1.9, help="line width of the simulation")
    ap.add_argument("--figure", required=True)
    ap.add_argument("--formats", default="png,svg,pdf")
    ap.add_argument("--baseline", default="model", choices=("model", "shared", "separate", "residual", "none"),
                    help="display baseline (zulf_processing.display_baseline): separate (data and model each "
                         "their own), shared (one from the data, on both), residual (from data - model, on the "
                         "data), none")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prob = fj.build_problem(args)
    load_fit(prob, json.load(open(args.fit)))
    entry = prob.series[0]
    fs = float((entry.get("record") or {}).get("sampling_rate_hz", 4000.0))
    phase0 = (entry.get("phasing") or {}).get("phase0_rad")
    lo, hi = (float(v) for v in args.view.split(","))
    f, spec, acq, phasing, _ = display_spectrum(args.fid, args.display_crop, args.display_window,
                                                args.display_zero_fill, fs=fs, phase0_rad=phase0)
    fit_ranges = [tuple(r) for r in entry.get("ranges", [(prob.lo, prob.hi)])]
    f, y, m, gains, lines = render_model(prob, prob.z0, f, spec, acq, phasing, (lo - 5, hi + 5), fit_ranges)
    phi = float(np.angle(gains[0]))
    phi -= np.pi * round(phi / np.pi)
    y, m = (y * np.exp(-1j * phi)).real, (m * np.exp(-1j * phi)).real
    yc, mc, _ = display_baseline(f, y, m, lines, 1.0, 1.5, 2.5, method=args.baseline)
    sel = (f >= lo) & (f <= hi)
    scale = 1000.0 / float(np.max(yc[sel]))
    w, h = (float(v) for v in args.size.split(","))
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 12})
    fig, ax = plt.subplots(figsize=(w, h), layout="constrained")
    ax.plot(f[sel], yc[sel] * scale, color=EXPERIMENT, lw=0.8, label="Experiment")
    ax.plot(f[sel], mc[sel] * scale, color=SIMULATION, lw=args.sim_width, label=args.sim_label)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("Frequency [Hz]", fontsize=14)
    ax.set_title("Signal [a.u.]", loc="left", fontsize=13, pad=6)
    ax.legend(loc="upper right", frameon=False, ncol=2, fontsize=12, bbox_to_anchor=(1.0, 1.13))
    ax.tick_params(labelsize=12)
    out = Path(args.figure)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in args.formats.split(","):
        fig.savefig(out.with_suffix("." + ext.strip()), dpi=args.dpi)
    print(out.with_suffix(".png"))


if __name__ == "__main__":
    main()
