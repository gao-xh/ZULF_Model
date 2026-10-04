"""Static pictures of a fit trace (fit_joint_series --trace): a contact sheet of a few frames and an animated GIF
of all frames (data, model, residual and the objective along the path).

    python scripts/trace_frames.py runs/processed/RUN [--band 112,140] [--spectrum 0] [--panels 6]
"""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--band", default="")
    ap.add_argument("--spectrum", type=int, default=0)
    ap.add_argument("--panels", type=int, default=6)
    args = ap.parse_args()
    run = Path(args.run)
    meta = json.loads((run / "trace.json").read_text())
    arr = np.load(run / "trace.npz")
    s = args.spectrum
    f, y, models = arr[f"f{s}"], arr[f"y{s}"], arr[f"model{s}"]
    lo, hi = (float(v) for v in args.band.split(",")) if args.band else (float(f.min()), float(f.max()))
    sel = (f >= lo) & (f <= hi)
    scale = np.abs(y[sel]).max()
    obj = np.array([fr["objective"] for fr in meta["frames"]])
    keys = list(meta["frames"][0]["J"])

    def frame(ax, k, compact=False):
        m = models[k]
        ax.plot(f[sel], y.real[sel] / scale, "k", lw=0.8, label="data")
        ax.plot(f[sel], m.real[sel] / scale, c="#c0392b", lw=1.0, label="model")
        ax.plot(f[sel], (y.real[sel] - m.real[sel]) / scale - 0.9, c="0.6", lw=0.7, label="residual (offset)")
        fr = meta["frames"][k]
        ax.set_title(f"frame {k + 1}/{len(meta['frames'])}, evaluation {fr['evaluation']}, {fr['stage']}, "
                     f"objective {fr['objective']:.4g}", fontsize=8 if compact else 9, loc="left")
        ax.set_yticks([])
        ax.set_xlim(lo, hi)

    picks = sorted(set(np.linspace(0, len(models) - 1, args.panels).round().astype(int).tolist()))
    fig, axes = plt.subplots(len(picks), 1, figsize=(11, 2.0 * len(picks)), sharex=True)
    for ax, k in zip(np.atleast_1d(axes), picks):
        frame(ax, k, compact=True)
    np.atleast_1d(axes)[0].legend(fontsize=7, frameon=False, ncol=3, loc="upper right")
    np.atleast_1d(axes)[-1].set_xlabel("frequency (Hz)")
    fig.suptitle(f"Fit trace ({meta['origin']}, {meta['evaluations']} evaluations), real part", fontsize=10)
    fig.tight_layout()
    fig.savefig(run / "trace_sheet.png", dpi=110)
    plt.close(fig)
    images = []
    for k in range(len(models)):
        fig = plt.figure(figsize=(10, 5.2))
        ax = fig.add_axes([0.05, 0.42, 0.93, 0.5])
        frame(ax, k)
        ax.legend(fontsize=7, frameon=False, ncol=3, loc="upper right")
        bx = fig.add_axes([0.05, 0.08, 0.42, 0.24])
        bx.semilogy(obj, c="#2563eb", lw=1)
        bx.plot(k, obj[k], "o", c="#c0392b")
        bx.set_title("objective along the path", fontsize=8, loc="left")
        bx.tick_params(labelsize=7)
        tx = fig.add_axes([0.55, 0.05, 0.43, 0.3])
        tx.axis("off")
        rows = [f"{n}: {meta['frames'][k]['J'][n][0]:8.3f}  (start {meta['frames'][0]['J'][n][0]:.2f}, "
                f"final {meta['frames'][-1]['J'][n][0]:.2f})" for n in keys]
        tx.text(0, 1, "\n".join(rows), va="top", family="monospace", fontsize=7)
        fig.canvas.draw()
        images.append(Image.frombuffer("RGBA", fig.canvas.get_width_height(), fig.canvas.buffer_rgba()).convert("P"))
        plt.close(fig)
    images[0].save(run / "trace.gif", save_all=True, append_images=images[1:], duration=180, loop=0)
    print(run / "trace_sheet.png", run / "trace.gif")


if __name__ == "__main__":
    main()
