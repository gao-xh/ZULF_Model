"""Diagnostic figures of the finished run_batch entries -> runs/figures/batch/<id><tag>.png

    python scripts/batch_figures.py [--tag=_v2] [--redo]

The figures show the finished batch entries: experiment and the stage-2 fit's own model (real part, no
display baseline), whole grid plus the strongest band, residual offset below. Skips figures that exist."""
import json, sys
from pathlib import Path
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import os
os.chdir(ROOT)
from zulf_studio.session import StudioSession
TAG = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--tag=")), "")
entries = json.load(open("configs/batch_2026-10-10.json"))["entries"]
for e in entries:
    eid = e["id"]; fit = Path(f"runs/processed/batch_{eid}{TAG}_s2"); out = Path(f"runs/figures/batch/{eid}{TAG}.png")
    if not (fit / "fit.json").exists() or (out.exists() and "--redo" not in sys.argv):
        continue
    s = StudioSession(log_file=False); s.load_spectrum(series=f"runs/series/batch_{eid}/series.json"); s.apply_fit(str(fit))
    lo, hi = (float(v) for v in (e.get("grid") or "20,380").split(","))
    s.set_view(lo, hi)
    s.fit_model_curve()                                    # the fit's own model (else simulate() has only a quick look)
    sim = s.simulate()
    if "fit_re" not in sim:
        print("no exact model for", eid, flush=True)
        continue
    f = np.asarray(sim["f"]); d = np.asarray(sim["data_re"]); m = np.asarray(sim["fit_re"])
    k = (f >= lo) & (f <= hi)
    ranges = json.load(open(f"runs/series/batch_{eid}/series.json"))[0]["ranges"]
    inside = np.zeros(len(f), bool)
    for a, b in ranges:                                   # fit ranges only (no mains or instrument lines)
        inside |= (f >= a) & (f <= b)
    kk = k & inside & (f > 90.0)                         # above the low-frequency drift
    c = f[kk][np.argmax(np.abs(d[kk]))]                   # strongest fitted point: zoom +-15 Hz around it
    score = json.load(open(fit / "fit.json"))["scores"][0]
    fig, axs = plt.subplots(2, 1, figsize=(11, 7))
    for ax, (a, b) in zip(axs, [(lo, hi), (c - 15, c + 15)]):
        q = (f >= a) & (f <= b); sc = np.abs(d[q]).max()
        ax.plot(f[q], d[q] / sc, "k", lw=0.7, label="experiment")
        ax.plot(f[q], m[q] / sc, color="#d4682e", lw=1.1, label="fit (stage 2)")
        ax.plot(f[q], (d[q] - m[q]) / sc - 0.6, color="0.5", lw=0.6, label="residual (offset)")
        ax.set_xlim(a, b); ax.set_xlabel("frequency (Hz)")
    axs[0].set_title(f"{eid}{TAG}: stage 2 score {score:.4f}  (conditional fit of the given structure)", loc="left", fontsize=10)
    axs[0].legend(frameon=False, fontsize=8, ncol=3)
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)
    print("figure", out, flush=True)
