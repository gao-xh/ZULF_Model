"""Check the automatic display baseline on several fits at once (overlay_panels.display_baseline).

    python scripts/baseline_check.py OUT.png RUN:lo,hi [RUN:lo,hi ...]

Per fit two rows: the experiment as fitted with the chosen baseline (blue) and the masked line regions (shaded),
then experiment and fit after the baseline. Printed per fit: the settings, the route, and three numbers that a
wrong baseline shows: the bend (rms distance of the baseline from the straight chord across each masked run),
the largest baseline excursion inside the masked runs and the negative area near the lines, all relative to the
largest data value. A baseline is suspect when the bend or the excursion exceeds about 1 % of the signal.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


def metrics(r):
    f, y, b = r["f"], r["data"], r["base"]
    top = float(np.abs(y).max()) or 1.0
    import overlay_panels as op
    mask = op.line_mask(f, r["lines"], r["widths"])
    edges = np.flatnonzero(np.diff(np.r_[0, mask.astype(int), 0]))
    bend, excursion = [], 0.0
    for a, c in zip(edges[0::2], edges[1::2]):
        if c - a < 3 or (a == 0 and c >= len(f)):
            continue
        lo_i, hi_i = max(a - 1, 0), min(c, len(f) - 1)
        if a == 0 or c >= len(f):          # a run at the view edge: compare with the level at its open side
            chord = np.full(hi_i - lo_i + 1, b[hi_i] if a == 0 else b[lo_i])
        else:
            chord = np.interp(f[lo_i:hi_i + 1], [f[lo_i], f[hi_i]], [b[lo_i], b[hi_i]])
        dev = b[lo_i:hi_i + 1] - chord
        if (f[hi_i] - f[lo_i]) > 15.0:     # long masked stretches may curve legitimately: report, do not score
            continue
        bend.append(np.sqrt(np.mean(dev ** 2)))
        excursion = max(excursion, float(np.abs(dev).max()))
    near = np.zeros(len(f), bool)
    for v in r["lines"]:
        near |= np.abs(f - v) < 6.0
    dc = (y - b)[near]
    neg = float(-dc[dc < 0].sum() / (np.abs(dc).sum() or 1.0))
    return dict(bend=float(np.mean(bend) / top) if bend else 0.0, excursion=excursion / top, negative=neg), mask


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import overlay_panels as op
    out, items = sys.argv[1], sys.argv[2:]
    fig, axes = plt.subplots(2 * len(items), 1, figsize=(15, 3.2 * len(items)))
    for i, item in enumerate(items):
        run, view = item.rsplit(":", 1)
        lo, hi = (float(v) for v in view.split(","))
        r = op.display_baseline(run, (lo, hi))
        mt, mask = metrics(r)
        f, y, m, b, bm = r["f"], r["data"], r["fit"], r["base"], r["base_fit"]
        top = np.abs(y).max()
        flag = "  SUSPECT" if mt["bend"] > 0.01 or mt["excursion"] > 0.01 else ""
        text = (f"{Path(run).name}: {op.describe_baseline(r['settings'])}; fit {r['route']}; bend {100 * mt['bend']:.2f} %, "
                f"excursion {100 * mt['excursion']:.2f} %, negative {mt['negative']:.3f}{flag}")
        print(text)
        a = axes[2 * i]
        a.plot(f, y / top, color="black", lw=0.6)
        a.plot(f, b / top, color="#3d7cc9", lw=1.2)
        a.fill_between(f, -1, 1, where=mask, color="#24476e", alpha=0.07, lw=0)
        a.set_ylim(-0.3, 0.3)
        a.set_title(text, loc="left", fontsize=9, color="#c8553d" if flag else "#24476e")
        c = axes[2 * i + 1]
        c.plot(f, (y - b) / top, color="black", lw=0.7)
        c.plot(f, (m - bm) / top, color="#d4682e", lw=1.0)
        c.axhline(0, color="0.8", lw=0.5)
        for ax in (a, c):
            ax.set_xlim(lo, hi)
            for side in ("top", "right"):
                ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=110)


if __name__ == "__main__":
    main()
