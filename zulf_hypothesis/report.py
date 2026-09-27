"""Automatic report of a search or structure fit: ranked table, J matrix of the leading fits, phased fit figure.

`write_report(result, observed, prefix)` writes
  prefix.json  - table, log, J matrices, phasing of the best fit;
  prefix.md    - ranked table and J matrices (Markdown);
  prefix.png   - the best fit phase-corrected with its own fitted phase (needs matplotlib; skipped otherwise).
Couplings built as 0 Hz and held there (unspecified in the structure) are marked as fixed.
"""
from __future__ import annotations

import json
from typing import Dict, List, Optional, Tuple

import numpy as np

from zulf_core.render.phasing import phase_correct

from .search import Evaluated, SearchResult


def j_matrix(e: Evaluated) -> dict:
    """Labels (labelled sites first, then proton groups, in fragment order) and the refined coupling matrix of
    one fit; None where the model has no such coupling; `fixed` marks couplings held at their built value."""
    model = e.model
    keys = list(model.coupling_names)
    members = []
    for k in keys:
        a, b = k[2:-1].split(",")
        members += [a, b]
    order = []
    if model.fragment is not None:
        order = [s.label for s in model.fragment.sites if s.label_isotopes()] + [p.label for p in model.fragment.protons]
    labels = [l for l in order if l in members] + sorted(set(members) - set(order))
    index = {l: i for i, l in enumerate(labels)}
    n = len(labels)
    values: List[List[Optional[float]]] = [[None] * n for _ in range(n)]
    fixed: List[List[bool]] = [[False] * n for _ in range(n)]
    held = set(model.unspecified)
    for k in keys:
        a, b = k[2:-1].split(",")
        if k not in e.couplings:
            continue
        i, j = index[a], index[b]
        values[i][j] = values[j][i] = round(float(e.couplings[k]), 3)
        fixed[i][j] = fixed[j][i] = k in held
    return {"hypothesis": e.key, "labels": labels, "J_hz": values, "fixed": fixed}


def j_matrix_markdown(m: dict) -> str:
    labels = m["labels"]
    lines = ["| | " + " | ".join(labels) + " |", "|---" * (len(labels) + 1) + "|"]
    for i, a in enumerate(labels):
        cells = []
        for j in range(len(labels)):
            v = m["J_hz"][i][j]
            cells.append("" if v is None else (f"{v:g} (fixed)" if m["fixed"][i][j] else f"{v:g}"))
        lines.append(f"| {a} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def fit_phasing(e: Evaluated) -> Tuple[float, float]:
    """(phase0_rad, delay_s) of a fit: the shared gain phase (strongest component) and the fitted phase delay.
    phase_correct(values, f, phase0, delay, acquisition) then gives the absorption spectrum."""
    gains = np.array([complex(*g) for g in e.summary["gains"]])
    phase0 = float(np.angle(gains[np.argmax(np.abs(gains))])) if len(gains) else 0.0
    return phase0, float(e.summary["parameters"].get("phase_delay", 0.0))


def _zoom_windows(f: np.ndarray, mask: np.ndarray, gap_hz: float = 6.0, pad_hz: float = 3.0,
                  limit: int = 4) -> List[Tuple[float, float]]:
    pts = np.sort(f[mask])
    if not len(pts):
        return []
    groups, start, prev = [], pts[0], pts[0]
    for x in pts[1:]:
        if x - prev > gap_hz:
            groups.append((start, prev))
            start = x
        prev = x
    groups.append((start, prev))
    groups = sorted(groups, key=lambda g: g[1] - g[0], reverse=True)[:limit]
    return sorted((lo - pad_hz, hi + pad_hz) for lo, hi in groups)


def plot_fit(result: SearchResult, observed, path: str, which: Optional[Evaluated] = None) -> Optional[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    e = which or result.best
    sel = observed.selected
    f, y, p = observed.frequencies_hz[sel], observed.values[sel], np.asarray(e.prediction)
    band = observed.band_index[sel]
    if observed.real_only:
        yd, pd_ = y.real, p.real
        phase_note = "phased input"
    else:
        phase0, delay = fit_phasing(e)
        ph = phase_correct(np.ones(len(f), complex), f, phase0, delay, observed.acquisition)
        yd, pd_ = (y * ph).real, (p * ph).real
        phase_note = f"phase0 {np.degrees(phase0):.0f} deg, delay {delay * 1e3:.2f} ms (from the fit)"
    s = result.yardstick.sigma
    windows = _zoom_windows(f, result.yardstick.mask)
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, max(len(windows), 1), height_ratios=[1, 1.2])
    ax = fig.add_subplot(gs[0, :])
    for b in np.unique(band):
        m = band == b
        ax.plot(f[m], yd[m] / s, color="#222222", lw=0.8, label="data" if b == band.min() else None)
        ax.plot(f[m], pd_[m] / s, color="#d1495b", lw=0.9, alpha=0.9, label="fit" if b == band.min() else None)
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.set_title("phase-corrected real part / noise sigma")
    ax.set_xlabel("frequency (Hz)")
    for k, (lo, hi) in enumerate(windows):
        ax = fig.add_subplot(gs[1, k])
        m = (f > lo) & (f < hi)
        off = 1.15 * max(np.abs(yd[m]).max() / s, 1.0)
        ax.plot(f[m], yd[m] / s, color="#222222", lw=1.0, label="data")
        ax.plot(f[m], pd_[m] / s, color="#d1495b", lw=1.0, label="fit")
        ax.plot(f[m], (yd[m] - pd_[m]) / s - off, color="#00798c", lw=0.9, label="residual (offset)")
        ax.set_title(f"{lo:.0f}-{hi:.0f} Hz")
        ax.set_xlabel("frequency (Hz)")
        if k == 0:
            ax.legend(frameon=False, fontsize=8)
    red = e.chi2 / max(e.n - e.k, 1)
    fig.suptitle(f"{e.key}: reduced chi2 {red:.1f} on the data cores, k {e.k}; {phase_note}\n"
                 f"findings: {', '.join(f.code for f in e.findings) or 'none'}", x=0.01, ha="left",
                 fontweight="bold", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


def write_report(result: SearchResult, observed, prefix: str, top: int = 3) -> Dict[str, str]:
    table = result.table()
    ranked = result.ranked()
    leaders = [result.best] + [e for e in ranked if e is not result.best][:top - 1]
    matrices = [j_matrix(e) for e in leaders if e is not None]
    best = result.best
    phase0, delay = fit_phasing(best) if best is not None else (0.0, 0.0)
    out = {"best": best.key if best else None, "c_hat": result.c_hat, "table": table, "log": result.log,
           "j_matrices": matrices, "best_phasing": {"phase0_rad": phase0, "delay_s": delay},
           "best_rates_per_s": {k.replace("log_rate", "rate"): float(np.exp(v))
                                for k, v in (best.summary["parameters"].items() if best else []) if "log_rate" in k},
           "best_gaussian_sigma_hz": {k: float(v) for k, v in (best.summary["parameters"].items() if best else [])
                                      if k.endswith(".sigma")},
           "best_boundary_hits": best.summary.get("boundary_hits") if best else None}
    paths = {"json": prefix + ".json", "md": prefix + ".md"}
    with open(paths["json"], "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    md = [f"# {best.key if best else 'no result'}", "",
          f"c_hat {result.c_hat:.2f}; best phase0 {np.degrees(phase0):.1f} deg, delay {delay * 1e3:.3f} ms; "
          f"rates (1/s) {', '.join(f'{v:.2f}' for v in out['best_rates_per_s'].values())}; Gaussian sigma (Hz) "
          f"{', '.join(f'{v:.3f}' for v in out['best_gaussian_sigma_hz'].values()) or 'none'}; parameters at bounds "
          f"{out['best_boundary_hits'] or 'none'}", "",
          "| rank | hypothesis | variant | status | delta | chi2 | k | findings |", "|---|---|---|---|---|---|---|---|"]
    for r in table:
        md.append(f"| {r['rank']} | {r['hypothesis']} | {r['variant']} | {r['status']} | {r['delta']} | {r['chi2']} | "
                  f"{r['k']} | {', '.join(r['findings'])}{' (demoted)' if r.get('demoted') else ''} |")
    for m in matrices:
        md += ["", f"## J matrix (Hz): {m['hypothesis']}", "", j_matrix_markdown(m)]
    with open(paths["md"], "w") as fh:
        fh.write("\n".join(md) + "\n")
    png = plot_fit(result, observed, prefix + ".png") if best is not None else None
    if png:
        paths["png"] = png
    return paths
