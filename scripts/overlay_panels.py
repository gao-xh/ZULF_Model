"""Stacked overlay figure of several fitted data sets, one row per data set: experiment and the fit's own model
(real part), the molecule, one small drawing per 13C isotopologue with its local spin topology (XA = CH,
XA2 = CH2, XA3 = CH3) and coloured bars over the bands its model lines occupy.

    python scripts/overlay_panels.py --group amines [--out runs/figures/overlay_amines.png]
    python scripts/overlay_panels.py --group metabolites | pyridines

Fits: the latest pass of each batch entry (runs/processed/batch_<id>[_v3|_v2]_s2, else _s1). Display baselines:
experiment and simulation each get their own AsLS baseline (owner 2026-10-10), chosen per curve (never one setting
for all) from a small grid by the least negative area near the lines plus the loss of peak height; both parameter
sets are written on each row. Fits are conditional numerical results of the given structures.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

GROUPS = {
    "amines": (["tert-butylamine", "isopropylamine", "triethylamine", "ethylenediamine", "diethylamine",
                "n-n-dimethylethylenediamine", "sec-butylamine", "dipropylamine"], (95.0, 285.0)),
    "metabolites": (["lactic-acid", "alanine", "glycerol", "serine", "methyl-2-bromopropionate"], (95.0, 310.0)),
    "pyridines": (["pyridine-02-mol", "pyridine-10-mol", "pyridine-33-mol", "pyridine-50-mol", "pyridine-66-mol",
                   "pyridine-75-mol", "pyridine-100-mol", "pyridinium-6m-acidic"], (145.0, 195.0)),
}
TOPOLOGY = {1: ("XA", "#2f8f5b", "J"), 2: ("XA2", "#3b6ea8", "3/2 J"), 3: ("XA3", "#b4469b", "J & 2J")}
# drawings for specifications without a recorded molecule (motifs, hand-written chains of the protonated
# carbons only): full SMILES and the site -> atom index map, for the figure only
DRAWINGS = {
    "lactic acid": ("CC(O)C(=O)O", {"C2": 0, "C1": 1}),
    "L-alanine": ("CC([NH3+])C(=O)[O-]", {"C2": 0, "C1": 1}),
    "triethylamine": ("CCN(CC)CC", {"C2": 0, "C1": 1}),
    "ethylenediamine": ("NCCN", {"C1": 1}),
    "pyridine": ("c1ccncc1", {"A4": 0, "A3": 1, "A2": 2}),
}
CONFIGS = ["configs/batch_2026-10-10_v3.json", "configs/batch_2026-10-10_v2.json", "configs/batch_2026-10-10.json"]


# the pass chosen per data set after comparing the passes (docs/analysis/2026-10-10_overnight_plan.md); others: the
# latest of v3 / v2 / v1
# a value may also be a run directory name (a refinement outside run_batch) with the config pass of its structure
RUNS = {"methyl-2-bromopropionate": ("batch_methyl-2-bromopropionate_v7_longrange", "v6")}
CHOSEN = {"methyl-2-bromopropionate": "_v6", "serine": "_v11", "n-n-dimethylethylenediamine": "_v13",
          "diethylamine": "_v13", "sec-butylamine": "_v13",
          "lactic-acid": "_v9", "alanine": "_v9", "glycerol": "_v12"}   # v9: window 1.5 1/s; v11, v12: background off (2026-10-11)


def latest_fit(eid):
    """(run directory, tag) of the chosen pass (CHOSEN) or the latest finished pass of a batch entry."""
    if eid in RUNS and (ROOT / "runs/processed" / RUNS[eid][0] / "fit.json").exists():
        return ROOT / "runs/processed" / RUNS[eid][0], RUNS[eid][1] + " refined"
    order = [CHOSEN[eid]] if eid in CHOSEN else ["_v3", "_v2", ""]
    for tag in order:
        for stage in ("s2", "s1"):
            run = ROOT / "runs/processed" / f"batch_{eid}{tag}_{stage}"
            if (run / "fit.json").exists():
                return run, (tag.strip("_") or "v1") + " " + stage
    return None, None


def fit_series(run):
    """The series file a fit was made on (its command line), so a pass with another window is drawn on its own data."""
    argv = json.loads((Path(run) / "monitor" / "status.json").read_text())["argv"]
    return argv[argv.index("--series") + 1]


def entry_spec(eid, tag):
    """Structure specification of the pass used: the config of that pass first, then the earlier configs."""
    version = tag.split()[0]
    own = "configs/batch_2026-10-10.json" if version == "v1" else f"configs/batch_2026-10-10_{version}.json"
    order = [own] + [c for c in CONFIGS if c != own]
    for cfg in order:
        for e in json.load(open(ROOT / cfg))["entries"]:
            if e["id"] == eid:
                return e["structure"], e
    raise KeyError(eid)


def auto_baseline(f, y, lines_hz):
    """AsLS baseline chosen for this spectrum: grid of (pre-smoothing, smoothness, p); score = negative area within
    6 Hz of the model lines (relative to the absolute area there) + loss of the five highest peaks. Returns
    (baseline, settings)."""
    from scipy.ndimage import gaussian_filter1d
    from zulf_processing import asls_baseline
    step = float(np.median(np.diff(f)))
    near = np.zeros(len(f), bool)
    for v in lines_hz:
        near |= np.abs(f - v) < 6.0
    if not near.any():
        near[:] = True
    raw_top = np.sort(y[near])[-5:].mean()
    best = None
    for pre in (0.0, 0.15):
        ys = gaussian_filter1d(y, pre / step) if pre else y
        for smooth in (0.25, 0.35, 0.5, 0.7):
            for p in (0.05, 0.15):
                b = asls_baseline(ys, f, smooth_hz=smooth, p=p)
                q = (f > f.max() - 30)
                b = b + np.median((y - b)[q])
                dc = y - b
                neg = float(-dc[near][dc[near] < 0].sum() / (np.abs(dc[near]).sum() or 1.0))
                loss = max(0.0, 1.0 - np.sort(dc[near])[-5:].mean() / (raw_top - np.median(y) or 1.0))
                score = neg + loss
                if best is None or score < best[0]:
                    best = (score, b, dict(pre=pre, smooth=smooth, p=p))
    return best[1], best[2]


def draw_molecule(spec, highlight=None, color=None, size=(260, 200)):
    """RGBA array of the molecule (RDKit); highlight = site label drawn as 13C with its protons, in colour."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    from zulf_hypothesis.molecule import molecule_from_structure
    if spec.get("compound") in DRAWINGS and not spec.get("molecule"):
        smiles, sites = DRAWINGS[spec["compound"]]
        mol = Chem.MolFromSmiles(smiles)
    else:
        mol, sites = molecule_from_structure(spec)
    mol = Chem.Mol(mol)
    hl_atoms, colors = [], {}
    if highlight is not None and highlight in sites:
        i = sites[highlight]
        mol.GetAtomWithIdx(i).SetIsotope(13)
        mol = Chem.AddHs(mol, onlyOnAtoms=(i,))
        rgb = tuple(int(color.lstrip("#")[k:k + 2], 16) / 255 for k in (0, 2, 4))
        hl_atoms = [i] + [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetSymbol() == "H"]
        colors = {a: rgb + (0.45,) for a in hl_atoms}
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = d.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.08
    opts.isotopeLabels = True
    d.DrawMolecule(mol, highlightAtoms=hl_atoms, highlightAtomColors=colors, highlightBonds=[])
    d.FinishDrawing()
    import matplotlib.image as mpimg
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")


def proton_counts(spec):
    from zulf_hypothesis.structure_spec import fragment_from_spec
    frag = fragment_from_spec(spec)
    n = {}
    for g in frag.protons:
        n[g.site] = n.get(g.site, 0) + g.size
    return n


def panel(ax, eid, view):
    import matplotlib.pyplot as plt  # noqa: F401
    from matplotlib.offsetbox import AnnotationBbox, OffsetImage
    from zulf_studio.session import StudioSession
    run, tag = latest_fit(eid)
    if run is None:
        ax.text(0.5, 0.5, f"{eid}: no fit", transform=ax.transAxes, ha="center")
        return
    spec, entry = entry_spec(eid, tag)
    s = StudioSession(log_file=False)
    s.load_spectrum(series=str(ROOT / fit_series(run)))
    s.apply_fit(str(run))
    s.set_view(*view)
    s.fit_model_curve()
    sim = s.simulate()
    f = np.asarray(sim["f"])
    y = np.asarray(sim["data_re"])
    m = np.nan_to_num(np.asarray(sim.get("fit_re", np.zeros_like(y))))
    lines = [r for r in s.lines(0.05) if view[0] <= r["frequency_hz"] <= view[1]]
    # experiment and simulation each get their own AsLS baseline, chosen separately (owner, 2026-10-10)
    b, settings = auto_baseline(f, y, [r["frequency_hz"] for r in lines])
    bm, settings_m = auto_baseline(f, m, [r["frequency_hz"] for r in lines])
    dc, mc = y - b, m - bm
    sc = float(np.max(dc)) or 1.0
    ax.plot(f, dc / sc, color="#1d1d1f", lw=0.7, zorder=2)
    ax.plot(f, mc / sc, color="#d4682e", lw=1.8, alpha=0.85, zorder=3)        # the fit drawn over the experiment
    ax.set_xlim(*view)
    ax.set_ylim(-0.15, 1.55)
    ax.set_yticks([0, 0.5, 1.0])
    ax.axhline(0, color="0.88", lw=0.5, zorder=0)
    # bands of each 13C isotopologue: clusters of its model lines (gap > 12 Hz splits J and 2J of CH3)
    counts = proton_counts(spec)
    placed = []
    for comp in sorted({r["component"] for r in lines}):
        mm = re.match(r"13C@([A-Za-z]+\d*)", comp)
        if not mm:
            continue
        site = mm.group(1)
        nh = counts.get(site, 0)
        if nh not in TOPOLOGY:
            continue
        name, color, _ = TOPOLOGY[nh]
        fr = np.sort([r["frequency_hz"] for r in lines if r["component"] == comp])
        clusters = np.split(fr, np.flatnonzero(np.diff(fr) > 12.0) + 1)
        for cl in clusters:
            lo, hi = cl.min() - 1.0, cl.max() + 1.0
            ax.add_patch(plt.Rectangle((lo, 1.42), hi - lo, 0.09, color=color, alpha=0.35, lw=0, zorder=1))
        main = clusters[-1]                    # CH3: its 2J band; CH2 / CH: their only band
        x0 = float(np.mean(main))
        if all(abs(x0 - p) > 0.09 * (view[1] - view[0]) for p in placed):
            img = draw_molecule(spec, site, color, size=(200, 150))
            ab = AnnotationBbox(OffsetImage(img, zoom=0.32), (x0, 1.12), frameon=False, zorder=4)
            ax.add_artist(ab)
            placed.append(x0)
    img = draw_molecule(spec, size=(220, 160))
    ax.add_artist(AnnotationBbox(OffsetImage(img, zoom=0.36), (0.07, 0.66), xycoords="axes fraction",
                                 frameon=True, bboxprops=dict(edgecolor="0.6", lw=0.6), zorder=4))
    fit = json.load(open(run / "fit.json"))
    ax.text(0.07, 0.22, f"{eid.replace('-', ' ')}", transform=ax.transAxes, fontsize=8, fontweight="bold",
            ha="center", va="top")
    # which pass is drawn, in plain sight (owner, 2026-10-11): the fit folder names the version and stage
    ax.text(0.07, 0.13, f"version: {run.name.removeprefix(f'batch_{eid}').strip('_') or 'v1'}", transform=ax.transAxes, fontsize=7.5,
            fontweight="bold", color="#8a3b12", ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#d4682e", lw=0.6))
    def describe(st):
        pre = f"pre {st['pre']} Hz, " if st["pre"] else ""
        return f"{pre}{st['smooth']} Hz, p {st['p']}"
    ax.text(0.995, 0.03, f"fit {tag}, score {fit['scores'][0]:.3f}; AsLS experiment {describe(settings)}; "
            f"simulation {describe(settings_m)}", transform=ax.transAxes, fontsize=6.5, color="0.4", ha="right",
            va="bottom")


def legend(ax):
    import matplotlib.pyplot as plt
    ax.axis("off")
    x = 0.02
    for nh, (name, color, center) in TOPOLOGY.items():
        ax.add_patch(plt.Rectangle((x, 0.35), 0.03, 0.3, color=color, alpha=0.45, transform=ax.transAxes))
        ax.text(x + 0.04, 0.5, f"{name} (13C-H{'' if nh == 1 else nh}): multiplet centre {center}",
                transform=ax.transAxes, va="center", fontsize=9)
        x += 0.27
    ax.plot([0.83, 0.86], [0.62, 0.62], color="#1d1d1f", lw=1, transform=ax.transAxes)
    ax.text(0.865, 0.62, "experiment", transform=ax.transAxes, va="center", fontsize=8)
    ax.plot([0.83, 0.86], [0.32, 0.32], color="#d4682e", lw=1.8, alpha=0.85, transform=ax.transAxes)
    ax.text(0.865, 0.32, "fit (own model)", transform=ax.transAxes, va="center", fontsize=8)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", required=True, choices=sorted(GROUPS))
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    ids, view = GROUPS[args.group]
    fig, axs = plt.subplots(len(ids) + 1, 1, figsize=(12, 1.9 * len(ids) + 0.6),
                            gridspec_kw={"height_ratios": [0.3] + [1] * len(ids)}, sharex=False)
    legend(axs[0])
    for ax, eid in zip(axs[1:], ids):
        panel(ax, eid, view)
        ax.tick_params(labelsize=7)
        if ax is not axs[-1]:
            ax.set_xticklabels([])
    axs[-1].set_xlabel("frequency (Hz)")
    fig.text(0.003, 0.5, "signal (real part, row max = 1)", rotation=90, va="center", fontsize=9)
    fig.tight_layout(rect=(0.012, 0, 1, 1), h_pad=0.25)
    out = Path(args.out or ROOT / f"runs/figures/overlay_{args.group}.png")
    fig.savefig(out, dpi=150)
    print("figure", out)


if __name__ == "__main__":
    main()
