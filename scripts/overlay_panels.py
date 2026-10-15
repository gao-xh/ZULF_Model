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
    "amines": (["tert-butylamine", "isopropylamine", "triethylamine", "ethylenediamine",
                "n-n-dimethylethylenediamine", "sec-butylamine", "dipropylamine"], (95.0, 285.0)),
    "metabolites": (["lactic-acid", "alanine", "glycerol", "serine", "methyl-2-bromopropionate"], (95.0, 310.0)),
    "pyridines": (["pyridine-02-mol", "pyridine-10-mol", "pyridine-33-mol", "pyridine-50-mol", "pyridine-66-mol",
                   "pyridine-75-mol", "pyridine-100-mol"], (145.0, 195.0)),     # pyridinium left out (owner, 2026-10-12)
}
# left out of the overlays, with the reason (owner's decision)
EXCLUDED = {"diethylamine": "experimental spectrum too poor: magnetometer offset at the ADC rail, half of the scans "
                            "pinned, weak signal (owner, 2026-10-12; docs/analysis/2026-10-11_amines_refits.md)"}
# figure width (inches) per group; default 12 (owner, 2026-10-12: amines narrower)
FIG_WIDTH = {"amines": 8.2}   # panel width : height about 4.5 : 1 (owner)
# display baselines chosen by hand for a spectrum (owner's request), used for experiment and simulation alike:
# eid -> (mask half-width Hz, smoothness Hz) of the line-masked smoother
BASELINE_OVERRIDES = {"lactic-acid": (0.6, 0.7),   # 2026-10-12: humps at 135-145 and 150-155 Hz removed
                      "glycerol": (0.6, 0.7)}      # 2026-10-12: crop-wing humps at 135, 150, 220-224 Hz removed
# groups whose 13C sites share one topology are coloured by site instead (owner, 2026-10-12: pyridine is all
# 13C-H, so the topology legend says nothing); site label -> (legend text, colour), one band row per site
SITE_STYLES = {"pyridines": {"A2": ("13C at C2 (ortho, next to N): its lines", "#2f8f5b"),
                             "A3": ("13C at C3 (meta): its lines", "#3b6ea8"),
                             "A4": ("13C at C4 (para): its lines", "#b4469b")}}
TOPOLOGY = {1: ("XA", "#2f8f5b", "J"), 2: ("XA2", "#3b6ea8", "3/2 J"), 3: ("XA3", "#b4469b", "J & 2J")}
# drawings for specifications without a recorded molecule (motifs, hand-written chains of the protonated
# carbons only): full SMILES and the site -> atom index map, for the figure only
DRAWINGS = {
    "lactic acid": ("CC(O)C(=O)O", {"C2": 0, "C1": 1}),
    "L-alanine": ("CC([NH3+])C(=O)[O-]", {"C2": 0, "C1": 1}),
    "triethylamine": ("CCN(CC)CC", {"C2": 0, "C1": 1}),
    "ethylenediamine": ("NCCN", {"C1": 1}),
    "pyridine": ("c1ccncc1", {"A4": 0, "A3": 1, "A2": 2}),
    "pyridinium": ("c1cc[nH+]cc1", {"A4": 0, "A3": 1, "A2": 2}),
}
CONFIGS = ["configs/batch_2026-10-10_v3.json", "configs/batch_2026-10-10_v2.json", "configs/batch_2026-10-10.json"]


# the pass chosen per data set after comparing the passes (docs/analysis/2026-10-10_overnight_plan.md); others: the
# latest of v3 / v2 / v1
# a value may also be a run directory name (a refinement outside run_batch) with the config pass of its structure
RUNS = {"methyl-2-bromopropionate": ("batch_methyl-2-bromopropionate_v7_longrange", "v6"),
        "n-n-dimethylethylenediamine": ("sign_scan/n-n-dimethylethylenediamine_v13_s2sign", "v13")}   # sign scan D65, 2026-10-15
CHOSEN = {"methyl-2-bromopropionate": "_v6", "serine": "_v11", "n-n-dimethylethylenediamine": "_v13",
          "diethylamine": "_v14", "sec-butylamine": "_v13", "tert-butylamine": "_v14", "ethylenediamine": "_v14",
          "isopropylamine": "_v16", "triethylamine": "_v16",
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


def fit_features(run):
    """Short labels of the non-default model and processing choices of a fit, read from its recorded command, its
    result and its series (owner, 2026-10-12: mark features such as slow N-H exchange on the overlay)."""
    run = Path(run)
    argv = json.loads((run / "monitor" / "status.json").read_text())["argv"]
    def opt(name, default=None):
        return argv[argv.index(name) + 1] if name in argv else default
    out = []
    if opt("--exchange") == "slow":
        out.append("slow N-H exchange" + (" (k fitted)" if opt("--nh-exchange") else ""))
    sp = next(iter(json.loads((run / "fit.json").read_text()).get("spectrum_parameters", {}).values()), {})
    if "--fit-field" in argv or opt("--field"):
        bt, bz = 1e3 * sp.get("field_transverse_ut", 0.0), 1e3 * sp.get("field_z_ut", 0.0)
        out.append(f"residual field {bt:.0f} / {bz:.0f} nT ({'fitted' if '--fit-field' in argv else 'fixed'})")
    fitj = json.loads((run / "fit.json").read_text())
    if (fitj.get("shape_refine") or {}).get("accepted"):
        out.append("shape refine (D67)")
    if (fitj.get("sign_scan") or {}).get("accepted"):
        out.append(f"sign scan: {len(fitj['sign_scan']['accepted'])} flipped")
    elif "sign_scan" in str(run):
        out.append("from a sign scan (D65)")
    if (fitj.get("peak_refine") or {}).get("accepted"):
        out.append("peak refine (D66)")
    if opt("--variant") == "free":
        out.append("free isotopologue amplitudes")
    if opt("--background-order") == "-1":
        out.append("background off")
    if opt("--band-weighting") == "equal":
        out.append("equal weighting")
    record = json.loads((ROOT / opt("--series")).read_text())[0].get("record", {})
    apod = record.get("apodization_rate_per_s")
    if apod is not None and abs(apod - 0.3) > 1e-9:
        out.append(f"window {apod:g} 1/s")
    return out


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


def auto_baseline(f, y, lines_hz, methods=("asls", "masked"), with_score=False):
    """AsLS baseline chosen for this spectrum: grid of (pre-smoothing, smoothness, p); score = negative area within
    6 Hz of the model lines (relative to the absolute area there) + the mean height loss over the line clusters
    (lines closer than 6 Hz form one cluster). Every cluster counts, not only the highest peaks: a weak, broad band
    (a 2J band next to a strong J band) must not be taken into the baseline (2026-10-11: the simulation's 2J band
    of triethylamine lost two thirds of its height when only the five highest peaks were checked). Returns
    (baseline, settings)."""
    from scipy.ndimage import gaussian_filter1d
    from zulf_processing import asls_baseline
    step = float(np.median(np.diff(f)))
    near = np.zeros(len(f), bool)
    for v in lines_hz:
        near |= np.abs(f - v) < 6.0
    if not near.any():
        near[:] = True
    clusters = []
    for v in sorted(lines_hz):
        if clusters and v - clusters[-1][1] < 6.0:
            clusters[-1][1] = v
        else:
            clusters.append([v, v])
    regions = [(f > lo - 1.5) & (f < hi + 1.5) for lo, hi in clusters] or [near]
    regions = [r for r in regions if r.any()]
    level = float(np.median(y))
    raw_height = [float(y[r].max() - level) for r in regions]
    mean_height = float(np.mean([h for h in raw_height if h > 0] or [1.0]))
    # line cores; humps beside them (crop wings, broad residual structure) are what the third term measures
    cores = {w: np.zeros(len(f), bool) for w in (0.6, 0.8, 1.2)}
    for v in lines_hz:
        for w in cores:
            cores[w] |= np.abs(f - v) < w
    beside = near & ~cores[0.8]

    def score_of(b):
        dc = y - b
        neg = float(-dc[near][dc[near] < 0].sum() / (np.abs(dc[near]).sum() or 1.0))
        loss = float(np.mean([max(0.0, 1.0 - dc[r].max() / (h or 1.0)) for r, h in zip(regions, raw_height)
                              if h > 0] or [0.0]))
        hump = float(np.mean(np.abs(dc[beside])) / mean_height) if beside.any() else 0.0
        # the hump term only breaks near-ties: weighted fully it let a flexible AsLS cut the DMEDA main line
        return neg + loss + 0.3 * hump

    best = (np.inf, None, None)
    for pre in ((0.0, 0.15) if "asls" in methods else ()):
        ys = gaussian_filter1d(y, pre / step) if pre else y
        for smooth in (0.25, 0.35, 0.5, 0.7):
            for p in (0.05, 0.15):
                b = asls_baseline(ys, f, smooth_hz=smooth, p=p)
                q = (f > f.max() - 30)
                b = b + np.median((y - b)[q])
                score = score_of(b)
                if score < best[0]:
                    best = (score, b, dict(pre=pre, smooth=smooth, p=p))
    # masked smoother (owner, 2026-10-12: remove the humps beside the lines): the line cores get weight 0, the
    # baseline is a Whittaker smoother through the rest, so it follows the humps but not the lines
    for w, core in (cores.items() if "masked" in methods else ()):
        if core.mean() > 0.8:                          # dense multiplets: nothing left to fit between the lines
            continue
        for smooth in (0.7, 1.0, 2.0, 3.0):
            b = masked_smoother(f, y, ~core, smooth)
            q = (f > f.max() - 30)
            b = b + np.median((y - b)[q])
            score = score_of(b)
            if score < best[0]:
                best = (score, b, dict(pre=0.0, smooth=smooth, p=None, mask=w))
    if best[1] is None:                                # no candidate of the requested kind: plain AsLS
        return auto_baseline(f, y, lines_hz, ("asls",), with_score)
    return (best[1], best[2], best[0]) if with_score else (best[1], best[2])


def masked_smoother(f, y, weight, smooth_hz):
    """Whittaker smoother (second differences, lambda = (smooth / step)^4) with 0 / 1 weights: fits only where
    weight is True and interpolates smoothly across the rest."""
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    n = len(y)
    step = float(np.median(np.diff(f)))
    lam = (smooth_hz / step) ** 4
    d = sparse.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    w = sparse.diags(np.asarray(weight, float))
    return spsolve(sparse.csc_matrix(w + lam * (d.T @ d)), np.asarray(weight, float) * y)


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


def panel(ax, eid, view, site_styles=None):
    k = min(1.0, ax.figure.get_figwidth() / 12)        # narrow figures: smaller drawings and labels
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
    # the kind of baseline (AsLS or lines masked) is chosen for experiment and simulation together, so both are
    # processed alike; the parameters of that kind are chosen for each separately (owner, 2026-10-10 / 10-12)
    hz = [r["frequency_hz"] for r in lines]
    if eid in BASELINE_OVERRIDES:
        w, sm = BASELINE_OVERRIDES[eid]
        core = np.zeros(len(f), bool)
        for v in hz:
            core |= np.abs(f - v) < w
        tail = f > f.max() - 30
        b = masked_smoother(f, y, ~core, sm)
        b = b + np.median((y - b)[tail])
        bm = masked_smoother(f, m, ~core, sm)
        bm = bm + np.median((m - bm)[tail])
        settings = settings_m = dict(pre=0.0, smooth=sm, p=None, mask=w)
    options = []
    for kind in ((("asls",), ("masked",)) if eid not in BASELINE_OVERRIDES else ()):
        b, settings, s_exp = auto_baseline(f, y, hz, kind, with_score=True)
        bm, settings_m, s_sim = auto_baseline(f, m, hz, kind, with_score=True)
        options.append((s_exp + s_sim, b, settings, bm, settings_m))
    if options:
        _, b, settings, bm, settings_m = min(options, key=lambda o: o[0])
    # two routes for the simulation (owner, 2026-10-13: "take both, keep the better"): the experiment's own baseline
    # subtracted from both (common: the drawn difference equals the true residual) or the simulation's own baseline
    # (separate). Keep the one whose drawn experiment - simulation is smaller near the lines.
    near = np.zeros(len(f), bool)
    for v in hz:
        near |= np.abs(f - v) < 6.0
    near = near if near.any() else np.ones(len(f), bool)
    rms_common = float(np.sqrt(np.mean((y - m)[near] ** 2)))
    rms_separate = float(np.sqrt(np.mean(((y - b) - (m - bm))[near] ** 2)))
    route = "common" if rms_common <= rms_separate else "separate"
    if route == "common":
        bm, settings_m = b, dict(settings, common=True)
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
        if site_styles and site in site_styles:
            color = site_styles[site][1]
            row = list(site_styles).index(site)
            y_band, h_band = 1.47 - 0.07 * row, 0.055
        elif nh in TOPOLOGY:
            name, color, _ = TOPOLOGY[nh]
            y_band, h_band = 1.42, 0.09
        else:
            continue
        fr = np.sort([r["frequency_hz"] for r in lines if r["component"] == comp])
        clusters = np.split(fr, np.flatnonzero(np.diff(fr) > 12.0) + 1)
        if site_styles and site in site_styles:
            # one row per site: a tick at each of its lines (relative >= 0.15), shade by strength
            own = [r for r in lines if r["component"] == comp and r["relative"] >= 0.15]
            for r in own:
                ax.plot([r["frequency_hz"]] * 2, [y_band, y_band + h_band], color=color, lw=1.6,
                        alpha=0.35 + 0.65 * min(1.0, r["relative"]), zorder=1, solid_capstyle="butt")
            x0 = max(own, key=lambda r: r["relative"])["frequency_hz"] if own else float(np.mean(fr))
            width = view[1] - view[0]
            while any(abs(x0 - q) < 0.09 * width for q in placed):        # every site gets its drawing
                x0 += 0.09 * width
            if x0 > view[1] - 0.05 * width:
                x0 = view[0] + 0.2 * width
            img = draw_molecule(spec, site, color, size=(200, 150))
            ax.add_artist(AnnotationBbox(OffsetImage(img, zoom=0.32 * k), (x0, 1.1), frameon=False, zorder=4))
            placed.append(x0)
            continue
        for cl in clusters:
            lo, hi = cl.min() - 1.0, cl.max() + 1.0
            ax.add_patch(plt.Rectangle((lo, y_band), hi - lo, h_band, color=color, alpha=0.45, lw=0, zorder=1))
        main = clusters[-1]                    # CH3: its 2J band; CH2 / CH: their only band
        x0 = float(np.mean(main))
        if all(abs(x0 - p) > 0.09 * (view[1] - view[0]) for p in placed):
            img = draw_molecule(spec, site, color, size=(200, 150))
            ab = AnnotationBbox(OffsetImage(img, zoom=0.32 * k), (x0, 1.12), frameon=False, zorder=4)
            ax.add_artist(ab)
            placed.append(x0)
    img = draw_molecule(spec, size=(220, 160))
    ax.add_artist(AnnotationBbox(OffsetImage(img, zoom=0.36 * k), (0.07, 0.66), xycoords="axes fraction",
                                 frameon=True, bboxprops=dict(edgecolor="0.6", lw=0.6), zorder=4))
    fit = json.load(open(run / "fit.json"))
    ax.text(0.07, 0.22, f"{eid.replace('-', ' ')}", transform=ax.transAxes, fontsize=8 * k ** 0.5, fontweight="bold",
            ha="center", va="top")
    # which pass is drawn, in plain sight (owner, 2026-10-11): the fit folder names the version and stage
    ax.text(0.07, 0.13, f"version: {run.name.removeprefix(f'batch_{eid}').removeprefix(eid).strip('_') or 'v1'}", transform=ax.transAxes, fontsize=7.5 * k ** 0.5,
            fontweight="bold", color="#8a3b12", ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#d4682e", lw=0.6))
    features = fit_features(run)
    if features:
        # one feature per line between the molecule box and the name (95-118 Hz: no lines in these spectra)
        ax.text(0.07, 0.405, "\n".join(features), transform=ax.transAxes, fontsize=6.8 * k ** 0.5, color="#1f4e8c",
                ha="center", va="top", fontweight="bold", linespacing=1.15, zorder=5)
    def describe(st):
        if st.get("common"):
            return "same baseline as the experiment"
        if st.get("mask"):
            return f"lines masked +-{st['mask']} Hz, {st['smooth']} Hz"
        pre = f"pre {st['pre']} Hz, " if st["pre"] else ""
        return f"{pre}{st['smooth']} Hz, p {st['p']}"
    ax.text(0.995, 0.03, f"fit {tag}, score {fit['scores'][0]:.3f}; AsLS experiment {describe(settings)}; "
            f"simulation {describe(settings_m)}", transform=ax.transAxes, fontsize=6.5 * k, color="0.4", ha="right",
            va="bottom")


def legend(ax, site_styles=None):
    import matplotlib.pyplot as plt
    ax.axis("off")
    fontsize = 9 * min(1.0, ax.figure.get_figwidth() / 12) ** 1.5      # narrow figures: smaller legend text
    x = 0.02
    entries = ([(text, color) for text, color in site_styles.values()] if site_styles else
               [(f"{name} (13C-H{'' if nh == 1 else nh}): multiplet centre {center}", color)
                for nh, (name, color, center) in TOPOLOGY.items()])
    for text, color in entries:
        ax.add_patch(plt.Rectangle((x, 0.35), 0.03, 0.3, color=color, alpha=0.45, transform=ax.transAxes))
        ax.text(x + 0.04, 0.5, text, transform=ax.transAxes, va="center", fontsize=fontsize)
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
    width = FIG_WIDTH.get(args.group, 12)
    fig, axs = plt.subplots(len(ids) + 1, 1, figsize=(width, 1.9 * len(ids) + 0.6),
                            gridspec_kw={"height_ratios": [0.3] + [1] * len(ids)}, sharex=False)
    styles = SITE_STYLES.get(args.group)
    legend(axs[0], styles)
    for ax, eid in zip(axs[1:], ids):
        panel(ax, eid, view, styles)
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
