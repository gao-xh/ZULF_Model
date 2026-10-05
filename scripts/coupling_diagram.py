"""Couplings drawn on the molecule: one curve per fitted coupling between the two groups it joins, width growing
with |J|, colour by sign, line style by reliability, labelled J +- sigma.

    python scripts/coupling_diagram.py --fit RUN/fit.json --smiles "CC(N)C" \\
        --atoms '{"C1": 1, "C2": 0, "N1": 2, "HC1": "H@1", "HC2": "H@0", "HC3": "H@3"}' \\
        [--variants RUN2/fit.json RUN3/fit.json ...] [--noise-scale 4] [--title NAME] [--figure OUT.png]

--atoms maps the group labels of the coupling keys (J(a,b)) to atoms of the SMILES (heavy-atom index) or to the
protons on a heavy atom ("H@i"). sigma = sqrt((noise_scale * linearised sigma)^2 + half range over --variants^2):
the linearised error of fit.json scaled for the residual correlation (--noise-scale, e.g. sqrt of the correlation
length in points; the budget of docs/WORKFLOW.md section 9) and the spread over refits with other model settings
(rate families, extra couplings, exchange). Classes: reliable (solid) sigma < 0.05 |J| and < 0.3 Hz; trend
(dashed) sigma < 0.5 |J|; not determined (dotted) otherwise. A refined coupling is a conditional result of the
model; the classes say how well this data set and these model variants fix it.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from figure_tools import GRID, INK, MUTED  # noqa: E402

POSITIVE, NEGATIVE = "#2a78d6", "#eb6834"      # sign of J (validated pair of figure_tools.PALETTE slots)


def coordinates(smiles):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(mol.GetNumAtoms())])
    return mol, xy


def group_atoms(mol, spec):
    if isinstance(spec, str) and spec.startswith("H@"):
        heavy = int(spec[2:])
        return [n.GetIdx() for n in mol.GetAtomWithIdx(heavy).GetNeighbors() if n.GetAtomicNum() == 1]
    return [int(spec)]


def reliability(j, sigma):
    if sigma < 0.05 * abs(j) and sigma < 0.3:
        return "reliable", "-"
    if sigma < 0.5 * abs(j):
        return "trend", (0, (5, 3))
    return "not determined", (0, (1, 2.5))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fit", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--atoms", required=True)
    ap.add_argument("--variants", nargs="*", default=[])
    ap.add_argument("--noise-scale", type=float, default=1.0)
    ap.add_argument("--names", default="{}", help="JSON {key: display name}, e.g. {\"J(C1,HC1)\": \"1J(C1,H1)\"}")
    ap.add_argument("--title", default="")
    ap.add_argument("--figure", default="")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath

    fit = json.load(open(args.fit))
    atoms = json.loads(args.atoms)
    names = json.loads(args.names)
    variants = [json.load(open(p)) for p in args.variants]
    mol, xy = coordinates(args.smiles)
    rows = []
    for key, c in fit["couplings"].items():
        a, b = key[2:-1].split(",")
        if a not in atoms or b not in atoms:
            continue
        j = float(c["J_at_x"][0])
        lin = float(c.get("J_std_at_x", [0.0])[0]) * args.noise_scale
        others = [float(v["couplings"][key]["J_at_x"][0]) for v in variants if key in v["couplings"]] + [j]
        spread = 0.5 * (max(others) - min(others)) if len(others) > 1 else 0.0
        sigma = float(np.hypot(lin, spread))
        rows.append((key, a, b, j, sigma, lin, spread, len(others)))

    # group graph: heavy atoms at their 2D positions (bond length about 1.6), every proton group as one node
    # placed outwards from its carbon (away from the carbon's heavy neighbours)
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    hxy = xy[heavy]
    bl = np.median([np.hypot(*(xy[b.GetBeginAtomIdx()] - xy[b.GetEndAtomIdx()])) for b in mol.GetBonds()
                    if mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomicNum() > 1
                    and mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomicNum() > 1])
    P = {i: (xy[i] - hxy.mean(axis=0)) * 2.0 / bl for i in heavy}
    node, node_label = {}, {}
    for g, spec in atoms.items():
        if isinstance(spec, str) and spec.startswith("H@"):
            c = int(spec[2:])
            nb = [P[n.GetIdx()] for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() > 1]
            # in the middle of the widest angular gap between the carbon's heavy-atom bonds
            ang = np.sort([np.arctan2(*(q - P[c])[::-1]) for q in nb]) if nb else np.array([0.0])
            gaps = np.diff(np.r_[ang, ang[0] + 2 * np.pi])
            k = int(np.argmax(gaps))
            t = ang[k] + 0.5 * gaps[k]
            node[g] = P[c] + 1.3 * np.array([np.cos(t), np.sin(t)])
            nh = len(group_atoms(mol, spec))
            node_label[g] = "H" if nh == 1 else f"H$_{nh}$"
        else:
            node[g] = P[int(spec)]
            node_label[g] = g
    fig, ax = plt.subplots(figsize=(10, 8))
    for b in mol.GetBonds():
        i, k = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in P and k in P:
            ax.plot(*np.array([P[i], P[k]]).T, color="#c9c9c5", lw=2.2, zorder=1, solid_capstyle="round")
    for g, spec in atoms.items():
        if isinstance(spec, str) and spec.startswith("H@"):
            ax.plot(*np.array([P[int(spec[2:])], node[g]]).T, color="#c9c9c5", lw=1.4, zorder=1)
    for i in heavy:                                          # heavy atoms without a group label (e.g. a second CH3)
        if i not in [v for v in atoms.values() if not isinstance(v, str)]:
            sym = mol.GetAtomWithIdx(i).GetSymbol()
            ax.text(*P[i], sym, ha="center", va="center", fontsize=12, color=MUTED, zorder=4,
                    bbox=dict(boxstyle="circle,pad=0.25", fc="white", ec="#c9c9c5"))
    for g in atoms:
        is_h = node_label[g].startswith("H")
        ax.text(*node[g], node_label[g], ha="center", va="center", fontsize=12, color=INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.3" if is_h else "circle,pad=0.3", fc="white",
                          ec=INK if not g.startswith("N") else MUTED, lw=1.0))
    # coupling curves
    jmax = max(abs(r[3]) for r in rows)
    centre = np.mean(list(node.values()), axis=0)
    labels = []
    for n, (key, a, b, j, sigma, lin, spread, nv) in enumerate(sorted(rows, key=lambda r: -abs(r[3]))):
        pa, pb = node[a], node[b]
        mid = 0.5 * (pa + pb)
        d = pb - pa
        length = max(np.hypot(*d), 1e-9)
        normal = np.array([-d[1], d[0]]) / length
        if np.dot(normal, mid - centre) < 0:
            normal = -normal
        direct = any(isinstance(atoms[g], str) and atoms[g] == f"H@{atoms[o]}" for g, o in ((a, b), (b, a))
                     if not isinstance(atoms[o], str))
        bow = 0.0 if direct else 0.22 * length + 0.25
        ctrl = mid + bow * normal
        width = 1.0 + 11.0 * np.sqrt(abs(j) / jmax)
        cls, style = reliability(j, sigma)
        path = MPath([pa, ctrl, pb], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
        ax.add_patch(PathPatch(path, fc="none", ec=POSITIVE if j > 0 else NEGATIVE, lw=width, ls=style,
                               alpha=0.7, capstyle="round", zorder=2 if not direct else 3))
        apex = 0.25 * pa + 0.5 * ctrl + 0.25 * pb
        labels.append([apex + (0.5 if direct else 0.38) * normal, f"{names.get(key, key)}\n{j:+.2f} $\\pm$ {sigma:.2f} Hz",
                       apex.copy()])
    # keep labels apart (simple repulsion), then draw
    for _ in range(200):
        moved = False
        for i in range(len(labels)):
            for k in range(i + 1, len(labels)):
                d = labels[i][0] - labels[k][0]
                if abs(d[0]) < 1.25 and abs(d[1]) < 0.62:
                    step = (d if np.hypot(*d) > 1e-6 else np.array([0.0, 1.0])) / max(np.hypot(*d), 1e-6) * 0.05
                    labels[i][0] = labels[i][0] + step
                    labels[k][0] = labels[k][0] - step
                    moved = True
        if not moved:
            break
    for pos, text, anchor in labels:
        if np.hypot(*(pos - anchor)) > 0.45:            # moved away from its curve: a leader line
            ax.plot(*np.array([anchor, pos]).T, color="#9a9a96", lw=0.7, zorder=6)
        ax.text(*pos, text, ha="center", va="center", fontsize=9, color=INK, zorder=7,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=GRID, lw=0.7, alpha=0.95))
    xy = np.array(list(node.values()) + [p for p, _, _ in labels])
    ax.set_aspect("equal")
    ax.axis("off")
    pad = 0.9
    ax.set_xlim(xy[:, 0].min() - pad, xy[:, 0].max() + pad)
    ax.set_ylim(xy[:, 1].min() - pad, xy[:, 1].max() + pad)
    handles = [Line2D([], [], color=POSITIVE, lw=3, alpha=0.75), Line2D([], [], color=NEGATIVE, lw=3, alpha=0.75),
               Line2D([], [], color=MUTED, lw=2.5, ls="-"), Line2D([], [], color=MUTED, lw=2.5, ls=(0, (5, 3))),
               Line2D([], [], color=MUTED, lw=2.5, ls=(0, (1, 2.5)))]
    for v in (1.0, 10.0, 130.0):
        handles.append(Line2D([], [], color=MUTED, lw=1.0 + 11.0 * np.sqrt(min(v, jmax) / jmax), alpha=0.7))
    ax.legend(handles, ["J > 0", "J < 0", "reliable", "trend", "not determined", "|J| = 1 Hz", "10 Hz", "130 Hz"],
              loc="lower left", bbox_to_anchor=(0.0, -0.02), ncol=4, frameon=False, fontsize=9)
    title = args.title or Path(args.fit).parent.name
    ax.set_title(f"{title}: fitted couplings (width ~ sqrt|J|)", fontsize=12, color=INK, loc="left")
    note = (f"sigma = sqrt(({args.noise_scale:g} x linearised)^2 + (half range over {len(variants)} model variants)^2); "
            "classes: reliable sigma < 5 % and < 0.3 Hz, trend sigma < 50 %, else not determined")
    fig.text(0.02, 0.01, note, fontsize=8, color=MUTED)
    out = Path(args.figure or Path(args.fit).parent / "coupling_diagram.png")
    fig.savefig(out, dpi=170, bbox_inches="tight")
    table = [{"key": r[0], "J": r[3], "sigma": r[4], "linearised_scaled": r[5], "variant_half_range": r[6],
              "fits": r[7], "class": reliability(r[3], r[4])[0]} for r in rows]
    out.with_suffix(".json").write_text(json.dumps(table, indent=1))
    for t in table:
        print(f"{t['key']:12s} {t['J']:+9.3f} +- {t['sigma']:.3f}  ({t['class']}; noise {t['linearised_scaled']:.3f}, "
              f"variants {t['variant_half_range']:.3f} over {t['fits']})")
    print(out)


if __name__ == "__main__":
    main()
