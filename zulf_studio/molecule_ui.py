"""The model drawn: a molecule (structure models; spin systems that carry a molecule, D61) or a spin network
(typed-in spin systems). Drawn with matplotlib in the plot style (RDKit only computes the 2D coordinates):

- molecule: atoms as text (CH3, CH2, OH; protons spelled out because they carry the couplings), bonds as lines
  (double and triple as parallel lines), the site label (C1, C2, ...) small and grey outside the molecule, and the
  labelled site of every isotopologue in a pill of its plot colour (the colours of the spectrum and Lines table);
- spin network: one node per group of equivalent spins, one edge per nonzero coupling with its value, width
  growing with |J|, dashed when negative.

MoleculeView is the thumbnail in the Model card (compact: no caption, a click opens the large view) and the large
view itself.
"""
from __future__ import annotations

import re

import matplotlib
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QComboBox, QDialog, QHBoxLayout, QInputDialog, QLabel, QMessageBox,
                               QVBoxLayout, QWidget)

from .theme import ISOTOPOLOGUE, matplotlib_style

# spin-network node colours by nucleus (isotopes of one element in shades of one hue); others by element, then grey
NUCLEUS_COLOURS = {"1H": "#8a96a8", "2H": "#4f5d73", "3H": "#2f3a4d", "13C": "#2f8f5b", "15N": "#3b6fb6",
                   "14N": "#7fa3d6", "17O": "#d1495b", "19F": "#9ac23c", "31P": "#e08a1e", "29Si": "#b8936a",
                   "77Se": "#b07cc6", "119Sn": "#7c8c99", "6Li": "#c26ab0", "7Li": "#a04f92", "11B": "#e0a3a3",
                   "10B": "#c27a7a", "23Na": "#8d6cc4", "27Al": "#9aa3ad", "33S": "#c8a400", "35Cl": "#3fae8f"}
ELEMENT_COLOURS = {"O": "#d1495b", "N": "#3b6fb6", "S": "#c8a400", "P": "#d1495b", "F": "#2f8f5b", "Cl": "#2f8f5b"}


def label_sites(label: str):
    """Sites of an isotopologue label: '13C@C1 (x2)' -> ['C1'], '13C@C1,C3+15N@N1' -> ['C1', 'C3', 'N1']."""
    out = []
    for part in label.split(" ")[0].split("+"):
        if "@" in part:
            out += [s for s in part.split("@", 1)[1].split(",") if s]
    return out


def _mix(colour, background, share):
    c = np.array(matplotlib.colors.to_rgb(colour))
    b = np.array(matplotlib.colors.to_rgb(background))
    return tuple(share * c + (1 - share) * b)


def draw_molecule(ax, spec, components, t, scale=1.0, explicit_h=True) -> dict:
    """Draw the molecule of a structure specification on ax; components: isotopologue labels (their sites are
    marked in the isotopologue colours, the labelled atom written 13C). explicit_h: every proton its own atom
    (exchangeable ones on O and S faint: they are not in the spin model); else grouped (CH3, CH2).
    Returns {"atoms": heavy atoms, "hydrogens": drawn H atoms, "marked": {site: colour}}."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from zulf_hypothesis.molecule import molecule_from_structure
    mol, sites = molecule_from_structure(spec)
    heavy = mol.GetNumAtoms()
    if explicit_h:
        mol = Chem.AddHs(mol)                          # H atoms appended after the heavy atoms
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    xy = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in range(n)])
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), str(b.GetBondType())) for b in mol.GetBonds()]
    heavy_bonds = [(i, j) for i, j, _ in bonds if i < heavy and j < heavy]
    ref = heavy_bonds or [(i, j) for i, j, _ in bonds]
    length = np.median([np.hypot(*(xy[i] - xy[j])) for i, j in ref]) if ref else 1.0
    xy = xy / length                                   # bond length 1
    label_of = {i: s for s, i in sites.items()}
    marked = {}
    for k, comp in enumerate(components):
        for site in label_sites(comp):
            if site in sites:
                marked[site] = ISOTOPOLOGUE[k % len(ISOTOPOLOGUE)]
    # font from the pixels per bond length, so labels never cover the bonds whatever the view size
    pad = 0.5
    span_x = np.ptp(xy[:, 0]) + 2 * pad
    span_y = np.ptp(xy[:, 1]) + 2 * pad
    w_in, h_in = ax.figure.get_size_inches()
    bbox = ax.get_position()
    px_per_unit = min(w_in * bbox.width / span_x, h_in * bbox.height / span_y) * 72.0   # in points
    fs = float(np.clip(0.34 * px_per_unit, 5.5, 15.0 * scale))

    def faint(i):                                      # exchangeable proton: not in the spin model
        a = mol.GetAtomWithIdx(i)
        return a.GetAtomicNum() == 1 and any(n.GetSymbol() in ("O", "S") for n in a.GetNeighbors())
    for i, j, kind in bonds:
        p, q = xy[i], xy[j]
        u = (q - p) / max(np.hypot(*(q - p)), 1e-9)
        nrm = np.array([-u[1], u[0]])
        offsets = {"DOUBLE": (-0.07, 0.07), "TRIPLE": (-0.11, 0.0, 0.11), "AROMATIC": (-0.07, 0.07)}.get(kind, (0,))
        for m, off in enumerate(offsets):
            a, b = p + off * nrm, q + off * nrm
            h_bond = i >= heavy or j >= heavy
            ax.plot([a[0], b[0]], [a[1], b[1]], color=t["muted"] if h_bond else t["ink"],
                    lw=(1.0 if h_bond else 1.5) * scale, solid_capstyle="round", zorder=1,
                    alpha=0.35 if (faint(i) or faint(j)) else 1.0,
                    ls=(0, (3, 2)) if kind == "AROMATIC" and m == 1 else "-")
    for i, atom in enumerate(mol.GetAtoms()):
        el, n_h = atom.GetSymbol(), (0 if explicit_h else atom.GetTotalNumHs())
        site = label_of.get(i)
        colour = marked.get(site)
        if el == "H":
            ax.text(*xy[i], "H", ha="center", va="center", fontsize=fs * 0.8, zorder=3, color=t["ink"],
                    alpha=0.35 if faint(i) else 0.85, bbox=dict(boxstyle="round,pad=0.12", fc=t["panel"], ec="none"))
            continue
        iso = "".join(c for c in comp_isotope(components, site) if c.isdigit()) if colour else ""
        text = (f"$^{{{iso}}}$" if iso else "") + el + ("H" if n_h else "") + (f"$_{{{n_h}}}$" if n_h > 1 else "")
        box = (dict(boxstyle="round,pad=0.22,rounding_size=0.5", fc=_mix(colour, t["panel"], 0.22), ec=colour,
                    lw=1.4) if colour else
               dict(boxstyle="round,pad=0.2", fc=t["panel"], ec="none"))
        ax.text(*xy[i], text, ha="center", va="center", fontsize=fs, zorder=3,
                color=ELEMENT_COLOURS.get(el, t["ink"]), fontweight="bold" if colour else "normal", bbox=box)
        if site:                                       # the site label in the freest direction around the atom
            nb = [(xy[a.GetIdx()] - xy[i]) / max(np.hypot(*(xy[a.GetIdx()] - xy[i])), 1e-9)
                  for a in atom.GetNeighbors()]
            cand = [np.array([np.cos(a), np.sin(a)]) for a in np.linspace(0, 2 * np.pi, 16, endpoint=False)]
            d = max(cand, key=lambda c: min([np.hypot(*(c - v)) for v in nb] or [2.0]))
            ax.text(*(xy[i] + 0.5 * d), site, ha="center", va="center", fontsize=fs * 0.58, color=t["muted"],
                    zorder=4)
    ax.plot(                                           # invisible corners: equal aspect without fixed limits
        [xy[:, 0].min() - pad, xy[:, 0].max() + pad], [xy[:, 1].min() - pad, xy[:, 1].max() + pad], alpha=0)
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale(tight=True)
    ax.set_axis_off()
    return {"atoms": heavy, "hydrogens": n - heavy, "marked": marked}


def comp_isotope(components, site) -> str:
    """The isotope labelling a site in the isotopologue labels ('13C@C1' -> '13C')."""
    for comp in components:
        for part in comp.split(" ")[0].split("+"):
            if "@" in part and site in part.split("@", 1)[1].split(","):
                return part.split("@", 1)[0]
    return ""


def nucleus_colour(nucleus: str) -> str:
    """Colour of a nucleus ('13C', '2H', ...) in the spin network."""
    if nucleus in NUCLEUS_COLOURS:
        return NUCLEUS_COLOURS[nucleus]
    el = re.sub(r"^\d+", "", nucleus)
    return ELEMENT_COLOURS.get(el, "#8a8a8f")


def draw_network(ax, spec, k, t, scale=1.0) -> dict:
    """Draw component k of a spin-system specification as a spin network on ax."""
    from zulf_hypothesis.spin_system import groups, numeric_matrix, tokens
    comp = spec["spin_system"]["components"][k]
    iso = list(comp["isotopes"])
    grp = groups(iso, tokens(comp))
    m = numeric_matrix(comp, spec["spin_system"].get("variables", {}))
    n = len(grp)
    ang = np.linspace(0.5 * np.pi, 2.5 * np.pi, n, endpoint=False)
    xy = np.c_[np.cos(ang), np.sin(ang)] if n > 1 else np.zeros((1, 2))
    js = [abs(m[ga[0], gb[0]]) for a, ga in enumerate(grp) for b, gb in enumerate(grp) if b > a]
    top = max(js + [1e-9])
    for a in range(n):
        for b in range(a + 1, n):
            j = m[grp[a][0], grp[b][0]]
            if j == 0:
                continue
            (x0, y0), (x1, y1) = xy[a], xy[b]
            ax.plot([x0, x1], [y0, y1], color=t["accent"] if j > 0 else t["muted"],
                    lw=(0.6 + 3.4 * np.sqrt(abs(j) / top)) * scale, alpha=0.85, zorder=1,
                    ls="-" if j > 0 else (0, (4, 2)))
            if n <= 8 or abs(j) >= 1.0:
                ax.text(0.6 * x0 + 0.4 * x1, 0.6 * y0 + 0.4 * y1, f"{j:.4g}", fontsize=7.5 * scale, ha="center",
                        va="center", color=t["ink"], zorder=3,
                        bbox=dict(boxstyle="round,pad=0.15", fc=t["panel"], ec="none", alpha=0.9))
    colours = {}
    for a, g in enumerate(grp):
        col = colours.setdefault(iso[g[0]], nucleus_colour(iso[g[0]]))
        text = iso[g[0]] + (f" x{len(g)}" if len(g) > 1 else "")
        ax.text(*xy[a], text, ha="center", va="center", fontsize=9 * scale, color=t["ink"], zorder=4,
                bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.8", fc=_mix(col, t["panel"], 0.3), ec=col,
                          lw=1.4))
    for m, (nuc, col) in enumerate(colours.items()):     # legend: the nuclei present
        ax.text(-1.45 + 0.62 * m, -1.3, f"\u25cf {nuc}", color=col, fontsize=7.5 * scale, ha="left", va="center")
    ax.plot([-1.5, 1.5], [-1.35, 1.35], alpha=0)
    ax.set_aspect("equal", adjustable="datalim")
    ax.autoscale(tight=True)
    ax.set_axis_off()
    return {"spins": len(iso), "groups": n, "name": comp.get("name", ""), "colours": colours}


class MoleculeView(QWidget):
    """Thumbnail (compact, in the Model card) or large view of the model drawing."""

    def __init__(self, session, window, parent=None, large=False):
        super().__init__(parent)
        self.session, self.window, self.large = session, window, large
        self.big = None
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.fig = Figure(figsize=(3.6, 2.0))
        self.fig.subplots_adjust(0.01, 0.01, 0.99, 0.99)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumHeight(420 if large else 120)
        if not large:
            self.canvas.setMaximumHeight(150)
            self.canvas.setCursor(Qt.PointingHandCursor)
            self.canvas.mpl_connect("button_press_event", lambda _e: self.open_large())
        self.canvas.mpl_connect("resize_event", lambda _e: self.refresh(force=True))   # fonts follow the size
        lay.addWidget(self.canvas, 1)
        self.caption = QLabel(objectName="hint", wordWrap=True)
        self.caption.setVisible(large)
        lay.addWidget(self.caption)
        row = self.row = QHBoxLayout()
        self.view = QComboBox()
        self.view.addItems(["molecule", "spin network"])
        self.view.setToolTip("a spin system with a molecule: draw the molecule or the spin network")
        self.view.currentIndexChanged.connect(lambda _i: self.refresh(force=True))
        self.component = QComboBox()
        self.component.setToolTip("component of the spin system drawn")
        self.component.currentIndexChanged.connect(lambda _i: self.refresh(force=True))
        self.hydrogens = QComboBox()
        self.hydrogens.addItems(["H atoms", "CH3 groups"])
        self.hydrogens.setToolTip("protons drawn as atoms, or grouped on their atom (CH3, CH2)")
        self.hydrogens.currentIndexChanged.connect(lambda _i: self.refresh(force=True))
        row.addWidget(self.view)
        row.addWidget(self.component)
        row.addWidget(self.hydrogens)
        row.addStretch(1)
        lay.addLayout(row)
        self._last = None

    def refresh(self, force=False):
        spec = self.session.spec
        t = self.window.t
        key = (repr(spec), self.component.currentIndex(), self.view.currentIndex(), self.hydrogens.currentIndex(),
               t["panel"])
        if key == self._last and not force:
            return
        self._last = key
        spin = "spin_system" in spec
        has_mol = bool(spec.get("molecule")) or not spin
        self.view.setVisible(spin and has_mol)
        network = spin and (not has_mol or self.view.currentIndex() == 1)
        self.hydrogens.setVisible(not network)
        scale = 1.5 if self.large else 0.85
        with matplotlib.rc_context(matplotlib_style(t)):
            self.fig.clear()
            self.fig.set_facecolor(t["panel"])
            ax = self.fig.add_subplot(111)
            try:
                if network:
                    comps = spec["spin_system"]["components"]
                    names = [c.get("name", f"component {i + 1}") for i, c in enumerate(comps)]
                    if [self.component.itemText(i) for i in range(self.component.count())] != names:
                        self.component.blockSignals(True)
                        self.component.clear()
                        self.component.addItems(names)
                        self.component.blockSignals(False)
                    self.component.setVisible(len(names) > 1)
                    r = draw_network(ax, spec, max(self.component.currentIndex(), 0), t, scale)
                    text = (f"{r['name']}: {r['spins']} spins in {r['groups']} groups of equivalent spins; edge "
                            "width ~ |J|, dashed: negative J (Hz)"
                            + ("" if has_mol else "; Edit model > Attach molecule draws it as a molecule"))
                else:
                    self.component.setVisible(False)
                    comps = ([c.get("name", "") for c in spec["spin_system"]["components"]] if spin
                             else [c["label"] for c in self.session.components()])
                    draw_molecule(ax, spec, comps, t, scale, explicit_h=self.hydrogens.currentIndex() == 0)
                    smiles = spec.get("molecule", {}).get("smiles")
                    text = ((f"SMILES {smiles}" if smiles else "skeleton (bond orders are not part of the "
                             "structure)") + "; ringed: the labelled site of each isotopologue, in its plot colour")
            except Exception as exc:
                ax.set_axis_off()
                ax.text(0.5, 0.5, "no drawing", ha="center", va="center", color=t["muted"], transform=ax.transAxes)
                text = f"no drawing: {type(exc).__name__}: {exc}"
        self.caption.setText(text)
        self.canvas.setToolTip(text + ("" if self.large else "\nclick: large view"))
        self.canvas.draw_idle()

    def ask_molecule(self, parent=None):
        parent = parent or self
        text, ok = QInputDialog.getMultiLineText(parent, "Model from a molecule",
                                                 "SMILES (e.g. CCO for ethanol) or the text of a mol file:")
        if not ok or not text.strip():
            return
        name, ok = QInputDialog.getText(parent, "Model from a molecule", "compound name (optional):")
        try:
            r = self.session.structure_from_molecule(text, compound=name.strip() or None)
        except Exception as exc:
            QMessageBox.warning(parent, "Molecule", f"{type(exc).__name__}: {exc}")
            return
        if r["notes"]:
            QMessageBox.information(parent, "Molecule", "\n".join(r["notes"]))

    def ask_attach(self, parent=None):
        parent = parent or self
        text, ok = QInputDialog.getMultiLineText(parent, "Draw as a molecule",
                                                 "SMILES (e.g. CCO for ethanol: C1, C2, O1) or the text of a mol "
                                                 "file; drawing only, the model is unchanged:")
        if not ok or not text.strip():
            return
        try:
            self.session.attach_molecule(text)
        except Exception as exc:
            QMessageBox.warning(parent, "Molecule", f"{type(exc).__name__}: {exc}")
            return
        self.view.setCurrentIndex(0)

    def open_large(self):
        if self.big is None:
            self.big = QDialog(self.window)
            self.big.setWindowTitle("ZULF Studio - model")
            lay = QVBoxLayout(self.big)
            self.big.view = MoleculeView(self.session, self.window, self.big, large=True)
            lay.addWidget(self.big.view)
            screen = self.window.screen() or QApplication.primaryScreen()
            if screen is not None:
                g = screen.availableGeometry()
                self.big.resize(int(g.width() * 0.5), int(g.height() * 0.6))
        self.big.view.component.setCurrentIndex(self.component.currentIndex())
        self.big.view.view.setCurrentIndex(self.view.currentIndex())
        self.big.view.hydrogens.setCurrentIndex(self.hydrogens.currentIndex())
        self.big.view.refresh(force=True)
        self.big.show()
        self.big.raise_()
