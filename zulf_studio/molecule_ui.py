"""Molecule card: the model drawn as a molecule (structure models, RDKit) or as a spin network (typed-in spin
systems). Structure: the heavy-atom skeleton with its protons (or the recorded SMILES of a molecule built with
"From molecule"), site labels, and the 13C (15N) position of every isotopologue in its plot colour. Spin system:
one node per group of equivalent spins, one edge per nonzero coupling with its value, width growing with |J|.
"""
from __future__ import annotations

import re

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QByteArray, Qt
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (QApplication, QComboBox, QDialog, QHBoxLayout, QInputDialog, QLabel, QMessageBox,
                               QPushButton, QStackedWidget, QVBoxLayout, QWidget)

from .theme import ISOTOPOLOGUE, matplotlib_style


def _rgb(hex_colour, alpha=None):
    h = hex_colour.lstrip("#")
    rgb = tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return rgb + ((alpha,) if alpha is not None else ())


def label_sites(label: str):
    """Sites of an isotopologue label: '13C@C1 (x2)' -> ['C1'], '13C@C1,C3+15N@N1' -> ['C1', 'C3', 'N1']."""
    out = []
    for part in label.split(" ")[0].split("+"):
        if "@" in part:
            out += [s for s in part.split("@", 1)[1].split(",") if s]
    return out


def structure_svg(spec, components, theme, size=(360, 240)) -> str:
    """SVG drawing of a structure specification (zulf_hypothesis.molecule.molecule_from_structure)."""
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdDepictor
    from zulf_hypothesis.molecule import molecule_from_structure
    mol, sites = molecule_from_structure(spec)
    rdDepictor.Compute2DCoords(mol)
    for label, i in sites.items():
        atom = mol.GetAtomWithIdx(i)
        atom.SetProp("atomNote", label)
        n_h = atom.GetTotalNumHs()
        el = atom.GetSymbol()
        if el == "C" or n_h:                         # carbons with their protons spelled out (CH3, CH2, C)
            atom.SetProp("_displayLabel", el + ("H" if n_h else "") + (f"<sub>{n_h}</sub>" if n_h > 1 else ""))
    colours, atoms = {}, []
    for k, comp in enumerate(components):
        for site in label_sites(comp):
            if site in sites:
                atoms.append(sites[site])
                # opaque: Qt's SVG renderer draws RDKit's #RRGGBBAA colours black, so blend with the background
                c, bg = np.array(_rgb(ISOTOPOLOGUE[k % len(ISOTOPOLOGUE)])), np.array(_rgb(theme["panel"]))
                colours[sites[site]] = tuple(float(v) for v in 0.45 * c + 0.55 * bg)
    d = rdMolDraw2D.MolDraw2DSVG(*size)
    o = d.drawOptions()
    o.setBackgroundColour(_rgb(theme["panel"], 1.0))
    o.updateAtomPalette({-1: _rgb(theme["ink"]), 6: _rgb(theme["ink"]), 1: _rgb(theme["ink"])})
    o.setAnnotationColour(_rgb(theme["muted"]))
    o.annotationFontScale = 0.6
    o.highlightRadius = 0.45
    o.fixedBondLength = 38
    o.clearBackground = True
    d.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors=colours, highlightBonds=[])
    d.FinishDrawing()
    return d.GetDrawingText()


class MoleculeView(QWidget):
    """The Molecule card (and its large window)."""

    def __init__(self, session, window, parent=None, large=False):
        super().__init__(parent)
        self.session, self.window, self.large = session, window, large
        self.big = None
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.stack = QStackedWidget()
        self.svg = QSvgWidget()
        self.svg.renderer().setAspectRatioMode(Qt.KeepAspectRatio)
        self.fig = Figure(figsize=(3.6, 2.4), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        for w in (self.svg, self.canvas):
            w.setMinimumHeight(420 if large else 170)
        self.stack.addWidget(self.svg)
        self.stack.addWidget(self.canvas)
        lay.addWidget(self.stack, 1)
        self.caption = QLabel(objectName="hint", wordWrap=True)
        lay.addWidget(self.caption)
        row = QHBoxLayout()
        self.view = QComboBox()
        self.view.addItems(["molecule", "spin network"])
        self.view.setToolTip("a spin system with a molecule: draw the molecule or the spin network")
        self.view.currentIndexChanged.connect(lambda _i: self.refresh())
        row.addWidget(self.view)
        self.component = QComboBox()
        self.component.setToolTip("component of the spin system drawn")
        self.component.currentIndexChanged.connect(lambda _i: self.refresh())
        row.addWidget(self.component)
        row.addStretch(1)
        if not large:
            self.b_smiles = QPushButton("From molecule ...", objectName="small")
            self.b_smiles.setToolTip("build the model from a SMILES string or a mol file (every heavy atom a site, "
                                     "symmetry found automatically, 1J guesses from the hybridisation)")
            self.b_smiles.clicked.connect(self.ask_molecule)
            self.b_attach = QPushButton("Attach molecule ...", objectName="small")
            self.b_attach.setToolTip("draw this spin system as a molecule: a SMILES string or a mol file whose "
                                     "atoms are labelled C1, C2, O1, ... in order (so 13C@C1 marks its atom); "
                                     "the model and its couplings are unchanged")
            self.b_attach.clicked.connect(self.ask_attach)
            row.addWidget(self.b_attach)
            self.b_large = QPushButton("Large", objectName="small")
            self.b_large.clicked.connect(self.open_large)
            row.addWidget(self.b_smiles)
            row.addWidget(self.b_large)
        lay.addLayout(row)
        self._last = None

    def refresh(self):
        spec = self.session.spec
        t = self.window.t
        key = (repr(spec), self.component.currentIndex(), self.view.currentIndex(), t["panel"])
        if key == self._last:
            return
        self._last = key
        spin = "spin_system" in spec
        has_mol = bool(spec.get("molecule")) or not spin
        self.view.setVisible(spin and has_mol)
        if not self.large:
            self.b_smiles.setVisible(not spin)
            self.b_attach.setVisible(spin)
        if spin and (not has_mol or self.view.currentIndex() == 1):
            comps = spec["spin_system"]["components"]
            names = [c.get("name", f"component {i + 1}") for i, c in enumerate(comps)]
            if [self.component.itemText(i) for i in range(self.component.count())] != names:
                self.component.blockSignals(True)
                self.component.clear()
                self.component.addItems(names)
                self.component.blockSignals(False)
            self.component.setVisible(len(names) > 1)
            self.stack.setCurrentWidget(self.canvas)
            self._draw_network(spec, max(self.component.currentIndex(), 0))
            if not has_mol:
                self.caption.setText(self.caption.text() + "; Attach molecule ... draws it as a molecule")
            return
        self.component.setVisible(False)
        self.stack.setCurrentWidget(self.svg)
        comps = ([c.get("name", "") for c in spec["spin_system"]["components"]] if spin
                 else [c["label"] for c in self.session.components()])
        try:
            self.svg.load(QByteArray(structure_svg(spec, comps, t, (720, 480) if self.large else (360, 240))
                                     .encode()))
            mol = spec.get("molecule", {}).get("smiles")
            self.caption.setText((f"SMILES {mol}" if mol else "skeleton: bond orders are not part of the "
                                  "structure") + "; coloured: the labelled site of each isotopologue")
        except Exception as exc:
            self.svg.load(QByteArray(b"<svg xmlns='http://www.w3.org/2000/svg'/>"))
            self.caption.setText(f"no drawing: {type(exc).__name__}: {exc}")

    def _draw_network(self, spec, k):
        from zulf_hypothesis.spin_system import groups, numeric_matrix, tokens
        import matplotlib
        t = self.window.t
        comp = spec["spin_system"]["components"][k]
        iso = list(comp["isotopes"])
        tok = tokens(comp)
        grp = groups(iso, tok)
        m = numeric_matrix(comp, spec["spin_system"].get("variables", {}))
        n = len(grp)
        with matplotlib.rc_context(matplotlib_style(t)):
            self.fig.clear()
            self.fig.set_facecolor(t["panel"])
            ax = self.fig.add_subplot(111)
            ax.set_axis_off()
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
                            lw=0.6 + 3.4 * np.sqrt(abs(j) / top), alpha=0.85, zorder=1,
                            ls="-" if j > 0 else (0, (4, 2)))
                    if n <= 8 or abs(j) >= 1.0:
                        ax.text(0.5 * (x0 + x1), 0.5 * (y0 + y1), f"{j:.4g}", fontsize=7.5, ha="center",
                                va="center", color=t["ink"], zorder=3,
                                bbox=dict(boxstyle="round,pad=0.15", fc=t["panel"], ec="none", alpha=0.85))
            for a, g in enumerate(grp):
                el = re.sub(r"^\d+", "", iso[g[0]])
                col = ISOTOPOLOGUE[0] if el == "C" else ISOTOPOLOGUE[2] if el == "N" else t["panel2"]
                ax.scatter(*xy[a], s=900 if self.large else 520, color=col, edgecolor=t["line2"], zorder=2)
                text = iso[g[0]] + (f" x{len(g)}" if len(g) > 1 else "")
                ax.text(*xy[a], text, ha="center", va="center", fontsize=8.5, color=t["ink"], zorder=4)
            ax.set_xlim(-1.45, 1.45)
            ax.set_ylim(-1.3, 1.3)
            ax.set_aspect("equal")
        self.canvas.draw_idle()
        self.caption.setText(f"{comp.get('name', '')}: {len(iso)} spins in {n} groups of equivalent spins; "
                             "edge width ~ |J|, dashed: negative J (Hz)")

    def ask_molecule(self):
        text, ok = QInputDialog.getMultiLineText(self, "Model from a molecule",
                                                 "SMILES (e.g. CCO for ethanol) or the text of a mol file:")
        if not ok or not text.strip():
            return
        name, ok = QInputDialog.getText(self, "Model from a molecule", "compound name (optional):")
        try:
            r = self.session.structure_from_molecule(text, compound=name.strip() or None)
        except Exception as exc:
            QMessageBox.warning(self, "Molecule", f"{type(exc).__name__}: {exc}")
            return
        if r["notes"]:
            QMessageBox.information(self, "Molecule", "\n".join(r["notes"]))

    def ask_attach(self):
        text, ok = QInputDialog.getMultiLineText(self, "Draw as a molecule",
                                                 "SMILES (e.g. CCO for ethanol: C1, C2, O1) or the text of a mol "
                                                 "file; drawing only, the model is unchanged:")
        if not ok or not text.strip():
            return
        try:
            self.session.attach_molecule(text)
        except Exception as exc:
            QMessageBox.warning(self, "Molecule", f"{type(exc).__name__}: {exc}")
            return
        self.view.setCurrentIndex(0)

    def open_large(self):
        if self.big is None:
            self.big = QDialog(self.window)
            self.big.setWindowTitle("ZULF Studio - molecule")
            lay = QVBoxLayout(self.big)
            self.big.view = MoleculeView(self.session, self.window, self.big, large=True)
            lay.addWidget(self.big.view)
            screen = self.window.screen() or QApplication.primaryScreen()
            if screen is not None:
                g = screen.availableGeometry()
                self.big.resize(int(g.width() * 0.5), int(g.height() * 0.6))
        self.big.view.component.setCurrentIndex(self.component.currentIndex())
        self.big.view.view.setCurrentIndex(self.view.currentIndex())
        self.big.view._last = None
        self.big.view.refresh()
        self.big.show()
        self.big.raise_()
