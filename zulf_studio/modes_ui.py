"""Right-hand panels of the Studio task modes that are not fits (D58): Simulate and Process.

Simulate works on the model alone (no data): how the spectrum is drawn (weighted sum or every component), what
the model contains, and export of the simulated spectrum. Process turns scans or an averaged FID into a spectrum
(no model); the full recipe editor with live preview and scan selection is PLAN 8b.
"""
from __future__ import annotations

from PySide6.QtWidgets import (QCheckBox, QGroupBox, QHBoxLayout, QLabel, QPushButton, QRadioButton, QVBoxLayout,
                               QWidget)

from .theme import ISOTOPOLOGUE


class SimulatePanel(QWidget):
    """Display and export of the simulation; the model itself is edited on the left."""

    def __init__(self, session, window, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        lay = QVBoxLayout(self)
        disp = QGroupBox("Display")
        dl = QVBoxLayout(disp)
        self.sum_only = QRadioButton("weighted sum of the components")
        self.each = QRadioButton("every component and the sum")
        self.sum_only.setChecked(True)
        for w in (self.sum_only, self.each):
            dl.addWidget(w)
            w.toggled.connect(lambda _: window.schedule())
        self.show_data = QCheckBox("show the loaded spectrum behind the simulation")
        self.show_data.toggled.connect(lambda _: window.schedule())
        dl.addWidget(self.show_data)
        lay.addWidget(disp)
        comp = QGroupBox("Model")
        cl = QVBoxLayout(comp)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        cl.addWidget(self.summary)
        note = QLabel("Quick look: complex Lorentzian lines with the decay rate of the left panel (FWHM = rate / pi). "
                      "Exact rendering through an acquisition (window, zero fill, Gaussian broadening) is planned "
                      "(PLAN 8e).", objectName="hint")
        note.setWordWrap(True)
        cl.addWidget(note)
        lay.addWidget(comp)
        out = QGroupBox("Export")
        ol = QHBoxLayout(out)
        exp = QPushButton("Export simulation ...", objectName="primary")
        exp.clicked.connect(window.export)
        fig = QPushButton("Plot as image ...")
        fig.clicked.connect(window.export_plot)
        ol.addWidget(exp)
        ol.addWidget(fig)
        ol.addStretch(1)
        lay.addWidget(out)
        go = QPushButton("Fit this model to data  \u2192")
        go.clicked.connect(lambda: self.session.set_mode("fit"))
        lay.addWidget(go)
        lay.addStretch(1)

    def per_component(self) -> bool:
        return self.each.isChecked()

    def refresh(self):
        s = self.session
        comps = s.components()
        lines = s.lines(min_relative=0.01)
        parts = []
        for i, c in enumerate(comps):
            n = sum(1 for r in lines if r["component"] == c["label"])
            ab = c.get("abundance", c.get("weight"))
            parts.append(f"<span style='color:{ISOTOPOLOGUE[i % 6]}'><b>{c['label']}</b></span>"
                         + (f" weight {ab:.3g}" if isinstance(ab, (int, float)) else "") + f", {n} lines")
        self.summary.setText("<br>".join(parts) + f"<br>{len(s.couplings())} couplings")


class ProcessPanel(QWidget):
    """Scans or an averaged FID to a spectrum (no model)."""

    def __init__(self, session, window, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        lay = QVBoxLayout(self)
        box = QGroupBox("Source")
        bl = QVBoxLayout(box)
        for text, fn in (("Import scan folder ...", window.import_scans_dialog),
                         ("Import averaged FID ...", window.import_fid_dialog),
                         ("Open a processed spectrum ...", window.open_dialog)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            bl.addWidget(b)
        hint = QLabel("An import inspects the data first (scans, sampling rate, record, sequence, cached averages, "
                      "preview), then averages and processes with the recipe of the dialog; it runs as a job.",
                      objectName="hint")
        hint.setWordWrap(True)
        bl.addWidget(hint)
        lay.addWidget(box)
        self.info = QLabel(objectName="hint")
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        go = QPushButton("Fit a model to this spectrum  \u2192")
        go.clicked.connect(lambda: self.session.set_mode("fit"))
        blind = QPushButton("Blind analysis of this FID  \u2192")
        blind.clicked.connect(lambda: self.session.set_mode("blind"))
        lay.addWidget(go)
        lay.addWidget(blind)
        lay.addStretch(1)

    def refresh(self):
        d = self.session.data
        self.info.setText("no spectrum yet" if d is None else
                          f"current spectrum: {d['label']}, {len(d['freq'])} points, {len(d['ranges'])} fit ranges"
                          + (f"<br>from {d['source_fid']}" if d.get("source_fid") else ""))
