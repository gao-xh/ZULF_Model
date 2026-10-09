"""Right-hand panel of the Studio Simulate mode (D58).

Simulate works on the model alone (no data): how the spectrum is drawn (weighted sum or every component), what
the model contains, and export of the simulated spectrum. Process is in process_ui.py.
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
        self.show_data = QCheckBox("show the loaded data")
        self.show_data.setToolTip("draw the loaded experimental spectrum behind the simulation")
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
        ol = QVBoxLayout(out)
        exp = QPushButton("Export simulation ...", objectName="primary")
        exp.clicked.connect(window.export)
        fig = QPushButton("Plot as image ...")
        fig.clicked.connect(window.export_plot)
        ol.addWidget(exp)
        ol.addWidget(fig)
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
