"""ZULF Studio window (PySide6): live sliders, plot, line table, fitting with a progress slider, log, terminal,
Python console and the AI API, all on one StudioSession.

    python scripts/zulf_studio.py [--structure JSON] [--series SERIES.json] [--fit RUN_DIR] [--api-port 8766]
"""
from __future__ import annotations

import code
import contextlib
import io
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, Qt, QTimer, Signal
from PySide6.QtCore import QUrl
from PySide6.QtGui import QAction, QDesktopServices, QFont, QKeySequence, QPixmap
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPlainTextEdit, QPushButton, QScrollArea, QSlider, QSpinBox, QSplitter,
                               QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from .api import TOOLS, StudioAPI, serve  # noqa: E402
from .session import ROOT, StudioSession  # noqa: E402

COMPONENT_COLORS = ["#d1495b", "#3b6fb6", "#2a9d5c", "#e08a2c", "#8e5ac8", "#7a7a7a"]
MONO = QFont("Menlo", 11)


class Bridge(QObject):
    """Session events (any thread: API server, fit reader) -> Qt signals handled in the GUI thread."""
    changed = Signal(str)
    logged = Signal(dict)


class ValueSlider(QWidget):
    """A labelled float slider (linear or log) with a spin box; emits value_changed(float) while dragging and
    committed(float) on release."""
    value_changed = Signal(float)

    def __init__(self, label, lo, hi, value, decimals=3, unit="", log=False, fine=0.0, parent=None):
        super().__init__(parent)
        self.lo, self.hi, self.log, self.fine_span, self._anchor = lo, hi, log, fine, value
        grid = QGridLayout(self)
        grid.setContentsMargins(0, 2, 0, 2)
        grid.setHorizontalSpacing(6)
        grid.setVerticalSpacing(0)
        self.name = QLabel(label)
        self.name.setMinimumWidth(110)
        self.spin = QDoubleSpinBox()
        self.spin.setDecimals(decimals)
        self.spin.setRange(min(lo, value), max(hi, value))
        self.spin.setSingleStep(10 ** -min(decimals, 2))
        self.spin.setKeyboardTracking(False)
        self.spin.setSuffix(f" {unit}" if unit else "")
        self.spin.setMinimumWidth(110)
        self.coarse = QSlider(Qt.Horizontal)
        self.coarse.setRange(0, 1000)
        grid.addWidget(self.name, 0, 0)
        grid.addWidget(self.coarse, 0, 1)
        grid.addWidget(self.spin, 0, 2)
        self.fine = None
        if fine:
            self.fine = QSlider(Qt.Horizontal)
            self.fine.setRange(-500, 500)
            fl = QLabel(f"fine +-{fine:g}")
            fl.setStyleSheet("color: gray; font-size: 10px")
            grid.addWidget(fl, 1, 0)
            grid.addWidget(self.fine, 1, 1)
            self.fine.sliderPressed.connect(lambda: setattr(self, "_anchor", self.spin.value()))
            self.fine.valueChanged.connect(self._fine_moved)
            self.fine.sliderReleased.connect(lambda: self._reset_fine())
        grid.setColumnStretch(1, 1)
        self.coarse.valueChanged.connect(self._coarse_moved)
        self.spin.valueChanged.connect(self._spin_changed)
        self.set_value(value)

    def _to_pos(self, v):
        if self.log:
            a, b = math.log10(self.lo), math.log10(self.hi)
            return int(round(1000 * (math.log10(max(v, self.lo)) - a) / (b - a)))
        return int(round(1000 * (v - self.lo) / (self.hi - self.lo)))

    def _from_pos(self, p):
        if self.log:
            a, b = math.log10(self.lo), math.log10(self.hi)
            return 10 ** (a + p / 1000 * (b - a))
        return self.lo + p / 1000 * (self.hi - self.lo)

    def _coarse_moved(self, p):
        v = self._from_pos(p)
        self._anchor = v
        self._emit(v)

    def _fine_moved(self, p):
        if self.fine.isSliderDown():
            self._emit(self._anchor + p / 500 * self.fine_span)

    def _reset_fine(self):
        self.fine.blockSignals(True)
        self.fine.setValue(0)
        self.fine.blockSignals(False)
        self._anchor = self.spin.value()

    def _spin_changed(self, v):
        self._anchor = v
        self.coarse.blockSignals(True)
        self.coarse.setValue(self._to_pos(v))
        self.coarse.blockSignals(False)
        self.value_changed.emit(float(v))

    def _emit(self, v):
        self.spin.blockSignals(True)
        if v < self.spin.minimum() or v > self.spin.maximum():
            self.spin.setRange(min(v, self.spin.minimum()), max(v, self.spin.maximum()))
        self.spin.setValue(v)
        self.spin.blockSignals(False)
        self.value_changed.emit(float(v))

    def set_value(self, v):
        for w in (self.spin, self.coarse):
            w.blockSignals(True)
        if v < self.spin.minimum() or v > self.spin.maximum():
            self.spin.setRange(min(v, self.spin.minimum()), max(v, self.spin.maximum()))
        self.spin.setValue(v)
        if not self.coarse.isSliderDown():
            self.coarse.setValue(self._to_pos(v))
        for w in (self.spin, self.coarse):
            w.blockSignals(False)
        if self.fine is None or not self.fine.isSliderDown():
            self._anchor = v

    def value(self):
        return self.spin.value()


class Console(QWidget):
    """Python console with `session`, `api` and `np` in its namespace."""

    def __init__(self, namespace, parent=None):
        super().__init__(parent)
        self.interp = code.InteractiveConsole(namespace)
        lay = QVBoxLayout(self)
        self.out = QPlainTextEdit(readOnly=True)
        self.out.setFont(MONO)
        self.inp = QLineEdit()
        self.inp.setFont(MONO)
        self.inp.setPlaceholderText(">>> python (session, api, np); e.g. session.set_field(40, 45)")
        lay.addWidget(self.out)
        lay.addWidget(self.inp)
        self.history, self.hpos = [], 0
        self.inp.returnPressed.connect(self.run)
        self.out.appendPlainText("Python console: `session` (StudioSession), `api` (StudioAPI), `np`.")

    def run(self):
        line = self.inp.text()
        self.inp.clear()
        self.history.append(line)
        self.hpos = len(self.history)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            more = self.interp.push(line)
        self.out.appendPlainText(("... " if more else ">>> ") + line)
        if buf.getvalue():
            self.out.appendPlainText(buf.getvalue().rstrip())

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key_Up, Qt.Key_Down) and self.history:
            self.hpos = max(0, min(len(self.history), self.hpos + (-1 if e.key() == Qt.Key_Up else 1)))
            self.inp.setText(self.history[self.hpos] if self.hpos < len(self.history) else "")
        else:
            super().keyPressEvent(e)


class Terminal(QWidget):
    """Shell commands in the repository directory, with the studio's Python first on PATH."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        lay = QVBoxLayout(self)
        self.out = QPlainTextEdit(readOnly=True)
        self.out.setFont(MONO)
        row = QHBoxLayout()
        self.inp = QLineEdit()
        self.inp.setFont(MONO)
        self.inp.setPlaceholderText(f"$ shell command in {ROOT} (zsh), e.g. git status, python scripts/fit_monitor.py runs/studio/fits --text")
        self.stop_btn = QPushButton("Stop")
        self.clear_btn = QPushButton("Clear")
        row.addWidget(self.inp, 1)
        row.addWidget(self.stop_btn)
        row.addWidget(self.clear_btn)
        lay.addWidget(self.out)
        lay.addLayout(row)
        self.proc = QProcess(self)
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PATH", str(Path(sys.executable).parent) + os.pathsep + env.value("PATH"))
        env.insert("PYTHONUNBUFFERED", "1")
        self.proc.setProcessEnvironment(env)
        self.proc.setWorkingDirectory(str(ROOT))
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(lambda code_, _s: self.out.appendPlainText(f"[exit {code_}]"))
        self.inp.returnPressed.connect(self.run)
        self.stop_btn.clicked.connect(self.proc.kill)
        self.clear_btn.clicked.connect(self.out.clear)
        self.history, self.hpos = [], 0

    def run(self):
        cmd = self.inp.text().strip()
        if not cmd:
            return
        if self.proc.state() != QProcess.NotRunning:
            self.out.appendPlainText("[busy: stop the running command first]")
            return
        self.inp.clear()
        self.history.append(cmd)
        self.hpos = len(self.history)
        self.out.appendPlainText(f"$ {cmd}")
        self.session.log(f"$ {cmd}", "terminal")
        self.proc.start("/bin/zsh", ["-c", cmd])

    def _read(self):
        self.out.appendPlainText(bytes(self.proc.readAllStandardOutput()).decode(errors="replace").rstrip())

    def keyPressEvent(self, e):
        if e.key() in (Qt.Key_Up, Qt.Key_Down) and self.history:
            self.hpos = max(0, min(len(self.history), self.hpos + (-1 if e.key() == Qt.Key_Up else 1)))
            self.inp.setText(self.history[self.hpos] if self.hpos < len(self.history) else "")
        else:
            super().keyPressEvent(e)


class StudioWindow(QMainWindow):
    def __init__(self, session: StudioSession, api_server=None):
        super().__init__()
        self.session, self.api_server = session, api_server
        self.api = StudioAPI(session)
        self.setWindowTitle("ZULF Studio")
        self.resize(1500, 950)
        self.bridge = Bridge()
        self.bridge.changed.connect(self.on_changed)
        self.bridge.logged.connect(self.on_logged)
        session.listeners.append(self.bridge.changed.emit)
        session.log_listeners.append(self.bridge.logged.emit)
        self.coupling_rows = {}
        self.redraw_timer = QTimer(self, singleShot=True, interval=25)
        self.redraw_timer.timeout.connect(self.redraw)
        self.status_timer = QTimer(self, interval=1000)
        self.status_timer.timeout.connect(self.update_fit_status)
        self.status_timer.start()

        left = QWidget()
        self.left_lay = QVBoxLayout(left)
        self.left_lay.addWidget(self._structure_box())
        self.coupling_box = QGroupBox("Couplings (Hz)")
        self.coupling_lay = QVBoxLayout(self.coupling_box)
        self.left_lay.addWidget(self.coupling_box)
        self.left_lay.addWidget(self._field_box())
        self.left_lay.addWidget(self._data_box())
        self.left_lay.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidget(left)
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(470)

        center = QWidget()
        cl = QVBoxLayout(center)
        bar = QHBoxLayout()
        self.view_lo = QDoubleSpinBox(decimals=2, maximum=5000.0, keyboardTracking=False)
        self.view_hi = QDoubleSpinBox(decimals=2, maximum=5000.0, keyboardTracking=False)
        self.part = QComboBox()
        self.part.addItems(["re", "im", "abs"])
        self.show_sticks = QCheckBox("lines", checked=True)
        self.show_trace = QCheckBox("fit trace", checked=True)
        self.lock_scale = QCheckBox("lock scale")
        self.lock_scale.setToolTip("freeze the display scale of the simulation (automatic: least squares on |spectrum|)")
        self.lock_scale.toggled.connect(lambda on: self._guard(self.session.lock_scale, on))
        self.field_label = QLabel()
        self.field_label.setStyleSheet("font-weight: 600")
        for w in (QLabel("view"), self.view_lo, QLabel("-"), self.view_hi, QLabel("Hz   part"), self.part,
                  self.show_sticks, self.show_trace, self.lock_scale):
            bar.addWidget(w)
        bar.addStretch(1)
        bar.addWidget(self.field_label)
        self.fig = Figure(figsize=(8, 6), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        cl.addLayout(bar)
        cl.addWidget(NavigationToolbar2QT(self.canvas, self))
        cl.addWidget(self.canvas, 1)
        self.view_lo.valueChanged.connect(self._view_from_spins)
        self.view_hi.valueChanged.connect(self._view_from_spins)
        self.part.currentTextChanged.connect(lambda p: self.session.set_display(part=p))
        self.show_sticks.toggled.connect(lambda _: self.schedule())
        self.show_trace.toggled.connect(lambda _: self.schedule())

        self.tabs = QTabWidget()
        self.tabs.addTab(self._lines_tab(), "Lines")
        self.tabs.addTab(self._fit_tab(), "Fit")
        self.tabs.addTab(self._figure_tab(), "Figure")
        self.log_view = QPlainTextEdit(readOnly=True)
        self.log_view.setFont(MONO)
        self.log_view.setMaximumBlockCount(5000)
        self.tabs.addTab(self.log_view, "Log")
        self.tabs.addTab(Terminal(session), "Terminal")
        self.tabs.addTab(Console({"session": session, "api": self.api, "np": np}), "Python")
        self.tabs.addTab(self._api_tab(), "AI API")

        right = QSplitter(Qt.Vertical)
        right.addWidget(center)
        right.addWidget(self.tabs)
        right.setSizes([620, 330])
        main = QSplitter(Qt.Horizontal)
        main.addWidget(scroll)
        main.addWidget(right)
        main.setSizes([480, 1020])
        self.setCentralWidget(main)
        self._menu()
        for e in session.read_log(200):
            self.on_logged(e)
        self.rebuild_couplings()
        self.refresh_widgets()
        self.schedule()

    # ---- left panels -----------------------------------------------------------------------
    def _structure_box(self):
        box = QGroupBox("Structure")
        lay = QVBoxLayout(box)
        row = QHBoxLayout()
        self.motif = QComboBox()
        self.motif.addItem("(custom JSON below)")
        for m in self.session.list_motifs():
            self.motif.addItem(m["name"])
        self.exchange = QComboBox()
        self.exchange.addItems(["fast", "slow"])
        self.exchange.setToolTip("exchangeable N-H / O-H protons: fast = decoupled, slow = kept")
        row.addWidget(QLabel("motif"))
        row.addWidget(self.motif, 1)
        row.addWidget(QLabel("N-H/O-H"))
        row.addWidget(self.exchange)
        self.spec_edit = QPlainTextEdit()
        self.spec_edit.setFont(MONO)
        self.spec_edit.setFixedHeight(92)
        build = QPushButton("Build")
        build.clicked.connect(self.build_structure)
        self.motif.currentTextChanged.connect(self._motif_chosen)
        self.exchange.currentTextChanged.connect(lambda m: self._guard(self.session.set_exchange, m))
        self.components_label = QLabel()
        self.components_label.setWordWrap(True)
        lay.addLayout(row)
        lay.addWidget(self.spec_edit)
        r2 = QHBoxLayout()
        r2.addWidget(self.components_label, 1)
        r2.addWidget(build)
        lay.addLayout(r2)
        return box

    def _motif_chosen(self, name):
        m = next((m for m in self.session.list_motifs() if m["name"] == name), None)
        if m is not None:
            self.spec_edit.setPlainText(json.dumps({"motif": name, "one_bond": m["one_bond"]}, indent=1))

    def build_structure(self):
        try:
            spec = json.loads(self.spec_edit.toPlainText())
        except json.JSONDecodeError as exc:
            QMessageBox.warning(self, "Structure", f"Invalid JSON: {exc}")
            return
        self._guard(self.session.set_structure, spec)

    def _field_box(self):
        box = QGroupBox("Static field and line width")
        lay = QVBoxLayout(box)
        s = self.session
        self.b_t = ValueSlider("B transverse", 0.0, 300.0, s.field_nt[0], 1, "nT", fine=5.0)
        self.b_z = ValueSlider("B z", 0.0, 300.0, s.field_nt[1], 1, "nT", fine=5.0)
        self.rate = ValueSlider("decay rate", 0.05, 30.0, s.rate_per_s, 3, "1/s", log=True)
        self.b_mag = QLabel()
        self.b_t.value_changed.connect(lambda v: self.session.set_field(transverse_nt=v))
        self.b_z.value_changed.connect(lambda v: self.session.set_field(z_nt=v))
        self.rate.value_changed.connect(self.session.set_linewidth)
        for w in (self.b_t, self.b_z, self.b_mag, self.rate):
            lay.addWidget(w)
        zero = QPushButton("Zero field")
        zero.clicked.connect(lambda: self.session.set_field(0.0, 0.0))
        lay.addWidget(zero)
        return box

    def _data_box(self):
        box = QGroupBox("Experimental spectrum")
        lay = QVBoxLayout(box)
        row = QHBoxLayout()
        self.data_label = QLabel("none")
        load = QPushButton("Load series.json ...")
        load.clicked.connect(self.load_series)
        row.addWidget(self.data_label, 1)
        row.addWidget(load)
        lay.addLayout(row)
        self.phase = ValueSlider("display phase", -180.0, 180.0, 0.0, 1, "deg")
        self.delay = ValueSlider("display delay", -10.0, 10.0, 0.0, 3, "ms", fine=0.2)
        self.phase.value_changed.connect(lambda v: self.session.set_display(phase_deg=v))
        self.delay.value_changed.connect(lambda v: self.session.set_display(delay_ms=v))
        lay.addWidget(self.phase)
        lay.addWidget(self.delay)
        hint = QLabel("The display phase only rotates the shown data; the simulation is a quick Lorentzian look. "
                      "Use Fit for the model rendered through the data processing.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: gray; font-size: 10px")
        lay.addWidget(hint)
        return box

    def rebuild_couplings(self):
        while self.coupling_lay.count():
            w = self.coupling_lay.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        self.coupling_rows = {}
        for c in self.session.couplings():
            v = c["value"]
            span = 3.0 if c["one_bond"] else 6.0
            row = ValueSlider(c["key"] + (" *" if c["unspecified"] else ""), v - span, v + span, v, 3, "Hz",
                              fine=0.3)
            row.setToolTip("* not specified by the structure (built as 0 Hz)" if c["unspecified"] else "")
            row.value_changed.connect(lambda val, k=c["key"]: self._guard(self.session.set_couplings, {k: val}))
            self.coupling_lay.addWidget(row)
            self.coupling_rows[c["key"]] = row
        add = QWidget()
        al = QHBoxLayout(add)
        al.setContentsMargins(0, 4, 0, 0)
        self.new_key = QLineEdit(placeholderText="J(a,b)")
        self.new_val = QDoubleSpinBox(decimals=3, minimum=-500.0, maximum=500.0)
        btn = QPushButton("Add / set")
        btn.clicked.connect(lambda: self._guard(self.session.set_couplings,
                                                {self.new_key.text().strip(): self.new_val.value()}))
        al.addWidget(self.new_key, 1)
        al.addWidget(self.new_val)
        al.addWidget(btn)
        self.coupling_lay.addWidget(add)
        comps = self.session.components()
        self.components_label.setText("isotopologues: " + ", ".join(
            f"<span style='color:{COMPONENT_COLORS[i % 6]}'>{c['label']}</span>" for i, c in enumerate(comps)))

    # ---- tabs ------------------------------------------------------------------------------
    def _lines_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        row = QHBoxLayout()
        self.min_rel = QDoubleSpinBox(decimals=3, minimum=0.0, maximum=1.0, value=0.02, singleStep=0.01)
        self.min_rel.valueChanged.connect(lambda _: self.fill_lines())
        export = QPushButton("Export (parameters, lines, spectrum, figure) ...")
        export.clicked.connect(self.export)
        row.addWidget(QLabel("min relative amplitude"))
        row.addWidget(self.min_rel)
        row.addStretch(1)
        row.addWidget(export)
        self.lines_table = QTableWidget(0, 4)
        self.lines_table.setHorizontalHeaderLabels(["isotopologue", "frequency (Hz)", "amplitude", "relative"])
        self.lines_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.lines_table.setSortingEnabled(True)
        lay.addLayout(row)
        lay.addWidget(self.lines_table)
        return w

    def _fit_tab(self):
        w = QWidget()
        lay = QHBoxLayout(w)
        form = QFormLayout()
        self.f_starts = QSpinBox(minimum=1, maximum=200, value=8)
        self.f_workers = QSpinBox(minimum=1, maximum=os.cpu_count() or 8, value=min(4, os.cpu_count() or 4))
        self.f_nfev = QSpinBox(minimum=10, maximum=5000, value=300, singleStep=50)
        self.f_trace = QSpinBox(minimum=0, maximum=400, value=40)
        self.f_field = QCheckBox("fit the static field (starts from the sliders)", checked=True)
        self.f_edges = QLineEdit(placeholderText="e.g. 135.5,137.2,200 or auto (empty: one rate per isotopologue)")
        self.f_rates = QLineEdit("0.2,15")
        self.f_extra = QLineEdit(placeholderText="extra fit_joint_series options, e.g. --model-line-passes 1")
        for label, wid in (("starts", self.f_starts), ("workers", self.f_workers), ("max evaluations", self.f_nfev),
                           ("trace frames", self.f_trace), ("", self.f_field), ("rate families", self.f_edges),
                           ("rate bounds (1/s)", self.f_rates), ("extra", self.f_extra)):
            form.addRow(label, wid)
        btns = QHBoxLayout()
        self.f_start = QPushButton("Start fit")
        self.f_start.setStyleSheet("font-weight: 600")
        self.f_stop = QPushButton("Stop")
        self.f_apply = QPushButton("Apply result")
        self.f_load = QPushButton("Load run ...")
        for b in (self.f_start, self.f_stop, self.f_apply, self.f_load):
            btns.addWidget(b)
        self.f_start.clicked.connect(self.start_fit)
        self.f_stop.clicked.connect(self.session.stop_fit)
        self.f_apply.clicked.connect(lambda: self._guard(self.session.apply_fit))
        self.f_load.clicked.connect(self.load_run)
        left = QVBoxLayout()
        left.addLayout(form)
        left.addLayout(btns)
        self.f_status = QLabel("no fit")
        self.f_status.setWordWrap(True)
        left.addWidget(self.f_status)
        left.addStretch(1)
        right = QVBoxLayout()
        right.addWidget(QLabel("<b>Fit progress</b> (trace of the fit: drag through the evaluations)"))
        self.trace_slider = QSlider(Qt.Horizontal)
        self.trace_slider.setEnabled(False)
        self.trace_slider.valueChanged.connect(lambda i: self._guard(self.session.trace_frame, i))
        self.trace_info = QPlainTextEdit(readOnly=True)
        self.trace_info.setFont(MONO)
        apply_frame = QPushButton("Copy this frame's couplings to the sliders")
        apply_frame.clicked.connect(lambda: self._guard(self.session.trace_frame, self.trace_slider.value(), True))
        right.addWidget(self.trace_slider)
        right.addWidget(self.trace_info, 1)
        right.addWidget(apply_frame)
        lay.addLayout(left, 1)
        lay.addLayout(right, 1)
        return w

    def _figure_tab(self):
        w = QWidget()
        lay = QHBoxLayout(w)
        form = QFormLayout()
        self.g_title = QLineEdit(placeholderText="title (default: the spectrum id)")
        self.g_wide = QLineEdit("5,300")
        self.g_segments = QLineEdit(placeholderText="lo,hi;lo,hi (default: the fit ranges)")
        self.g_gains = QLineEdit(placeholderText="e.g. 1,2.5")
        self.g_window = QDoubleSpinBox(decimals=2, minimum=0.0, maximum=5.0, value=0.1, singleStep=0.05)
        self.g_dpi = QSpinBox(minimum=72, maximum=1200, value=300, singleStep=50)
        fmt = QHBoxLayout()
        fmt.setContentsMargins(0, 0, 0, 0)
        self.g_formats = {f: QCheckBox(f, checked=True) for f in ("png", "pdf", "svg")}
        for c in self.g_formats.values():
            fmt.addWidget(c)
        fmt_w = QWidget()
        fmt_w.setLayout(fmt)
        self.g_colors = QLineEdit(placeholderText='{"13C@C1": "#2f8f5b"}')
        self.g_insets = QLineEdit(placeholderText="insets JSON (paper_figure --insets)")
        try:
            import rdkit  # noqa: F401
        except ImportError:
            self.g_insets.setEnabled(False)
            self.g_insets.setPlaceholderText("structure insets need RDKit (not installed)")
        for label, wid in (("title", self.g_title), ("whole range (Hz)", self.g_wide), ("detail panels", self.g_segments),
                           ("panel gains", self.g_gains), ("display window (1/s)", self.g_window), ("formats", fmt_w),
                           ("PNG dpi", self.g_dpi), ("colours", self.g_colors), ("insets", self.g_insets)):
            form.addRow(label, wid)
        btns = QHBoxLayout()
        self.g_make = QPushButton("Generate figure")
        self.g_make.setStyleSheet("font-weight: 600")
        self.g_export = QPushButton("Export ...")
        self.g_open = QPushButton("Open folder")
        for b in (self.g_make, self.g_export, self.g_open):
            btns.addWidget(b)
        self.g_make.clicked.connect(self.make_figure)
        self.g_export.clicked.connect(self.export_figure)
        self.g_open.clicked.connect(lambda: self.session.figure and QDesktopServices.openUrl(
            QUrl.fromLocalFile(self.session.figure["directory"])))
        self.g_status = QLabel("Draws the applied fit; after moving sliders, the current parameters, labelled "
                               "'manual parameters (not a fit)'.")
        self.g_status.setWordWrap(True)
        left = QVBoxLayout()
        left.addLayout(form)
        left.addLayout(btns)
        left.addWidget(self.g_status)
        left.addStretch(1)
        self.g_preview = QLabel("no figure yet")
        self.g_preview.setAlignment(Qt.AlignCenter)
        self.g_preview.setMinimumSize(200, 120)
        scroll = QScrollArea()
        scroll.setWidget(self.g_preview)
        scroll.setWidgetResizable(True)
        lay.addLayout(left, 1)
        lay.addWidget(scroll, 2)
        return w

    def make_figure(self):
        opts = {"title": self.g_title.text().strip(), "wide": self.g_wide.text().strip(),
                "segments": self.g_segments.text().strip(), "gains": self.g_gains.text().strip(),
                "display_window": self.g_window.value(), "dpi": self.g_dpi.value(),
                "formats": ",".join(f for f, c in self.g_formats.items() if c.isChecked()) or "png",
                "colors": self.g_colors.text().strip(), "insets": self.g_insets.text().strip()}
        if self.session.data is not None and not self.session.data.get("source_fid"):
            fid, _ = QFileDialog.getOpenFileName(self, "Averaged FID of this spectrum (not in the series file)",
                                                 str(Path.home() / "research"), "FID (*.npy)")
            if not fid:
                return
            self.session.data["source_fid"] = fid
        if self._guard(self.session.make_figure, **opts) is not None:
            self.g_status.setText("generating ...")
            self.g_make.setEnabled(False)

    def export_figure(self):
        st = self.session.figure_status()
        if not st.get("files"):
            self.statusBar().showMessage("no figure yet", 5000)
            return
        d = QFileDialog.getExistingDirectory(self, "Export figure to", st["directory"])
        if d:
            r = self._guard(self.session.export_figure, d)
            if r:
                self.statusBar().showMessage(f"figure exported to {r['directory']}: {', '.join(r['files'])}", 8000)

    def show_figure(self):
        st = self.session.figure_status()
        self.g_make.setEnabled(not st.get("running"))
        if st.get("returncode") not in (0, None):
            self.g_status.setText(f"figure failed (code {st['returncode']}); see the Log tab")
            return
        if st.get("png"):
            pix = QPixmap(st["png"])
            self.g_preview.setPixmap(pix.scaledToWidth(max(self.g_preview.width() - 20, 400),
                                                       Qt.SmoothTransformation))
            self.g_status.setText(("MANUAL PARAMETERS (not a fit). " if st["manual"] else "Applied fit. ")
                                  + f"{st['directory']}\n" + ", ".join(st["files"]))

    def _api_tab(self):
        w = QPlainTextEdit(readOnly=True)
        w.setFont(MONO)
        if self.api_server is None:
            w.setPlainText("API server off (start with --api-port PORT).")
            return w
        host, port = self.api_server.server_address[:2]
        base = f"http://{host}:{port}"
        lines = [f"AI API on {base}", "",
                 f"  GET  {base}/api/tools            (?format=anthropic or openai for tool schemas)",
                 f"  GET  {base}/api/state",
                 f"  POST {base}/api/call   {{\"tool\": NAME, \"args\": {{...}}}}",
                 f"  POST {base}/api/NAME   {{...args...}}", "",
                 "Example:", f"  curl -s -X POST {base}/api/set_field -d '{{\"transverse_nt\": 37, \"z_nt\": 45}}'", "",
                 "Tools:"]
        lines += [f"  {n:16s} {d}" for n, (_, d, _) in TOOLS.items()]
        w.setPlainText("\n".join(lines))
        return w

    def _menu(self):
        m = self.menuBar().addMenu("&File")
        for text, key, fn in (("Load series ...", "Ctrl+O", self.load_series), ("Load fit run ...", "Ctrl+R", self.load_run),
                              ("Export ...", "Ctrl+E", self.export), ("Generate figure", "Ctrl+G", self.make_figure)):
            a = QAction(text, self)
            a.setShortcut(QKeySequence(key))
            a.triggered.connect(fn)
            m.addAction(a)

    # ---- actions ---------------------------------------------------------------------------
    def _guard(self, fn, *args, **kw):
        try:
            return fn(*args, **kw)
        except Exception as exc:
            self.session.log(f"error: {type(exc).__name__}: {exc}", "session")
            self.statusBar().showMessage(f"{type(exc).__name__}: {exc}", 8000)

    def load_series(self):
        path, _ = QFileDialog.getOpenFileName(self, "Series file", str(ROOT / "runs" / "series"), "series (*.json)")
        if path:
            self._guard(self.session.load_spectrum, series=path)

    def load_run(self):
        d = QFileDialog.getExistingDirectory(self, "Fit run directory", str(ROOT / "runs"))
        if d:
            if (Path(d) / "fit.json").exists():
                self._guard(self.session.apply_fit, d)
            elif (Path(d) / "trace.json").exists():
                self._guard(self.session.load_trace, d)

    def start_fit(self):
        extra = self.f_extra.text().split()
        self._guard(self.session.start_fit, starts=self.f_starts.value(), workers=self.f_workers.value(),
                    max_nfev=self.f_nfev.value(), trace=self.f_trace.value(), fit_field=self.f_field.isChecked(),
                    family_edges=self.f_edges.text().strip(), rate_bounds=self.f_rates.text().strip() or "0.2,15",
                    extra_args=extra)
        self.tabs.setCurrentIndex(1)

    def export(self):
        d = QFileDialog.getExistingDirectory(self, "Export to", str(self.session.workspace))
        if d:
            r = self._guard(self.session.export, d)
            if r:
                self.fig.savefig(Path(d) / "figure.png", dpi=160)
                self.statusBar().showMessage(f"exported to {d}", 6000)

    def _view_from_spins(self):
        lo, hi = self.view_lo.value(), self.view_hi.value()
        if hi > lo:
            self._guard(self.session.set_view, lo, hi)

    # ---- session events ----------------------------------------------------------------------
    def on_changed(self, event):
        if event == "structure":
            self.rebuild_couplings()
        elif event == "couplings" and set(self.coupling_rows) != {c["key"] for c in self.session.couplings()}:
            self.rebuild_couplings()
        self.refresh_widgets()
        if event in ("fit_done", "trace", "fit_started"):
            self.update_fit_status()
        if event == "figure_done":
            self.show_figure()
        self.schedule()

    def on_logged(self, e):
        self.log_view.appendPlainText(f"{e['time']} [{e['source']}] {e['message']}")

    def refresh_widgets(self):
        s = self.session
        for c in s.couplings():
            row = self.coupling_rows.get(c["key"])
            if row is not None and abs(row.value() - c["value"]) > 1e-9:
                row.set_value(c["value"])
        self.b_t.set_value(s.field_nt[0])
        self.b_z.set_value(s.field_nt[1])
        self.rate.set_value(s.rate_per_s)
        bt, bz = s.field_nt
        txt = "zero field" if bt == 0 and bz == 0 else \
            f"B\u22a5 {bt:.1f} nT   Bz {bz:.1f} nT   |B| {math.hypot(bt, bz):.1f} nT"
        self.b_mag.setText(f"|B| = {math.hypot(bt, bz):.1f} nT   (FWHM {s.rate_per_s / math.pi:.3f} Hz)")
        self.field_label.setText(txt)
        for spin, v in ((self.view_lo, s.view[0]), (self.view_hi, s.view[1])):
            spin.blockSignals(True)
            spin.setValue(v)
            spin.blockSignals(False)
        self.part.blockSignals(True)
        self.part.setCurrentText(s.display)
        self.part.blockSignals(False)
        self.lock_scale.blockSignals(True)
        self.lock_scale.setChecked(s.scale_lock is not None)
        self.lock_scale.blockSignals(False)
        self.data_label.setText(s.data["label"] if s.data else "none")
        self.exchange.blockSignals(True)
        self.exchange.setCurrentText(s.exchange)
        self.exchange.blockSignals(False)
        if not self.spec_edit.hasFocus():
            self.spec_edit.setPlainText(json.dumps(s.spec, indent=1))
        if s.trace is not None:
            n = len(s.trace["meta"]["frames"])
            self.trace_slider.blockSignals(True)
            self.trace_slider.setEnabled(True)
            self.trace_slider.setRange(0, n - 1)
            self.trace_slider.setValue(max(s.trace_index, 0))
            self.trace_slider.blockSignals(False)
            fr = s.trace["meta"]["frames"][max(s.trace_index, 0)]
            self.trace_info.setPlainText(
                f"run {s.trace['run']}\nframe {s.trace_index + 1}/{n}  evaluation {fr['evaluation']}  "
                f"stage {fr['stage']}\nobjective {fr['objective']:.6g}\n\n" +
                "\n".join(f"{k:14s} {v[0]:10.4f} Hz" for k, v in fr["J"].items()))
        else:
            self.trace_slider.setEnabled(False)

    def update_fit_status(self):
        st = self.session.fit_status()
        if st.get("out_dir") is None:
            return
        best = st.get("best_objective")
        self.f_status.setText(
            f"{'running' if st['running'] else 'finished (code %s)' % st['returncode']}: {st['seconds']} s, "
            f"{st['starts_finished']} starts finished, best objective {best if best is None else f'{best:.5g}'}\n"
            f"{st['out_dir']}\n{st['last_line'][:160]}")
        self.f_start.setEnabled(not st["running"])

    def fill_lines(self, lines=None):
        lines = lines if lines is not None else self.session.lines(self.min_rel.value(), self.session.view)
        lines = [r for r in lines if r["relative"] >= self.min_rel.value()]
        t = self.lines_table
        t.setSortingEnabled(False)
        t.setRowCount(len(lines))
        for i, r in enumerate(lines):
            for j, v in enumerate((r["component"], f"{r['frequency_hz']:.4f}", f"{r['amplitude']:.4g}",
                                   f"{r['relative']:.4f}")):
                item = QTableWidgetItem(v)
                if j:
                    item.setData(Qt.DisplayRole, float(v))
                t.setItem(i, j, item)
        t.setSortingEnabled(True)

    # ---- plot --------------------------------------------------------------------------------
    def schedule(self):
        self.redraw_timer.start()

    def redraw(self):
        s = self.session
        try:
            sim = s.simulate()
        except Exception as exc:
            self.statusBar().showMessage(f"simulation: {exc}", 8000)
            return
        part = sim["part"]
        pick = (lambda re, im: np.asarray(re) if part == "re" else np.asarray(im) if part == "im"
                else np.hypot(re, im))
        f = np.asarray(sim["f"])
        self.fig.clear()
        has_data = "data_re" in sim
        rows = [3, 1, 1] if has_data else [3, 1]
        axes = self.fig.subplots(len(rows), 1, sharex=True, gridspec_kw={"height_ratios": rows})
        ax = axes[0]
        if has_data:
            d = pick(sim["data_re"], sim["data_im"])
            ax.plot(f, d, color="black", lw=0.8, label=s.data["label"])
        m = pick(sim["sim_re"], sim["sim_im"])
        ax.plot(f, m, color="#d1495b", lw=1.0, label="simulation (Lorentzian)")
        tr = s.trace if self.show_trace.isChecked() else None
        if tr is not None and s.trace_index >= 0:
            ft = tr["f"]
            sel = (ft >= s.view[0]) & (ft <= s.view[1])
            ph = np.exp(1j * (np.radians(s.data_phase_deg) + 2 * np.pi * ft * 1e-3 * s.data_delay_ms))
            mod = tr["models"][s.trace_index] * ph
            mt = mod.real if part == "re" else mod.imag if part == "im" else np.abs(mod)
            ax.plot(ft[sel], mt[sel], color="#2a9d5c", lw=1.0,
                    label=f"fit trace frame {s.trace_index + 1} (objective "
                          f"{tr['meta']['frames'][s.trace_index]['objective']:.4g})")
        ax.legend(loc="upper left", fontsize=8, frameon=False)
        bt, bz = s.field_nt
        ax.set_title(f"{s.spec.get('compound', s.spec.get('motif', ''))}:  "
                     + ("zero field" if bt == bz == 0 else f"$B_\\perp$ {bt:.1f} nT, $B_z$ {bz:.1f} nT")
                     + f",  rate {s.rate_per_s:.3g} 1/s", fontsize=10, loc="left")
        if has_data:
            axes[1].plot(f, d - m, color="gray", lw=0.7)
            axes[1].axhline(0, color="#cccccc", lw=0.6)
            axes[1].set_ylabel("residual", fontsize=8)
        sax = axes[-1]
        if self.show_sticks.isChecked():
            labels = [c["label"] for c in s.components()]
            segs, cols = [], []
            for r in sim["lines"]:
                segs.append([(r["frequency_hz"], 0), (r["frequency_hz"], r["relative"])])
                cols.append(COMPONENT_COLORS[labels.index(r["component"]) % 6] if r["component"] in labels else "k")
            sax.add_collection(LineCollection(segs, colors=cols, linewidths=1.2))
            sax.set_ylim(0, 1.05)
        sax.set_ylabel("lines", fontsize=8)
        sax.set_xlim(*s.view)
        sax.set_xlabel("frequency (Hz)")
        for a in axes:
            a.tick_params(labelsize=8)
        self.canvas.draw_idle()
        self.fill_lines(sim["lines"])
        msg = (f"simulation {1e3 * sim['seconds']:.0f} ms, {len(sim['lines'])} lines in view, scale {sim['scale']:.4g}"
               + (" (locked)" if s.scale_lock is not None else ""))
        if has_data:
            msg += f", rms residual {sim['residual_rms']:.3g}"
        self.statusBar().showMessage(msg)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="ZULF Studio: real-time zero-field spectrum simulator with fitting")
    ap.add_argument("--structure", default="", help="structure JSON (motif or chain); default acetonitrile")
    ap.add_argument("--series", default="", help="series.json of a processed spectrum (scripts/make_series_entry.py)")
    ap.add_argument("--fit", default="", help="fit run directory to apply (couplings, field, rate, trace)")
    ap.add_argument("--api-port", type=int, default=8766, help="AI API port (0: any free port, -1: off)")
    ap.add_argument("--workspace", default="runs/studio")
    ap.add_argument("--no-gui", action="store_true", help="API server only (for agents; Ctrl+C to quit)")
    args = ap.parse_args(argv)
    session = StudioSession(json.loads(args.structure) if args.structure else None, workspace=args.workspace)
    if args.series:
        session.load_spectrum(series=args.series)
    if args.fit:
        session.apply_fit(args.fit)
    server = None if args.api_port < 0 else serve(session, port=args.api_port, background=not args.no_gui)
    if args.no_gui:
        print(f"ZULF Studio API on http://127.0.0.1:{server.server_address[1]}/api/tools", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
        return
    app = QApplication.instance() or QApplication(sys.argv)
    win = StudioWindow(session, server)
    win.show()
    sys.exit(app.exec())
