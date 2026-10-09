"""ZULF Studio window (PySide6): live sliders, plot, line table, fitting with a progress slider, log, terminal,
Python console and the AI API, all on one StudioSession.

    python scripts/run_studio.py [--structure JSON] [--series SERIES.json] [--fit RUN_DIR] [--api-port 8766]
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

for _var in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):     # one BLAS thread: the
    os.environ.setdefault(_var, "1")                  # window should not compete with the fits it starts

import numpy as np
from PySide6.QtCore import QObject, QProcess, QProcessEnvironment, QSettings, Qt, QTimer, Signal
from PySide6.QtCore import QUrl
from PySide6.QtGui import QAction, QDesktopServices, QFont, QKeySequence, QPixmap
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPlainTextEdit, QPushButton, QScrollArea, QSlider, QSpinBox, QSplitter,
                               QTableWidget, QTableWidgetItem, QTabWidget, QTextBrowser, QVBoxLayout, QWidget)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from .api import TOOLS, StudioAPI, serve  # noqa: E402
from .jobs_ui import ActivityIndicator, AnalysisPanel, JobsPanel, MachineLabel, confirm_workers, job_line  # noqa: E402
from .session import SESSION_SUFFIX  # noqa: E402
from .theme import DARK, ISOTOPOLOGUE, LIGHT, matplotlib_style, stylesheet  # noqa: E402
from .session import ROOT, StudioSession  # noqa: E402

COMPONENT_COLORS = ISOTOPOLOGUE
SPIN_TEXT = "\u25d0\u25d3\u25d1\u25d2"

AI_SETUP_GUIDE = """
<h3>Setting up the AI assistant</h3>
<p>The assistant sends your request to a language model, which then operates this session with the Studio
tools (every call is shown here and in the Log tab). Using it needs API access, billed per token and separate
from a Claude or ChatGPT subscription. Ask your lab or university first: an organisation may already have an
API account.</p>
<h4>1. Get a key</h4>
<ul>
<li><b>Anthropic (Claude)</b>: platform.claude.com - sign in, add a payment method under Billing, create a key
under API Keys (it starts with <code>sk-ant-</code>). Default model <code>claude-opus-5-5</code>
($4 / $20 per million input / output tokens; a session of a few dozen tool calls costs about 1-3 USD).</li>
<li><b>OpenAI (e.g. Codex)</b>: create an API key in the OpenAI platform and pick a model your account offers;
type its name in the model field (or set <code>OPENAI_MODEL</code>).</li>
</ul>
<h4>2. Give the key to Studio (one of these)</h4>
<ul>
<li><b>Easiest</b>: paste it into the key field below, tick <b>Remember (macOS Keychain)</b> and press
<b>Use</b>. Studio stores it encrypted in the macOS Keychain and loads it every time it starts: enter it once.
Without Remember it stays in memory until Studio closes. Never written to a file or the log; Forget removes it
(also from the Keychain).</li>
<li><b>Before starting Studio</b>, in iTerm2:<br>
<code>export ANTHROPIC_API_KEY=sk-ant-...</code> (or <code>export OPENAI_API_KEY=...</code>), then start Studio
from that same terminal: <code>python scripts/run_studio.py ...</code></li>
<li><b>Anthropic login profile</b>: <code>ant auth login</code> stores a profile that the Claude library reads
on its own (also while Studio is open). It needs the <code>ant</code> command-line tool, which is not installed
on this Mac.</li>
</ul>
<p><b>Not</b> in the Terminal tab of Studio: each command there runs in its own shell, so an
<code>export</code> there does not reach Studio.</p>
<h4>3. Use it</h4>
<p>Choose the provider and model, write a request (Ctrl+Enter sends), for example
<i>"auto-phase the data, then scan J(C1,HC1) from 136.2 to 136.4 Hz and set the value with the smallest rms
residual"</i>. Stop ends the loop after the current step; New conversation forgets the history; "max steps" limits
the tool calls per request.</p>
<p>Never put a key into code, configuration files or the repository (it is public). Full guide:
docs/STUDIO.md.</p>
"""
MONO = QFont("Menlo", 11)


class Bridge(QObject):
    """Session events (any thread: API server, fit reader, AI assistant) -> Qt signals handled in the GUI thread."""
    changed = Signal(str)
    logged = Signal(dict)
    ai_message = Signal(str, str)


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
        self.name.setMinimumWidth(96)
        self.spin = QDoubleSpinBox()
        self.spin.setDecimals(decimals)
        self.spin.setRange(min(lo, value), max(hi, value))
        self.spin.setSingleStep(10 ** -min(decimals, 2))
        self.spin.setKeyboardTracking(False)
        self.spin.setSuffix(f" {unit}" if unit else "")
        self.spin.setMinimumWidth(104)
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
            fl.setObjectName("fine")
            self.fine.setObjectName("fineSlider")
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
        self.out.setObjectName("mono")
        self.inp = QLineEdit()
        self.inp.setObjectName("mono")
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
        self.out.setObjectName("mono")
        row = QHBoxLayout()
        self.inp = QLineEdit()
        self.inp.setObjectName("mono")
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
        self.theme_mode = "auto"
        self.t = LIGHT
        self.apply_theme()
        self.bridge = Bridge()
        self.bridge.changed.connect(self.on_changed)
        self.bridge.logged.connect(self.on_logged)
        self.bridge.ai_message.connect(self.on_ai_message)
        self.assistant = None
        session.listeners.append(self.bridge.changed.emit)
        session.log_listeners.append(self.bridge.logged.emit)
        self.coupling_rows = {}
        self.dirty = False
        self.redraw_timer = QTimer(self, singleShot=True, interval=25)
        self.redraw_timer.timeout.connect(self.redraw)
        self.status_timer = QTimer(self, interval=1000)
        self.status_timer.timeout.connect(self.update_fit_status)
        self.status_timer.timeout.connect(self.update_jobs)
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
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(500)

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
        self.field_label.setObjectName("badge")
        for w in (QLabel("View"), self.view_lo, QLabel("to"), self.view_hi, QLabel("Hz"), QLabel("   Part"), self.part,
                  QLabel("  "), self.show_sticks, self.show_trace, self.lock_scale):
            bar.addWidget(w)
        bar.addStretch(1)
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
        self.analysis_panel = AnalysisPanel(self.session, self)
        self.tabs.addTab(self.analysis_panel, "Analysis")
        self.jobs_panel = JobsPanel(self.session)
        self.tabs.addTab(self.jobs_panel, "Jobs")
        self.log_view = QPlainTextEdit(readOnly=True)
        self.log_view.setObjectName("mono")
        self.log_view.setMaximumBlockCount(5000)
        self.tabs.addTab(self.log_view, "Log")
        self.tabs.addTab(self._ai_tab(), "AI assistant")
        self.settings = self._settings_dialog()          # AI configuration, API, appearance (menu: Settings)
        self.tools = self._tools_window()                # terminal and Python console (menu: Tools)

        right = QSplitter(Qt.Vertical)
        right.addWidget(center)
        right.addWidget(self.tabs)
        right.setSizes([620, 330])
        main = QSplitter(Qt.Horizontal)
        main.addWidget(scroll)
        main.addWidget(right)
        main.setSizes([510, 990])
        self.splitters = {"right": right, "main": main}
        header = QWidget(objectName="header")
        hl = QHBoxLayout(header)
        hl.setContentsMargins(16, 10, 16, 10)
        title = QLabel("ZULF Studio", objectName="title")
        self.subtitle = QLabel(objectName="subtitle")
        hl.addWidget(title)
        hl.addSpacing(12)
        hl.addWidget(self.subtitle)
        hl.addStretch(1)
        self.busy_pill = QLabel(objectName="pill")
        self.busy_pill.hide()
        hl.addWidget(self.busy_pill)
        hl.addSpacing(8)
        hl.addWidget(self.field_label)
        root = QWidget(objectName="root")
        rl = QVBoxLayout(root)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.setSpacing(0)
        rl.addWidget(header)
        body = QWidget()
        bl = QVBoxLayout(body)
        bl.setContentsMargins(10, 8, 10, 6)
        bl.addWidget(main)
        rl.addWidget(body, 1)
        self.setCentralWidget(root)
        self.activity = ActivityIndicator(self.session)
        self.machine = MachineLabel(self.session)
        self.statusBar().addPermanentWidget(self.machine)
        self.statusBar().addPermanentWidget(self.activity)
        self.settings_store = QSettings("ZULF", "Studio")
        self.setAcceptDrops(True)
        self._menu()
        self._restore_layout()
        self.update_title()
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
        self.spec_edit.setObjectName("mono")
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
        auto = QHBoxLayout()
        self.auto_method = QComboBox()
        self.auto_method.addItem("match the simulation", "model")
        self.auto_method.addItem("data only (peak phases)", "data")
        self.auto_delay = QCheckBox("also delay", checked=True)
        auto_btn = QPushButton("Auto phase")
        auto_btn.setToolTip("model: phase and delay that best match the current simulation (lines roughly in place); "
                            "data: model-free (strongest peaks), delay within +-0.5 ms")
        auto_btn.clicked.connect(self.auto_phase)
        zero = QPushButton("Reset")
        zero.clicked.connect(lambda: self.session.set_display(phase_deg=0.0, delay_ms=0.0))
        for wdg in (auto_btn, self.auto_method, self.auto_delay, zero):
            auto.addWidget(wdg)
        lay.addLayout(auto)
        hint = QLabel("The display phase only rotates the shown data; the simulation is a quick Lorentzian look. "
                      "Use Fit for the model rendered through the data processing.")
        hint.setWordWrap(True)
        hint.setObjectName("hint")
        lay.addWidget(hint)
        return box

    def auto_phase(self):
        r = self._guard(self.session.auto_phase, self.auto_method.currentData(), self.auto_delay.isChecked())
        if r:
            self.statusBar().showMessage(f"auto phase ({r['method']}): {r['phase_deg']:.1f} deg, {r['delay_ms']:.3f} ms, "
                                         f"match {r['match']:.3f}", 8000)

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
        self.f_workers = QSpinBox(minimum=1, maximum=os.cpu_count() or 8,
                                  value=min(4, self.session.machine_status()["suggested_workers"]))
        self.f_workers.setToolTip("worker processes (one start each at a time, one BLAS thread each); the default "
                                  "is what the machine has free (status bar: analysis workers / cores)")
        self.f_nfev = QSpinBox(minimum=10, maximum=5000, value=300, singleStep=50)
        self.f_trace = QSpinBox(minimum=0, maximum=400, value=40)
        self.f_precision = QComboBox()
        for text, value in (("0.1 Hz", 0.1), ("0.01 Hz (default)", 0.01), ("0.001 Hz", 0.001), ("0.0001 Hz", 0.0001),
                            ("run to the tolerances", 0.0)):
            self.f_precision.addItem(text, value)
        self.f_precision.setCurrentIndex(1)
        self.f_precision.setToolTip("the fit stops when every coupling changes by less than this (and the fit no "
                                    "longer improves); also the decimals reported")
        self.f_field = QCheckBox("fit the static field (starts from the sliders)", checked=True)
        self.f_edges = QLineEdit(placeholderText="e.g. 135.5,137.2,200 or auto (empty: one rate per isotopologue)")
        self.f_rates = QLineEdit("0.2,15")
        self.f_extra = QLineEdit(placeholderText="extra fit_joint_series options, e.g. --model-line-passes 1")
        for label, wid in (("starts", self.f_starts), ("workers", self.f_workers), ("max evaluations", self.f_nfev),
                           ("trace frames", self.f_trace), ("coupling precision", self.f_precision), ("", self.f_field), ("rate families", self.f_edges),
                           ("rate bounds (1/s)", self.f_rates), ("extra", self.f_extra)):
            form.addRow(label, wid)
        btns = QHBoxLayout()
        self.f_start = QPushButton("Start fit")
        self.f_start.setObjectName("primary")
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
        self.trace_info.setObjectName("mono")
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
        self.g_make.setObjectName("primary")
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

    def _ai_tab(self):
        """The chat: the request box, the transcript, and one status line; configuration is in Settings."""
        w = QWidget()
        lay = QVBoxLayout(w)
        self._ai_config_widgets()
        row = QHBoxLayout()
        self.ai_status = QLabel()
        settings_btn = QPushButton("Settings ...")
        settings_btn.clicked.connect(lambda: self.open_settings("AI assistant"))
        row.addWidget(self.ai_status, 1)
        row.addWidget(settings_btn)
        lay.addLayout(row)
        self.ai_log = QPlainTextEdit(readOnly=True)
        self.ai_log.setObjectName("mono")
        self.ai_prompt = QPlainTextEdit()
        self.ai_prompt.setFixedHeight(64)
        self.ai_prompt.setPlaceholderText("e.g. auto-phase the data, then scan J(C1,HC1) from 136.2 to 136.4 Hz and set the "
                                          "value with the smallest rms residual  (Ctrl+Enter to send)")
        btns = QHBoxLayout()
        self.ai_send = QPushButton("Send")
        self.ai_send.setObjectName("primary")
        stop = QPushButton("Stop")
        new = QPushButton("New conversation")
        self.ai_send.clicked.connect(self.ai_ask)
        stop.clicked.connect(lambda: self.assistant and self.assistant.stop())
        new.clicked.connect(self._ai_new)
        for b in (self.ai_send, stop, new):
            btns.addWidget(b)
        btns.addStretch(1)
        send_key = QAction(self)
        send_key.setShortcut(QKeySequence("Ctrl+Return"))
        send_key.triggered.connect(self.ai_ask)
        self.ai_prompt.addAction(send_key)
        lay.addWidget(self.ai_log, 1)
        lay.addWidget(self.ai_prompt)
        lay.addLayout(btns)
        self._ai_provider_changed()
        return w

    def _ai_config_widgets(self):
        """Widgets of the AI configuration (shown in the Settings window)."""
        self.ai_provider = QComboBox()
        self.ai_provider.addItem("Anthropic (Claude)", "anthropic")
        self.ai_provider.addItem("OpenAI (e.g. Codex)", "openai")
        self.ai_model = QLineEdit("claude-opus-5-5")
        self.ai_model.setToolTip("Claude: claude-opus-5-5 (default). OpenAI: a model your account offers (or OPENAI_MODEL)")
        self.ai_provider.currentIndexChanged.connect(self._ai_provider_changed)
        self.ai_steps = QSpinBox(minimum=1, maximum=200, value=30)
        self.ai_creds = QLabel()
        self.ai_key = QLineEdit()
        self.ai_key.setEchoMode(QLineEdit.Password)
        self.ai_key.setPlaceholderText("API key (kept in memory; tick Remember to keep it in the macOS Keychain)")
        self.ai_remember = QCheckBox("Remember (macOS Keychain)")
        self.ai_remember.setToolTip("store the key encrypted in the macOS Keychain (service zulf-studio); Studio loads it "
                                    "at start. Forget removes it. Never written to a file or the log.")
        from . import credentials
        self.ai_remember.setEnabled(credentials.available())
        use_key = QPushButton("Use")
        forget_key = QPushButton("Forget")
        check = QPushButton("Check")
        use_key.clicked.connect(self._ai_use_key)
        forget_key.clicked.connect(self._ai_forget_key)
        check.clicked.connect(self._ai_provider_changed)
        self.ai_key_buttons = (use_key, forget_key, check)
        self.ai_guide = QTextBrowser()
        self.ai_guide.setHtml(AI_SETUP_GUIDE)

    def _settings_dialog(self):
        dlg = QDialog(self)
        dlg.setWindowTitle("Settings")
        dlg.resize(820, 640)
        tabs = QTabWidget()
        ai = QWidget()
        form = QFormLayout(ai)
        form.addRow("provider", self.ai_provider)
        form.addRow("model", self.ai_model)
        form.addRow("max steps per request", self.ai_steps)
        key_row = QHBoxLayout()
        key_row.addWidget(self.ai_key, 1)
        for b in self.ai_key_buttons:
            key_row.addWidget(b)
        form.addRow("API key", key_row)
        form.addRow("", self.ai_remember)
        form.addRow("status", self.ai_creds)
        form.addRow(self.ai_guide)
        tabs.addTab(ai, "AI assistant")
        tabs.addTab(self._api_tab(), "AI API")
        look = QWidget()
        lf = QFormLayout(look)
        theme = QComboBox()
        for text, mode in (("follow the system", "auto"), ("light", "light"), ("dark", "dark")):
            theme.addItem(text, mode)
        theme.currentIndexChanged.connect(lambda _: self.apply_theme(theme.currentData()))
        lf.addRow("theme", theme)
        tabs.addTab(look, "Appearance")
        lay = QVBoxLayout(dlg)
        lay.addWidget(tabs)
        close = QPushButton("Close")
        close.clicked.connect(dlg.hide)
        row = QHBoxLayout()
        row.addStretch(1)
        row.addWidget(close)
        lay.addLayout(row)
        self.settings_tabs = tabs
        return dlg

    def open_settings(self, page: str = ""):
        tabs = self.settings_tabs
        names = [tabs.tabText(i) for i in range(tabs.count())]
        if page in names:
            tabs.setCurrentIndex(names.index(page))
        self.settings.show()
        self.settings.raise_()

    def _tools_window(self):
        win = QWidget(self, Qt.Window)
        win.setWindowTitle("ZULF Studio tools")
        win.resize(900, 560)
        tabs = QTabWidget()
        tabs.addTab(Terminal(self.session), "Terminal")
        tabs.addTab(Console({"session": self.session, "api": self.api, "np": np}), "Python")
        QVBoxLayout(win).addWidget(tabs)
        self.tools_tabs = tabs
        return win

    def open_tools(self, page: str = ""):
        names = [self.tools_tabs.tabText(i) for i in range(self.tools_tabs.count())]
        if page in names:
            self.tools_tabs.setCurrentIndex(names.index(page))
        self.tools.show()
        self.tools.raise_()

    def _ai_key_variable(self):
        return "ANTHROPIC_API_KEY" if self.ai_provider.currentData() == "anthropic" else "OPENAI_API_KEY"

    def _ai_has_credentials(self):
        if self.ai_provider.currentData() == "anthropic":
            return bool(os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")
                        or Path.home().joinpath(".config", "anthropic").exists())
        return bool(os.environ.get("OPENAI_API_KEY"))

    def _ai_use_key(self):
        key = self.ai_key.text().strip()
        if not key:
            return
        var = self._ai_key_variable()
        os.environ[var] = key                            # this process; the Keychain only when Remember is ticked
        self.ai_key.clear()
        where = "this session (not saved)"
        if self.ai_remember.isChecked():
            from . import credentials
            where = ("this session and the macOS Keychain" if credentials.save(var, key)
                     else "this session (Keychain save failed)")
        self.session.log(f"{var} set for {where}", "ai")
        self._ai_provider_changed()

    def _ai_forget_key(self):
        var = self._ai_key_variable()
        os.environ.pop(var, None)
        from . import credentials
        removed = credentials.available() and credentials.delete(var)
        self.session.log(f"{var} removed from this session" + (" and the macOS Keychain" if removed else ""), "ai")
        self._ai_provider_changed()

    def _ai_provider_changed(self):
        prov = self.ai_provider.currentData()
        if prov == "anthropic":
            if not self.ai_model.text() or not self.ai_model.text().startswith("claude"):
                self.ai_model.setText("claude-opus-5-5")
        elif self.ai_model.text().startswith("claude"):
            self.ai_model.setText(os.environ.get("OPENAI_MODEL", ""))
        ok = self._ai_has_credentials()
        hint = ("credentials found" if ok else
                f"no credentials: set {self._ai_key_variable()} (see the guide below)")
        self.ai_creds.setText(hint)
        self.ai_creds.setStyleSheet("color: %s" % (self.t["good"] if ok else self.t["bad"]))
        if hasattr(self, "ai_status"):
            name = self.ai_provider.currentText().split(" (")[0]
            model = self.ai_model.text() or "no model set"
            self.ai_status.setText(f"{name} \u00b7 {model} \u00b7 " +
                                   ("ready" if ok else "no API key yet: open Settings to set it up"))
            self.ai_status.setStyleSheet("color: %s" % (self.t["muted"] if ok else self.t["bad"]))
        self.assistant = None                      # a new provider or model starts a new conversation

    def _ai_new(self):
        self.assistant = None
        self.ai_log.appendPlainText("--- new conversation ---")

    def ai_ask(self):
        prompt = self.ai_prompt.toPlainText().strip()
        if not prompt:
            return
        try:
            from .assistant import StudioAssistant
            if self.assistant is None or self.assistant.model != self.ai_model.text().strip() \
                    or self.assistant.provider != self.ai_provider.currentData():
                self.assistant = StudioAssistant(self.session, self.ai_provider.currentData(),
                                                 self.ai_model.text().strip(), self.ai_steps.value(),
                                                 on_message=self.bridge.ai_message.emit)
            self.assistant.max_steps = self.ai_steps.value()
            self.assistant.ask(prompt, background=True)
        except Exception as exc:
            self.on_ai_message("error", f"{type(exc).__name__}: {exc}")
            return
        self.ai_prompt.clear()
        self.ai_send.setEnabled(False)

    def on_ai_message(self, role, text):
        prefix = {"user": "you", "assistant": "AI", "tool": "  tool", "error": "error"}.get(role, role)
        self.ai_log.appendPlainText(f"{prefix}: {text}" + ("\n" if role in ("assistant", "error") else ""))
        if role in ("assistant", "error"):
            self.ai_send.setEnabled(True)

    def _api_tab(self):
        w = QPlainTextEdit(readOnly=True)
        w.setObjectName("mono")
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

    def apply_theme(self, mode=None):
        """Light or dark style sheet and plot colours; "auto" follows the system appearance."""
        if mode is not None:
            self.theme_mode = mode
        dark = self.theme_mode == "dark"
        if self.theme_mode == "auto":
            try:
                dark = QApplication.styleHints().colorScheme() == Qt.ColorScheme.Dark
            except AttributeError:
                dark = False
        self.t = DARK if dark else LIGHT
        QApplication.instance().setStyleSheet(stylesheet(self.t))
        if hasattr(self, "fig"):
            self.fig.set_facecolor(self.t["panel"])
            if hasattr(self, "ai_creds"):
                self._ai_provider_changed()
            self.schedule()

    def _menu(self):
        def action(menu, text, slot, shortcut=None, role=None):
            a = QAction(text, self)
            if shortcut is not None:
                a.setShortcut(QKeySequence(shortcut))
            if role is not None:
                a.setMenuRole(role)
            a.triggered.connect(lambda _=False: slot())
            menu.addAction(a)
            return a
        m = self.menuBar().addMenu("&File")
        action(m, "Open ...", self.open_dialog, QKeySequence.Open)
        self.recent_menu = m.addMenu("Open recent")
        self._fill_recent()
        m.addSeparator()
        action(m, "Import averaged FID (.npy) ...", self.import_fid_dialog, "Ctrl+I")
        action(m, "Import scan folder (n.dat, n.ini) ...", self.import_scans_dialog)
        action(m, "Load series.json ...", self.load_series)
        action(m, "Load fit run ...", self.load_run, "Ctrl+R")
        m.addSeparator()
        action(m, "Save session", self.save_session, QKeySequence.Save)
        action(m, "Save session as ...", lambda: self.save_session(ask=True), QKeySequence.SaveAs)
        m.addSeparator()
        ex = m.addMenu("Export")
        action(ex, "Parameters, couplings, lines and spectrum (JSON, CSV) ...", self.export, "Ctrl+E")
        action(ex, "Plot as image (PNG, PDF, SVG) ...", self.export_plot)
        action(ex, "Publication figure (paper_figure) ...", self.export_figure)
        action(m, "Generate publication figure", self.make_figure, "Ctrl+G")
        m.addSeparator()
        action(m, "Settings ...", lambda: self.open_settings(), QKeySequence.Preferences, QAction.PreferencesRole)
        v = self.menuBar().addMenu("&View")
        for i in range(self.tabs.count()):
            action(v, self.tabs.tabText(i), lambda i=i: self.tabs.setCurrentIndex(i), f"Ctrl+{i + 1}")
        r = self.menuBar().addMenu("&Run")
        action(r, "Start fit", self.start_fit, "Ctrl+Return")
        action(r, "Blind analysis ...", lambda: self.tabs.setCurrentWidget(self.analysis_panel), "Ctrl+B")
        action(r, "Stop running job", self.activity._stop, "Ctrl+.")
        action(r, "Jobs", lambda: self.tabs.setCurrentWidget(self.jobs_panel), "Ctrl+J")
        t = self.menuBar().addMenu("&Tools")
        for text, key, page in (("Terminal", "Ctrl+Shift+T", "Terminal"), ("Python console", "Ctrl+Shift+P", "Python")):
            action(t, text, lambda pg=page: self.open_tools(pg), key)

    # ---- files: open, import, save, export, recent, drag and drop ------------------------------
    def _fill_recent(self):
        self.recent_menu.clear()
        paths = [p for p in (self.settings_store.value("recent", []) or []) if Path(p).exists()]
        for p in paths:
            a = QAction(p.replace(str(Path.home()), "~"), self)
            a.triggered.connect(lambda _=False, pp=p: self.open_any(pp))
            self.recent_menu.addAction(a)
        self.recent_menu.setEnabled(bool(paths))
        if paths:
            self.recent_menu.addSeparator()
            clear = QAction("Clear list", self)
            clear.triggered.connect(lambda: (self.settings_store.setValue("recent", []), self._fill_recent()))
            self.recent_menu.addAction(clear)

    def _remember(self, path):
        paths = [str(path)] + [p for p in (self.settings_store.value("recent", []) or []) if p != str(path)]
        self.settings_store.setValue("recent", paths[:12])
        self._fill_recent()

    def _start_dir(self, key, default):
        return self.settings_store.value(f"dir/{key}", str(default))

    def _remember_dir(self, key, path):
        self.settings_store.setValue(f"dir/{key}", str(Path(path).parent if Path(path).is_file() else path))

    def open_any(self, path):
        """Open a session, series, fit run, FID (imported) or scan folder (averaged, then imported)."""
        if not self._discard_ok():
            return
        r = self._guard(self.session.open_path, path)
        if r:
            self._remember(path)
            self.statusBar().showMessage(f"opened {r['opened']}: {path}", 6000)
            if r["opened"] in ("FID", "scan folder"):
                self.tabs.setCurrentWidget(self.jobs_panel)

    def open_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open", self._start_dir("open", ROOT / "runs"),
            f"ZULF files (*{SESSION_SUFFIX} *.json *.npy);;Studio session (*{SESSION_SUFFIX});;"
            "Series or fit (*.json);;Averaged FID (*.npy);;All files (*)")
        if path:
            self._remember_dir("open", path)
            self.open_any(path)

    def import_fid_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Averaged FID", self._start_dir("fid", Path.home() / "research"),
                                              "Averaged FID (*.npy)")
        if path:
            self._remember_dir("fid", path)
            self.open_any(path)

    def import_scans_dialog(self):
        d = QFileDialog.getExistingDirectory(self, "Instrument run folder (n.dat, n.ini; read-only)",
                                             self._start_dir("scans", Path.home() / "research"))
        if d:
            self._remember_dir("scans", d)
            self.open_any(d)

    def save_session(self, ask=False):
        path = self.session.session_file
        if ask or not path:
            path, _ = QFileDialog.getSaveFileName(self, "Save session", path or self._start_dir(
                "session", self.session.workspace / f"session{SESSION_SUFFIX}"), f"Studio session (*{SESSION_SUFFIX})")
        if path:
            r = self._guard(self.session.save_session, path)
            if r:
                self._remember(r["path"])
                self._remember_dir("session", r["path"])
                self.statusBar().showMessage(f"session saved: {r['path']}", 6000)
                return True
        return False

    def export_plot(self):
        path, _ = QFileDialog.getSaveFileName(self, "Plot as image", self._start_dir("plot", self.session.workspace /
                                                                                  "plot.png"),
                                              "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if path:
            if Path(path).suffix.lower() not in (".png", ".pdf", ".svg"):
                path += ".png"
            self.fig.savefig(path, dpi=200)
            self._remember_dir("plot", path)
            self.statusBar().showMessage(f"plot saved: {path}", 6000)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls() and all(u.isLocalFile() for u in e.mimeData().urls()):
            e.acceptProposedAction()

    def dropEvent(self, e):
        for u in e.mimeData().urls():
            self.open_any(u.toLocalFile())

    # ---- title, unsaved changes, layout ---------------------------------------------------------
    def update_title(self):
        name = Path(self.session.session_file).stem if self.session.session_file else "untitled"
        self.setWindowTitle(f"{name}{' \u2022' if self.dirty else ''} \u2014 ZULF Studio")
        self.setWindowModified(self.dirty)

    def _discard_ok(self):
        """Ask before replacing unsaved couplings (not in offscreen tests)."""
        if not self.dirty or os.environ.get("QT_QPA_PLATFORM") == "offscreen":
            return True
        r = QMessageBox.question(self, "Unsaved changes", "Save the current session first?",
                                 QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        if r == QMessageBox.Save:
            return self.save_session()
        return r == QMessageBox.Discard

    def _restore_layout(self):
        geo = self.settings_store.value("geometry")
        if geo is not None:
            self.restoreGeometry(geo)
        for name, sp in self.splitters.items():
            state = self.settings_store.value(f"splitter/{name}")
            if state is not None:
                sp.restoreState(state)
        screen = QApplication.primaryScreen()
        if screen is not None and not screen.availableGeometry().intersects(self.frameGeometry()):
            self.resize(min(1500, screen.availableGeometry().width()), min(950, screen.availableGeometry().height()))
            self.move(screen.availableGeometry().topLeft())

    def closeEvent(self, e):
        running = [st for st in self.session.job_list() if st["running"]]
        if running and os.environ.get("QT_QPA_PLATFORM") != "offscreen":
            r = QMessageBox.question(self, "Jobs running", f"{len(running)} job(s) still running ("
                                     + ", ".join(st["title"] for st in running) + "). They keep running after "
                                     "Studio closes. Close anyway?", QMessageBox.Yes | QMessageBox.No)
            if r != QMessageBox.Yes:
                e.ignore()
                return
        if not self._discard_ok():
            e.ignore()
            return
        self.settings_store.setValue("geometry", self.saveGeometry())
        for name, sp in self.splitters.items():
            self.settings_store.setValue(f"splitter/{name}", sp.saveState())
        super().closeEvent(e)

    def update_jobs(self):
        self.activity.refresh()
        self.machine.refresh()
        self.jobs_panel.refresh()
        self.analysis_panel.refresh()
        running = self.activity.running
        if running:
            self.busy_pill.setText(SPIN_TEXT[self.activity.tick % len(SPIN_TEXT)] + "  " + job_line(running[0]))
            self.busy_pill.show()
        else:
            self.busy_pill.hide()

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
        if not confirm_workers(self, self.session, self.f_workers.value()):
            return
        extra = self.f_extra.text().split()
        self._guard(self.session.start_fit, starts=self.f_starts.value(), workers=self.f_workers.value(),
                    max_nfev=self.f_nfev.value(), trace=self.f_trace.value(), fit_field=self.f_field.isChecked(),
                    precision=self.f_precision.currentData(),
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
        if event in ("structure", "couplings", "field", "linewidth"):
            self.dirty = True
        elif event in ("session_saved", "session_opened"):
            self.dirty = False
        if event in ("fit_started", "blind_started", "import_started"):
            self.update_jobs()
        if hasattr(self, "busy_pill"):
            self.update_title()
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
        name = s.spec.get("compound") or s.spec.get("motif") or ""
        comps = ", ".join(c["label"] for c in s.components())
        self.subtitle.setText(f"{name}   \u00b7   {comps}" + (f"   \u00b7   data: {s.data['label']}" if s.data else ""))
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
        self.phase.set_value(s.data_phase_deg)
        self.delay.set_value(s.data_delay_ms)
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
        t = self.t
        style = matplotlib.rc_context(matplotlib_style(t))
        style.__enter__()
        self.fig.clear()
        self.fig.set_facecolor(t["panel"])
        has_data = "data_re" in sim
        rows = [3, 1, 1] if has_data else [3, 1]
        axes = self.fig.subplots(len(rows), 1, sharex=True, gridspec_kw={"height_ratios": rows})
        ax = axes[0]
        if has_data:
            d = pick(sim["data_re"], sim["data_im"])
            ax.plot(f, d, color=t["data"], lw=0.8, label=f"experiment ({s.data['label']})")
        m = pick(sim["sim_re"], sim["sim_im"])
        ax.plot(f, m, color=t["sim"], lw=1.3, alpha=0.9, label="simulation (quick look)")
        tr = s.trace if self.show_trace.isChecked() else None
        if tr is not None and s.trace_index >= 0:
            ft = tr["f"]
            sel = (ft >= s.view[0]) & (ft <= s.view[1])
            ph = np.exp(1j * (np.radians(s.data_phase_deg) + 2 * np.pi * ft * 1e-3 * s.data_delay_ms))
            mod = tr["models"][s.trace_index] * ph
            mt = mod.real if part == "re" else mod.imag if part == "im" else np.abs(mod)
            ax.plot(ft[sel], mt[sel], color=t["trace"], lw=1.2,
                    label=f"fit trace frame {s.trace_index + 1} (objective "
                          f"{tr['meta']['frames'][s.trace_index]['objective']:.4g})")
        ax.legend(loc="upper right")
        ax.set_ylabel("signal")
        if has_data:
            axes[1].plot(f, d - m, color=t["resid"], lw=0.7)
            axes[1].axhline(0, color=t["line2"], lw=0.6)
            axes[1].set_ylabel("residual")
        sax = axes[-1]
        if self.show_sticks.isChecked():
            labels = [c["label"] for c in s.components()]
            segs, cols = [], []
            for r in sim["lines"]:
                segs.append([(r["frequency_hz"], 0), (r["frequency_hz"], r["relative"])])
                cols.append(COMPONENT_COLORS[labels.index(r["component"]) % 6] if r["component"] in labels else t["muted"])
            sax.add_collection(LineCollection(segs, colors=cols, linewidths=1.6))
            sax.set_ylim(0, 1.08)
        sax.grid(False)
        sax.set_yticks([])
        sax.spines["left"].set_visible(False)
        sax.set_ylabel("lines")
        sax.set_xlim(*s.view)
        sax.set_xlabel("frequency (Hz)")
        style.__exit__(None, None, None)
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
    from . import credentials
    loaded = credentials.load_into_environment()
    if loaded:
        session.log("API keys loaded from the macOS Keychain: " + ", ".join(loaded), "ai")
    if args.series:
        session.load_spectrum(series=args.series)
    if args.fit:
        session.apply_fit(args.fit)
    server = None
    if args.api_port >= 0:
        try:
            server = serve(session, port=args.api_port, background=not args.no_gui)
        except OSError as exc:                       # port taken (another Studio): any free port instead
            session.log(f"API port {args.api_port} unavailable ({exc}); using a free port", "api")
            server = serve(session, port=0, background=not args.no_gui)
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
