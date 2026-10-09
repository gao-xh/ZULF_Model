"""Live fit monitor inside Studio (the native view of scripts/fit_monitor.py).

A fit records itself in OUT/monitor (fit_joint_series --monitor, on by default): status.json, one jsonl per start
and console.log. This page reads those files every two seconds while it is visible: the objective of every start,
the data and the model at the best point so far of the selected start (rendered here, in a background thread,
from the run's recorded command; the fit itself is not touched), the starts table, the couplings of the selected
start and the console. A run whose records do not exist yet (a fit that is still setting up) shows a waiting line
instead of an error.
"""
from __future__ import annotations

import sys
import threading
import time
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QHBoxLayout, QHeaderView, QLabel, QPlainTextEdit,
                               QPushButton, QSlider, QSplitter, QTabWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from .flow import FlowLayout, flow_policy
from .session import ROOT, _resolve
from .theme import ISOTOPOLOGUE, matplotlib_style

sys.path.insert(0, str(ROOT / "scripts"))
from fit_monitor import SpectrumCache, find_runs, read_point, read_run  # noqa: E402

RUNS_ROOT = ROOT / "runs"
START_COLORS = list(ISOTOPOLOGUE) + ["#8c564b", "#17becf", "#bcbd22", "#7f7f7f", "#e377c2", "#1f3b73"]


def _short(path) -> str:
    p = str(path)
    return p[len(str(ROOT)) + 1:] if p.startswith(str(ROOT) + "/") else p


def _fmt(v):
    return "-" if v is None or not np.isfinite(v) else f"{v:.5g}"


class _Relay(QObject):
    spectrum = Signal(object)


class MonitorPanel(QWidget):
    """Objective chart, spectrum at the best point, starts, couplings and console of one fit run."""

    STARTS = ("#", "state", "stage", "evals", "now", "best", "final", "s")

    def __init__(self, session, window, parent=None, large=False):
        super().__init__(parent)
        self.session, self.window, self.large = session, window, large
        self.popout = None
        self.cache = SpectrumCache()
        self.relay = _Relay()
        self.relay.spectrum.connect(self._got_spectrum)
        self.run = None
        self.data = None
        self.selected = None
        self.spec = None
        self.spec_busy = False
        self.spec_key = None
        self.point = None                                  # evaluation shown (None: the latest best point)
        self.ax_obj = None
        lay = QVBoxLayout(self)
        top = QHBoxLayout()
        self.runs = QComboBox()
        self.runs.setMinimumContentsLength(12)
        self.runs.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.runs.setToolTip("fit runs with a monitor record below runs/, newest first")
        self.runs.activated.connect(lambda _i: self.set_run(self.runs.currentData()))
        self.b_reload = QPushButton("Reload list", objectName="small")
        self.b_reload.clicked.connect(self.fill_runs)
        self.follow = QCheckBox("follow", checked=True)
        self.follow.setToolTip("re-render the spectrum whenever the best point of the selected start improves")
        self.part = QComboBox()
        for p in ("real", "imaginary", "magnitude"):
            self.part.addItem(p)
        self.part.currentIndexChanged.connect(lambda _i: self._draw())
        top.addWidget(QLabel("run"))
        top.addWidget(self.runs, 1)
        top.addWidget(self.b_reload)
        tools_box = flow_policy(QWidget())
        tools = FlowLayout(tools_box, spacing=12)
        top.addWidget(tools_box)
        tools.addWidget(self.follow)
        tools.addWidget(self.part)
        if not large:
            self.b_large = QPushButton("Large window", objectName="small")
            self.b_large.setToolTip("open this monitor in its own large window (resize or maximize it freely)")
            self.b_large.clicked.connect(self.open_large)
            tools.addWidget(self.b_large)
        lay.addLayout(top)
        self.status = QLabel("", wordWrap=True)
        lay.addWidget(self.status)
        prow = QHBoxLayout()                               # which point of which start is drawn
        self.point_slider = QSlider(Qt.Horizontal)
        self.point_slider.setMinimum(1)
        self.point_slider.setToolTip("evaluation of the selected start whose spectrum is drawn; or click any curve "
                                     "of the objective plot (that start, that evaluation)")
        self.point_timer = QTimer(self, singleShot=True, interval=150)
        self.point_timer.timeout.connect(lambda: self.set_point(self.point_slider.value()))
        self.point_slider.valueChanged.connect(lambda _v: None if self._quiet else self.point_timer.start())
        self.point_label = QLabel("latest best point")
        self.point_label.setMinimumWidth(120)
        self.b_latest = QPushButton("Latest best", objectName="small")
        self.b_latest.setToolTip("back to the best point so far of the selected start (follows the run)")
        self.b_latest.clicked.connect(lambda: self.set_point(None))
        prow.addWidget(QLabel("point"))
        prow.addWidget(self.point_slider, 1)
        prow.addWidget(self.point_label)
        prow.addWidget(self.b_latest)
        lay.addLayout(prow)
        self._quiet = False
        split = QSplitter(Qt.Horizontal)                  # plots | starts, couplings, console
        self.fig = Figure(figsize=(6, 4.2), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumSize(200, 200)
        split.addWidget(self.canvas)
        self.canvas.mpl_connect("button_press_event", self._clicked)
        tables = QTabWidget()
        self.starts = QTableWidget(0, len(self.STARTS))
        self.starts.setHorizontalHeaderLabels(self.STARTS)
        self.couplings = QTableWidget(0, 3)
        self.couplings.setHorizontalHeaderLabels(("coupling", "J (Hz)", "from start"))
        for t in (self.starts, self.couplings):
            t.verticalHeader().setVisible(False)
            t.verticalHeader().setDefaultSectionSize(22)
            t.setShowGrid(False)
            t.setAlternatingRowColors(True)
            t.setSelectionBehavior(QTableWidget.SelectRows)
            t.setEditTriggers(QTableWidget.NoEditTriggers)
            t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
            t.horizontalHeader().setStretchLastSection(True)
            t.setMinimumHeight(60)
        self.starts.itemSelectionChanged.connect(self._pick_start)
        tables.addTab(self.couplings, "Couplings")
        tables.addTab(self.starts, "Starts")
        self.console = QPlainTextEdit(readOnly=True)
        self.console.setObjectName("mono")
        self.console.setMaximumBlockCount(400)
        self.console.setMinimumHeight(40)
        tables.addTab(self.console, "Console")
        tables.setMinimumWidth(200)
        split.addWidget(tables)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)
        lay.addWidget(split, 1)
        self.timer = QTimer(self, interval=2000)
        self.timer.timeout.connect(self.refresh)
        self.timer.start()
        self.fill_runs()

    def open_large(self):
        """The same monitor in a separate window of 85 % of the screen, on the run shown here."""
        if self.popout is None:
            self.popout = MonitorWindow(self.session, self.window)
        if self.run:
            self.popout.panel.set_run(self.run)
        self.popout.show()
        self.popout.raise_()
        self.popout.activateWindow()

    # ---- run choice ------------------------------------------------------------------------
    def fill_runs(self):
        current = self.run
        self.runs.blockSignals(True)
        self.runs.clear()
        paths = [str(p.resolve()) for p in find_runs(RUNS_ROOT)]
        if current and current not in paths:
            paths.insert(0, current)
        for p in paths:
            self.runs.addItem(_short(p), p)
        if current:
            self.runs.setCurrentIndex(paths.index(current))
        self.runs.blockSignals(False)
        if not current and paths:
            self.set_run(paths[0])

    def set_run(self, run):
        """Show this run directory (it may not have a monitor record yet)."""
        run = str(_resolve(run).resolve()) if run else None
        if run != self.run:
            self.run, self.data, self.selected, self.spec, self.spec_key = run, None, None, None, None
            self.point = None
            self.starts.setRowCount(0)
            self.couplings.setRowCount(0)
            self.console.clear()
        if run and self.runs.findData(run) < 0:
            self.runs.insertItem(0, _short(run), run)
        if run:
            self.runs.setCurrentIndex(self.runs.findData(run))
        self.refresh(force=True)

    # ---- polling ---------------------------------------------------------------------------
    def refresh(self, force=False):
        if not self.run or (not force and not self.isVisible()):
            return
        if not (Path(self.run) / "monitor" / "status.json").exists():
            self.status.setText(f"waiting for the monitor record of {_short(self.run)} (the fit writes it once its "
                                "set-up and component search are done; this page updates by itself)")
            self._draw()
            return
        if self.data and self.data["status"].get("phase") == "finished" and not force:
            return
        try:
            self.data = read_run(self.run)
        except Exception as exc:                          # a file being rewritten: try again next tick
            self.status.setText(f"reading the record: {type(exc).__name__}: {exc}")
            return
        st, starts = self.data["status"], self.data["starts"]
        end = st.get("finished_at") or self.data["now"]
        best = min((v["best"] for v in starts.values() if v["best"] is not None), default=None)
        alive = {True: "process alive", False: "process gone", None: ""}[self.data["alive"]]
        self.status.setText(f"phase <b>{st['phase']}</b> &middot; {(end - st['started']) / 60:.1f} min &middot; "
                            f"starts finished {len(st['finished'])}/{st['starts']} &middot; best {_fmt(best)}"
                            + (f" &middot; {alive}" if alive and st["phase"] != "finished" else ""))
        self._fill_starts()
        if self.selected not in starts and starts:
            self.selected = min(starts, key=lambda k: np.inf if starts[k]["best"] is None else starts[k]["best"])
        self._fill_couplings()
        self._sync_slider()
        log = Path(self.run) / "monitor" / "console.log"
        text = log.read_text(errors="replace")[-6000:] if log.exists() else ""
        if text != self.console.toPlainText():
            self.console.setPlainText(text)
            self.console.verticalScrollBar().setValue(self.console.verticalScrollBar().maximum())
        self._request_spectrum()
        self._draw()

    def _fill_starts(self):
        st, starts = self.data["status"], self.data["starts"]
        names = list(starts)
        self.starts.blockSignals(True)
        self.starts.setRowCount(len(names))
        for i, name in enumerate(names):
            s = starts[name]
            fin = st["finished"].get(name.replace("start_", "").lstrip("0") or "0")
            vals = (name.replace("start_", ""), "running" if s["running"] else "done", s["label"],
                    str(s["evaluations"]), _fmt(s["cost"]), _fmt(s["best"]), _fmt(fin), f"{s['seconds']:.0f}")
            for j, v in enumerate(vals):
                item = self.starts.item(i, j)
                if item is None:
                    item = QTableWidgetItem()
                    self.starts.setItem(i, j, item)
                if item.text() != v:
                    item.setText(v)
            if name == self.selected:
                self.starts.selectRow(i)
        self.starts.blockSignals(False)

    def _pick_start(self):
        r = self.starts.currentRow()
        names = list(self.data["starts"]) if self.data else []
        if 0 <= r < len(names) and names[r] != self.selected:
            self.selected, self.spec, self.spec_key, self.point = names[r], None, None, None
            self._fill_couplings()
            self._request_spectrum()
            self._draw()

    def _fill_couplings(self):
        s = self.data["starts"].get(self.selected) if self.data else None
        J = (s or {}).get("J") or {}
        if self.point is not None:                         # the couplings of the drawn point
            J = (self.spec or {}).get("J") or {}
        J0 = (s or {}).get("J0") or self.data["status"].get("start_couplings") or {} if self.data else {}
        rows = []
        for key, v in J.items():
            v0 = J0.get(key)
            rows.append((key, ", ".join(f"{x:.3f}" for x in v),
                         ", ".join(f"{a - b:+.3f}" for a, b in zip(v, v0)) if v0 else ""))
        self.couplings.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, v in enumerate(row):
                item = self.couplings.item(i, j)
                if item is None:
                    item = QTableWidgetItem()
                    self.couplings.setItem(i, j, item)
                if item.text() != v:
                    item.setText(v)

    # ---- the point shown: click on the objective plot, slider, Latest best ---------------------
    def set_point(self, n, start=None):
        """Draw the spectrum at evaluation n of the start (default: the selected one); None: the latest best."""
        if start is not None and start != self.selected:
            self.selected = start
            self._fill_starts()
        self.point = None if n is None else int(n)
        self.spec, self.spec_key = None, None
        self._sync_slider()
        self._request_spectrum()
        self._fill_couplings()
        self._draw()

    def _clicked(self, event):
        """Click near a curve of the objective plot: that start, the nearest evaluation."""
        if self.ax_obj is None or event.inaxes is not self.ax_obj or not self.data or event.xdata is None:
            return
        best = None
        for name, s in self.data["starts"].items():
            h = np.asarray(s["history"], float)
            if not len(h):
                continue
            i = int(np.argmin(np.abs(h[:, 0] - event.xdata)))
            for col in (1, 2):                         # thin (this evaluation) and thick (best so far)
                px = self.ax_obj.transData.transform((h[i, 0], h[i, col]))
                d = np.hypot(px[0] - event.x, px[1] - event.y)
                if best is None or d < best[0]:
                    best = (d, name, int(h[i, 0]))
        if best is not None and best[0] < 40:          # pixels
            self.set_point(best[2], best[1])

    def _sync_slider(self):
        if self.point_slider.isSliderDown():
            return
        s = self.data["starts"].get(self.selected) if self.data else None
        self._quiet = True
        self.point_slider.setMaximum(max(1, (s or {}).get("evaluations", 1)))
        self.point_slider.setValue(self.point or self.point_slider.maximum())
        self._quiet = False
        if self.point is None:
            self.point_label.setText("latest best point")
        else:
            self.point_label.setText(f"evaluation {self.point} / {(s or {}).get('evaluations', '?')}")

    # ---- spectrum at the best point (background thread) ------------------------------------
    def _request_spectrum(self):
        if not self.data or self.selected not in self.data["starts"] or self.spec_busy:
            return
        s = self.data["starts"][self.selected]
        point = self.point
        if point is None and s["z"] is None:
            return
        key = (self.run, self.selected, ("best", s["best"]) if point is None else ("point", point))
        if key == self.spec_key or (point is None and self.spec is not None and not self.follow.isChecked()):
            return
        self.spec_busy, self.spec_key = True, key
        run, start = self.run, self.selected
        keys = self.data["status"].get("keys", [])

        def work():
            try:
                rp = None if point is None else read_point(run, start, point)
                z = None if rp is None else rp["z"]
                out = self.cache.spectrum(run, start, z=z)
                if rp is not None:
                    joint = self.cache.get(run).joint
                    zz = np.asarray(rp["z"], float)
                    out["J"] = {k: joint.coupling_values(zz, i).tolist() for i, k in enumerate(keys)}
                out["point"] = rp
                out["key"] = key
            except BaseException as exc:                 # argparse exits with SystemExit: report, never hang
                out = {"key": key, "error": f"{type(exc).__name__}: {exc}"}
            self.relay.spectrum.emit(out)
        threading.Thread(target=work, daemon=True, name="monitor-spectrum").start()

    def _got_spectrum(self, out):
        self.spec_busy = False
        if out["key"] != self.spec_key:                    # the run, start or point changed meanwhile: render that
            self.spec_key = None
            self._request_spectrum()
            self._draw()
            return
        self.spec = out
        self._fill_couplings()
        self._draw()

    # ---- drawing ---------------------------------------------------------------------------
    def _draw(self):
        t = self.window.t
        with matplotlib.rc_context(matplotlib_style(t)):
            self.fig.clear()
            self.fig.set_facecolor(t["panel"])
            ax1, ax2 = self.fig.subplots(2, 1, height_ratios=[1, 1.6])
            self.ax_obj = ax1
            starts = self.data["starts"] if self.data else {}
            for i, (name, s) in enumerate(starts.items()):
                h = np.asarray(s["history"], float)
                if not len(h):
                    continue
                col = START_COLORS[i % len(START_COLORS)]
                top = name == self.selected
                ax1.plot(h[:, 0], h[:, 1], color=col, lw=0.6, alpha=0.35, zorder=3 if top else 1)
                ax1.plot(h[:, 0], h[:, 2], color=col, lw=2.0 if top else 1.2, zorder=4 if top else 2,
                         label=f"#{name.replace('start_', '')}" if len(starts) <= 8 else None)
            sel = starts.get(self.selected)
            if sel is not None and self.point is not None and len(sel["history"]):
                h = np.asarray(sel["history"], float)
                i = int(np.argmin(np.abs(h[:, 0] - self.point)))
                ax1.axvline(self.point, color=t["accent"], lw=0.9, alpha=0.8, zorder=5)
                ax1.plot([h[i, 0]], [h[i, 1]], "o", ms=6, mfc="none", mec=t["accent"], mew=1.6, zorder=6)
            if starts:
                ax1.set_yscale("log")
                ax1.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
                ax1.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:.3g}"))
                ax1.yaxis.set_minor_formatter(NullFormatter())
                if len(starts) <= 8:
                    ax1.legend(loc="upper right", ncols=min(len(starts), 4))
            ax1.set_xlabel("evaluation")
            ax1.set_ylabel("objective")
            ax1.set_title("objective per start (thick: best so far) - click a curve to draw that point", loc="left")
            if self.spec and "error" not in self.spec:
                sp = self.spec
                f = np.asarray(sp["f"])
                d = np.asarray(sp["data_re"]) + 1j * np.asarray(sp["data_im"])
                m = np.asarray(sp["model_re"]) + 1j * np.asarray(sp["model_im"])
                part = self.part.currentText()
                pick = {"real": np.real, "imaginary": np.imag, "magnitude": np.abs}[part]
                for lo, hi in sp["ranges"]:
                    if lo > sp["extent"][0] or hi < sp["extent"][1]:
                        ax2.axvspan(lo, hi, color=t["grid"], alpha=0.35, lw=0)
                ax2.plot(f, pick(d), color=t["ink"], lw=0.8, label="data")
                names = list(starts)
                col = START_COLORS[names.index(self.selected) % len(START_COLORS)] if self.selected in names \
                    else ISOTOPOLOGUE[1]                    # the colour of the start's curve above
                ax2.plot(f, pick(m), color=col, lw=1.2, label="model")
                ax2.plot(f, pick(d) - pick(m), color=t["muted"], lw=0.6, alpha=0.8, label="residual")
                ax2.legend(loc="upper right", ncols=3)
                pt = sp.get("point")
                where = ("best point so far" if pt is None else
                         f"evaluation {pt['n']}" if pt["kind"] == "evaluation" else
                         f"best point up to evaluation {self.point} (older record: evaluation {pt['n']})")
                ax2.set_title(f"start {self.selected.replace('start_', '')}, {where}: "
                              f"objective {sp['objective']:.5g}", loc="left")
            else:
                msg = (self.spec or {}).get("error") or ("rendering the spectrum ..." if self.spec_busy else
                                                         "no best point recorded yet")
                ax2.text(0.5, 0.5, msg, transform=ax2.transAxes, ha="center", va="center", color=t["muted"])
                ax2.set_title("spectrum at the best point so far", loc="left")
            ax2.set_xlabel("frequency (Hz)")
            ax2.set_yticks([])
        self.canvas.draw_idle()

    def showEvent(self, e):
        super().showEvent(e)
        self.refresh(force=True)


class MonitorWindow(QWidget):
    """A separate, large window with its own MonitorPanel (Monitor tab > Large window)."""

    def __init__(self, session, window):
        super().__init__(None, Qt.Window)
        self.setWindowTitle("ZULF Studio - fit monitor")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(10, 10, 10, 10)
        self.panel = MonitorPanel(session, window, self, large=True)
        lay.addWidget(self.panel)
        screen = (window.screen() if window is not None else None) or QApplication.primaryScreen()
        if screen is not None:
            g = screen.availableGeometry()
            self.resize(int(g.width() * 0.85), int(g.height() * 0.85))
            self.move(g.x() + (g.width() - self.width()) // 2, g.y() + (g.height() - self.height()) // 2)
