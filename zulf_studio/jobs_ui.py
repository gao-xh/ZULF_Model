"""Background-job widgets of ZULF Studio: the activity indicator of the status bar and header, the machine-load
label, the Jobs tab and the Analysis tab (blind analysis of an averaged FID).

Every job is a subprocess (session.FitJob: a fit, a blind analysis, an import, a figure). The widgets poll the
session once a second; a running job shows a moving bar (filled by finished starts when the job has a start
count, otherwise indeterminate), its stage, the elapsed time and a Stop button, on every tab.
"""
from __future__ import annotations

import os
import time
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QFormLayout, QGridLayout, QHBoxLayout,
                               QHeaderView, QLabel,
                               QLineEdit, QMessageBox, QProgressBar, QPushButton, QSpinBox, QTableWidget,
                               QTableWidgetItem, QTextBrowser, QVBoxLayout, QWidget)

SPINNER = "\u25d0\u25d3\u25d1\u25d2"
from .flow import FlowLayout, flow_policy

KIND_NAMES = {"fit": "Fit", "blind": "Blind analysis", "figure": "Figure", "import": "Import"}


def elapsed(seconds: float) -> str:
    s = int(seconds)
    return f"{s // 3600}:{s // 60 % 60:02d}:{s % 60:02d}" if s >= 3600 else f"{s // 60}:{s % 60:02d}"


def job_line(st: dict) -> str:
    """One line for a job: kind, title, progress, stage, elapsed, best objective."""
    parts = [KIND_NAMES.get(st["kind"], st["kind"])]
    if st.get("starts_total"):
        parts.append(f"{st['starts_finished']}/{st['starts_total']} starts")
    if st.get("phase") and st["running"]:
        parts.append(st["phase"])
    parts.append(elapsed(st["seconds"]))
    if st.get("best_objective") is not None:
        parts.append(f"best {st['best_objective']:.4g}")
    return "  \u00b7  ".join(parts)


def open_folder(path):
    if path:
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))


class ActivityIndicator(QWidget):
    """Status-bar widget: spinner, running job, thin progress bar, Stop. Hidden text 'idle' when nothing runs."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.tick = 0
        lay = QHBoxLayout(self)
        lay.setContentsMargins(6, 0, 6, 0)
        lay.setSpacing(8)
        self.icon = QLabel("\u25cb", objectName="idle")
        self.text = QLabel("idle", objectName="idle")
        self.bar = QProgressBar(textVisible=False)
        self.bar.setFixedWidth(140)
        self.stop = QPushButton("Stop", objectName="small")
        self.stop.clicked.connect(self._stop)
        for w in (self.icon, self.text, self.bar, self.stop):
            lay.addWidget(w)
        self.running = []
        self.refresh()

    def _stop(self):
        own = [st for st in self.running if not st.get("external")]
        if own:
            self.session.stop_job(own[0]["index"])

    def refresh(self):
        self.running = [st for st in self.session.job_list() if st["running"]]
        self.tick += 1
        if not self.running:
            done = next(iter(self.session.job_list()), None)
            self.icon.setText("\u25cb")
            self.icon.setObjectName("idle")
            self.text.setObjectName("idle")
            self.text.setText("idle" if done is None else
                              f"idle  \u00b7  last: {KIND_NAMES.get(done['kind'], done['kind'])} "
                              f"{'finished' if done['returncode'] == 0 else 'stopped (code %s)' % done['returncode']}"
                              f" in {elapsed(done['seconds'])}")
            self.bar.hide()
            self.stop.hide()
        else:
            st = self.running[0]
            self.icon.setText(SPINNER[self.tick % len(SPINNER)])
            self.icon.setObjectName("busy")
            self.text.setObjectName("busy")
            more = f"  (+{len(self.running) - 1} more)" if len(self.running) > 1 else ""
            self.text.setText(job_line(st) + more)
            if st.get("starts_total"):
                self.bar.setRange(0, int(st["starts_total"]))
                self.bar.setValue(int(st["starts_finished"]))
            else:
                self.bar.setRange(0, 0)               # indeterminate: Qt animates it
            self.bar.show()
            self.stop.setVisible(not st.get("external"))
        for w in (self.icon, self.text):
            w.style().unpolish(w)
            w.style().polish(w)


class MachineLabel(QLabel):
    """Cores, load and the analysis workers running on this machine (all programs); red when overcommitted."""

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        self.status = None
        self.last = 0.0

    def refresh(self, force=False):
        if not force and time.time() - self.last < 5.0:          # ps every 5 s at most
            return self.status
        self.last = time.time()
        m = self.status = self.session.machine_status()
        c = m["cores"]
        cores = f"{c['logical']} cores" + (f" ({c['performance']}P+{c['efficiency']}E)" if c["performance"] else "")
        load = f"load {m['load_average'][0]:.1f}" if m["load_average"] else ""
        over = m["busy_workers"] > c["logical"]
        self.setText(f"CPU {m['busy_workers']}/{c['logical']}" + (f"  \u00b7  {load}" if load else ""))
        self.setMinimumWidth(self.sizeHint().width())
        self.setObjectName("warn" if over else "idle")
        self.setToolTip(f"{cores}; analysis workers in use / cores (all programs on this machine)\n" +
                        ("\n".join(f"pid {p['pid']}: {p['script']}, {p['workers']} workers"
                                    f"{' (this Studio)' if p['own'] else ''}  {p['out']}" for p in m["processes"])
                         or "no analysis processes running"))
        self.style().unpolish(self)
        self.style().polish(self)
        return m


def confirm_workers(parent, session, workers: int) -> bool:
    """Warn before a job would put more analysis workers on the machine than it has cores."""
    m = session.machine_status()
    cores = m["cores"]["logical"]
    if m["busy_workers"] + workers <= cores:
        return True
    box = QMessageBox(parent)
    box.setIcon(QMessageBox.Warning)
    box.setWindowTitle("The machine is busy")
    box.setText(f"{m['busy_workers']} analysis workers already run on {cores} cores; {workers} more would "
                f"slow every run down.")
    box.setInformativeText(f"Suggested: {m['suggested_workers']} worker(s), or wait for a run to finish.\n\n" +
                           "\n".join(f"{p['script']} ({p['workers']} workers) {p['out']}" for p in m["processes"]))
    box.setStandardButtons(QMessageBox.Cancel)
    start = box.addButton("Start anyway", QMessageBox.AcceptRole)
    box.exec()
    return box.clickedButton() is start


class JobsPanel(QWidget):
    """Every job of the session with its state; Stop, Open folder, Apply (fits), Report (blind)."""

    COLUMNS = ("Job", "State", "Progress", "Stage", "Elapsed", "Best", "Output")

    def __init__(self, session, parent=None):
        super().__init__(parent)
        self.session = session
        lay = QVBoxLayout(self)
        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(22)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(len(self.COLUMNS) - 1, QHeaderView.Stretch)
        lay.addWidget(self.table, 1)
        row_box = flow_policy(QWidget())
        row = FlowLayout(row_box, spacing=12)
        self.b_stop = QPushButton("Stop")
        self.b_open = QPushButton("Open folder")
        self.b_apply = QPushButton("Apply fit result")
        self.b_report = QPushButton("Open report")
        self.b_monitor = QPushButton("Live monitor")
        self.b_monitor.setToolTip("live fit monitor of this job in the Monitor tab (objective per start, current "
                                  "simulated spectrum, couplings, console); double-click a job does the same")
        for b in (self.b_stop, self.b_monitor, self.b_open, self.b_apply, self.b_report):
            row.addWidget(b)
        lay.addWidget(row_box)
        self.table.setToolTip("Jobs run as separate processes with one BLAS thread per worker; fits started outside "
                              "Studio are listed too. The status bar shows the running one on every tab.")
        self.b_stop.clicked.connect(lambda: self._with(lambda st: self.session.stop_job(st["index"])))
        self.b_open.clicked.connect(lambda: self._with(lambda st: open_folder(st["out_dir"])))
        self.b_apply.clicked.connect(lambda: self._with(lambda st: self.session.apply_fit(st["out_dir"])))
        self.b_report.clicked.connect(lambda: self._with(self._report))
        self.b_monitor.clicked.connect(lambda: self._with(self._monitor))
        self.table.cellDoubleClicked.connect(lambda *_: self._with(self._monitor))
        self.rows = []
        self.on_monitor = None

    def _monitor(self, st):
        if self.on_monitor is not None:                   # Studio: the Monitor tab of the drawer
            self.on_monitor(st["out_dir"])
            return
        r = self.session.monitor_url(st["out_dir"])
        QDesktopServices.openUrl(QUrl(r["url"]))
        if not r["has_monitor"]:
            self.session.log(f"no monitor record in {st['out_dir']} yet; the page lists every run below runs/",
                             "studio")

    def _report(self, st):
        out = Path(st["out_dir"])
        report = next((p for p in (out / "blind.md", out / "structure.md", out / "RUN_LOG.md") if p.exists()), None)
        open_folder(report or out)

    def _with(self, fn):
        r = self.table.currentRow()
        if 0 <= r < len(self.rows):
            try:
                fn(self.rows[r])
            except Exception as exc:
                QMessageBox.warning(self, "Job", f"{type(exc).__name__}: {exc}")

    def refresh(self):
        self.rows = self.session.job_list()
        self.table.setRowCount(len(self.rows))
        for i, st in enumerate(self.rows):
            state = "running" if st["running"] else ("finished" if st["returncode"] == 0 else
                                                     f"stopped / failed ({st['returncode']})")
            prog = f"{st['starts_finished']}/{st['starts_total']} starts" if st.get("starts_total") else ""
            best = "" if st.get("best_objective") is None else f"{st['best_objective']:.5g}"
            vals = (f"{KIND_NAMES.get(st['kind'], st['kind'])}: {st['title']}", state, prog,
                    st.get("phase", "") if st["running"] else "", elapsed(st["seconds"]), best, st["out_dir"])
            for j, v in enumerate(vals):
                item = self.table.item(i, j)
                if item is None:
                    item = QTableWidgetItem()
                    self.table.setItem(i, j, item)
                if item.text() != v:
                    item.setText(v)
        sel = self.rows[self.table.currentRow()] if 0 <= self.table.currentRow() < len(self.rows) else None
        self.b_stop.setEnabled(bool(sel and sel["running"] and not sel.get("external")))
        self.b_apply.setEnabled(bool(sel and sel["kind"] == "fit" and sel["has_result"]))
        self.b_report.setEnabled(bool(sel and sel["kind"] in ("blind", "fit")))
        self.b_open.setEnabled(sel is not None)
        self.b_monitor.setEnabled(sel is not None)


class AnalysisPanel(QWidget):
    """Blind analysis (scripts/analyze_sample.py) of an averaged FID: no structure, or the current structure."""

    def __init__(self, session, window, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.fid = QLineEdit(placeholderText="default: the FID of the loaded series")
        browse = QPushButton("...", objectName="small")
        browse.clicked.connect(self._browse)
        fid_row = QHBoxLayout()
        fid_row.addWidget(self.fid, 1)
        fid_row.addWidget(browse)
        self.mode = QComboBox()
        self.mode.addItem("blind (no structure)", "blind")
        self.mode.addItem("known: current structure", "structure")
        self.labeling = QComboBox()
        for v in ("natural", "15N", "2H-exchange", "unknown"):
            self.labeling.addItem(v)
        self.workers = QSpinBox(minimum=1, maximum=os.cpu_count() or 8, value=4)
        form.addRow("FID", fid_row)
        form.addRow("mode", self.mode)
        form.addRow("labelling", self.labeling)
        form.addRow("workers", self.workers)
        buttons = QHBoxLayout()
        self.run = QPushButton("Run analysis", objectName="primary")
        self.stop = QPushButton("Stop")
        self.report = QPushButton("Open report")
        buttons = QGridLayout()
        buttons.addWidget(self.run, 0, 0, 1, 2)
        buttons.addWidget(self.stop, 1, 0)
        buttons.addWidget(self.report, 1, 1)
        self.run.clicked.connect(self._run)
        self.stop.clicked.connect(lambda: self.session.blind_job and self.session.blind_job.stop())
        self.report.clicked.connect(self._open_report)
        left = QVBoxLayout()
        left.addLayout(form)
        left.addLayout(buttons)
        self.status = QLabel("no analysis yet", objectName="hint")
        self.status.setWordWrap(True)
        left.addWidget(self.status)
        hint = QLabel("Blind: processing per dataset, hypotheses from the line pattern, search with every exchange "
                      "regime, ranked on one scale (docs/WORKFLOW.md W5; skills/zulf-blind-analysis). The sampling "
                      "rate comes from scans.json or the .ini next to the FID. The ranking is a list of conditional "
                      "candidates, not an assignment.", objectName="hint")
        hint.setWordWrap(True)
        left.addWidget(hint)
        self.view = QTextBrowser()
        self.view.setOpenExternalLinks(True)
        self.view.setPlaceholderText("the report appears here")
        lay.addLayout(left)
        lay.addWidget(self.view, 1)
        self.shown = None

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(self, "Averaged FID", str(Path.home() / "research"), "FID (*.npy)")
        if path:
            self.fid.setText(path)

    def _run(self):
        workers = self.workers.value()
        if not confirm_workers(self, self.session, workers):
            return
        structure = self.session.spec if self.mode.currentData() == "structure" else None
        try:
            self.session.start_blind(fid=self.fid.text().strip() or None, workers=workers, structure=structure,
                                     labeling=self.labeling.currentText() if structure else "")
        except Exception as exc:
            QMessageBox.warning(self, "Analysis", f"{type(exc).__name__}: {exc}")

    def _open_report(self):
        st = self.session.blind_status()
        open_folder(st.get("report") or st.get("out_dir"))

    def refresh(self):
        st = self.session.blind_status()
        if not self.fid.text() and self.session.data and self.session.data.get("source_fid"):
            self.fid.setPlaceholderText(f"default: {self.session.data['source_fid']}")
        self.run.setEnabled(not st["running"])
        self.stop.setEnabled(st["running"])
        self.report.setEnabled(bool(st.get("out_dir")))
        if st.get("out_dir") is None:
            return
        self.status.setText(("running: " if st["running"] else
                             "finished: " if st["returncode"] == 0 else f"stopped / failed ({st['returncode']}): ")
                            + job_line(st) + f"\n{st['out_dir']}\n{st['last_line'][:160]}")
        report = st.get("report")
        if report and report != self.shown:
            self.shown = report
            self.view.setMarkdown(Path(report).read_text())
