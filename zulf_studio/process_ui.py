"""Process mode of ZULF Studio (PLAN 8b): an averaged FID and its recipe with a live preview, the scans of its
run with keep checkboxes, and saving the processed spectrum (a series) with its recipe.

The recipe is applied by zulf_processing.series_spectrum, the same steps scripts/make_series_entry.py uses for a
saved series; slider moves are collected for 150 ms before the preview is recomputed. Re-averaging chosen scans
runs average_scans.py --keep as a job; nothing is written into the run folder.
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import (QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton, QSpinBox, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


class RecipeBox(QWidget):
    """Left card of Process: the source FID and the processing recipe (live)."""

    def __init__(self, session, window, slider_cls, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        self.pending = {}
        self.timer = QTimer(self, singleShot=True, interval=150)
        self.timer.timeout.connect(self._apply)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.source = QLabel("no FID: open an averaged FID or a scan folder", objectName="hint")
        self.source.setWordWrap(True)
        lay.addWidget(self.source)
        row = QHBoxLayout()
        for text, fn in (("FID ...", self._open_fid), ("Scan folder ...", self._open_scans)):
            b = QPushButton(text, objectName="small")
            b.clicked.connect(fn)
            row.addWidget(b)
        row.addStretch(1)
        lay.addLayout(row)
        self.crop = slider_cls("crop start", 0.0, 1.0, 0.1, 3, "s")
        self.record = slider_cls("record", 0.5, 30.0, 7.5, 2, "s")
        self.apod = slider_cls("window", 0.02, 5.0, 0.3, 3, "1/s", log=True)
        self.sg = slider_cls("drift filter", 0.0, 0.3, 0.0, 3, "s")
        self.sg.setToolTip("Savitzky-Golay drift filter length; 0 = the plan's default (201 samples at 4 kHz)")
        self.phase = slider_cls("phase0", -180.0, 180.0, 0.0, 1, "deg")
        self.delay = slider_cls("delay", -10.0, 2.0, -3.4, 3, "ms", fine=0.2)
        for key, w in (("crop_s", self.crop), ("record_s", self.record), ("apodization_per_s", self.apod),
                       ("sg_window_s", self.sg), ("phase0_deg", self.phase), ("delay_ms", self.delay)):
            w.value_changed.connect(lambda v, k=key: self._queue(k, float(v)))
            lay.addWidget(w)
        cal = QHBoxLayout()
        self.cal = QCheckBox("phase from the calibration and the switching edge", checked=True)
        self.cal.toggled.connect(self._calibration)
        cal.addWidget(self.cal)
        lay.addLayout(cal)
        form = QFormLayout()
        self.zf = QSpinBox(minimum=0, maximum=6, value=3)
        self.zf.valueChanged.connect(lambda v: self._queue("zero_fill", int(v)))
        self.grid_lo = QDoubleSpinBox(decimals=1, minimum=0.0, maximum=5000.0, value=20.0)
        self.grid_hi = QDoubleSpinBox(decimals=1, minimum=1.0, maximum=5000.0, value=380.0)
        grid = QHBoxLayout()
        grid.addWidget(self.grid_lo)
        grid.addWidget(QLabel("to"))
        grid.addWidget(self.grid_hi)
        for w in (self.grid_lo, self.grid_hi):
            w.valueChanged.connect(lambda _: self._queue("grid", [self.grid_lo.value(), self.grid_hi.value()]))
        self.exclude = QLineEdit("81.5,86")
        self.exclude.editingFinished.connect(lambda: self._queue("exclude", self.exclude.text().strip()))
        form.addRow("zero fill", self.zf)
        form.addRow("grid (Hz)", grid)
        form.addRow("not fitted (Hz)", self.exclude)
        lay.addLayout(form)
        self.info = QLabel(objectName="hint")
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        self._calibration(True)

    def _open_fid(self):
        path, _ = QFileDialog.getOpenFileName(self, "Averaged FID", self.window._start_dir(
            "fid", Path.home() / "research"), "Averaged FID (*.npy)")
        if path:
            self.window._remember_dir("fid", path)
            self.window._guard(self.session.set_process_source, path)

    def _open_scans(self):
        d = QFileDialog.getExistingDirectory(self, "Instrument run folder (n.dat, n.ini; read-only)",
                                             self.window._start_dir("scans", Path.home() / "research"))
        if d:
            self.window._remember_dir("scans", d)
            self.window._guard(self.session.import_scans, d, then="process")
            self.window.show_page(self.window.jobs_panel)

    def _calibration(self, on):
        for w in (self.phase, self.delay):
            w.setEnabled(not on)
        if self.session.process is not None:
            if on:
                self._queue("phase0_deg", None)
                self._queue("delay_ms", None)
            else:
                pv = self.session.process["preview"] or {}
                self._queue("phase0_deg", pv.get("phase0_deg", 0.0))
                self._queue("delay_ms", pv.get("delay_ms", -3.4))

    def _queue(self, key, value):
        if key in ("phase0_deg", "delay_ms") and self.cal.isChecked() and value is not None:
            return
        self.pending[key] = value
        self.timer.start()

    def _apply(self):
        if self.session.process is None or not self.pending:
            self.pending = {}
            return
        values, self.pending = self.pending, {}
        self.window._guard(self.session.set_recipe, **values)

    def refresh(self):
        p = self.session.process
        if p is None:
            return
        r = p["recipe"]
        st = self.session.process_status()
        self.source.setObjectName("dataInfo")
        self.source.setText(f"<b>{Path(p['fid']).parent.name}/{Path(p['fid']).name}</b><br>{len(p['y'])} points, "
                            f"{p['fs']:g} Hz ({Path(p['fs_source']).name}), record {len(p['y']) / p['fs']:.2f} s")
        for w, k in ((self.crop, "crop_s"), (self.record, "record_s"), (self.apod, "apodization_per_s"),
                     (self.sg, "sg_window_s")):
            if not w.spin.hasFocus():
                w.set_value(float(r[k] or 0.0))
        if self.cal.isChecked():
            self.phase.set_value(st["phase0_deg"] or 0.0)
            self.delay.set_value(st["delay_ms"] or 0.0)
        self.info.setText(f"switching edge {st['edge_ms']:.3f} ms; phase0 {st['phase0_deg']:.1f} deg, delay "
                          f"{st['delay_ms']:.3f} ms" + (" (calibration)" if self.cal.isChecked() else " (by hand)"))


class ScansPanel(QWidget):
    """The scans of the source run: per-scan deviation and late noise (average_scans.py), keep checkboxes, a
    deviation plot; re-average the chosen scans."""

    def __init__(self, session, window, theme, parent=None):
        super().__init__(parent)
        self.session, self.window, self.t = session, window, theme
        self.keep = None
        lay = QHBoxLayout(self)
        left = QVBoxLayout()
        self.fig = Figure(figsize=(5, 2.2), layout="constrained")
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.setMinimumHeight(150)
        left.addWidget(self.canvas, 1)
        ctl = QHBoxLayout()
        self.z = QDoubleSpinBox(decimals=1, minimum=0.5, maximum=100.0, value=5.0)
        by_z = QPushButton("keep z <=", objectName="small")
        by_z.clicked.connect(self._by_z)
        self.rng = QLineEdit(placeholderText="e.g. 0-5989")
        self.rng.setMaximumWidth(110)
        by_r = QPushButton("keep range", objectName="small")
        by_r.clicked.connect(self._by_range)
        all_ = QPushButton("all", objectName="small")
        all_.clicked.connect(lambda: self._set(np.ones(self._n(), bool)))
        for w in (by_z, self.z, self.rng, by_r, all_):
            ctl.addWidget(w)
        ctl.addStretch(1)
        left.addLayout(ctl)
        self.go = QPushButton("Average the kept scans", objectName="primary")
        self.go.clicked.connect(self._average)
        go_row = QHBoxLayout()
        go_row.addWidget(self.go)
        go_row.addStretch(1)
        left.addLayout(go_row)
        self.status = QLabel("no scan record: process an average made from a scan folder", objectName="hint")
        left.addWidget(self.status)
        lay.addLayout(left, 3)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["keep", "scan", "deviation", "z", "late noise"])
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(22)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.itemChanged.connect(self._item_changed)
        lay.addWidget(self.table, 2)
        self.shown = None

    def _n(self):
        return len(self.rows)

    @property
    def rows(self):
        p = self.session.process
        return (p or {}).get("scans", {}) and p["scans"].get("scans", []) or []

    def _set(self, keep):
        self.keep = np.asarray(keep, bool)
        self._fill_table()
        self._plot()

    def _by_z(self):
        z = np.array([r["deviation_z"] for r in self.rows])
        self._set(z <= self.z.value())

    def _by_range(self):
        import sys
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
        from average_scans import parse_scan_list
        try:
            chosen = parse_scan_list(self.rng.text())
        except ValueError:
            return
        ids = np.array([r["scan"] for r in self.rows])
        self._set(np.isin(ids, sorted(chosen)) & (self.keep if self.keep is not None else True))

    def _item_changed(self, item):
        if item.column() == 0 and self.keep is not None:
            self.keep[item.row()] = item.checkState() == Qt.Checked
            self._plot()

    def _average(self):
        if self.keep is None or not self.keep.any():
            return
        ids = [r["scan"] for r, k in zip(self.rows, self.keep) if k]
        try:
            self.session.average_selection(ids)
            self.window.show_page(self.window.jobs_panel)
        except Exception as exc:
            QMessageBox.warning(self, "Average", f"{type(exc).__name__}: {exc}")

    def _fill_table(self):
        rows = self.rows
        self.table.blockSignals(True)
        self.table.setRowCount(len(rows))
        for i, (r, k) in enumerate(zip(rows, self.keep)):
            c = QTableWidgetItem()
            c.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            c.setCheckState(Qt.Checked if k else Qt.Unchecked)
            self.table.setItem(i, 0, c)
            for j, v in enumerate((f"{r['scan']}", f"{r['deviation']:.1f}", f"{r['deviation_z']:.2f}",
                                   f"{r['late_noise']:.2f}")):
                self.table.setItem(i, j + 1, QTableWidgetItem(v))
        self.table.blockSignals(False)

    def _plot(self):
        rows = self.rows
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        t = self.t
        self.fig.set_facecolor(t["panel"])
        ax.set_facecolor(t["panel"])
        if rows:
            ids = np.array([r["scan"] for r in rows])
            dev = np.array([r["deviation"] for r in rows])
            ax.plot(ids[self.keep], dev[self.keep], ".", ms=2, color=t["data"], label="kept")
            ax.plot(ids[~self.keep], dev[~self.keep], "x", ms=3, color=t["bad"], label="left out")
            ax.set_yscale("log")
            ax.set_xlabel("scan", fontsize=8, color=t["muted"])
            ax.set_ylabel("deviation", fontsize=8, color=t["muted"])
            ax.tick_params(labelsize=7, colors=t["muted"])
            ax.legend(fontsize=7, frameon=False)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        self.canvas.draw_idle()
        n = int(self.keep.sum()) if self.keep is not None else 0
        self.status.setText(f"{n} of {len(rows)} scans kept" if rows else
                            "no scan record: process an average made from a scan folder")
        self.go.setEnabled(bool(rows) and n > 0)

    def refresh(self):
        p = self.session.process
        key = (p or {}).get("fid")
        if key == self.shown:
            return
        self.shown = key
        rows = self.rows
        self.keep = np.array([bool(r.get("kept", True)) for r in rows], bool)
        self._fill_table()
        self._plot()


class ProcessPanel(QWidget):
    """Right panel of Process: save the processed spectrum and its recipe; go on to Fit or Blind analysis."""

    def __init__(self, session, window, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        lay = QVBoxLayout(self)
        box = QGroupBox("Save the processed spectrum")
        bl = QFormLayout(box)
        self.name = QLineEdit(placeholderText="series name (default: folder and file name)")
        self.out = QLineEdit(placeholderText=str(session.workspace / "series" / "<name>"))
        bl.addRow("name", self.name)
        bl.addRow("folder", self.out)
        save = QPushButton("Save spectrum and recipe", objectName="primary")
        save.clicked.connect(self._save)
        bl.addRow(save)
        hint = QLabel("Writes series.json, the spectrum and recipe.json (make_series_entry.py with this recipe; fit "
                      "ranges: the whole grid minus power-line harmonics, instrument lines and the bands that are "
                      "not fitted). The spectrum is loaded for Fit when it is done.", objectName="hint")
        hint.setWordWrap(True)
        bl.addRow(hint)
        lay.addWidget(box)
        other = QGroupBox("Other sources")
        ol = QVBoxLayout(other)
        for text, fn in (("Import with the dialog (FID or scan folder) ...", window.import_fid_dialog),
                         ("Open a processed spectrum ...", window.open_dialog)):
            b = QPushButton(text)
            b.clicked.connect(fn)
            ol.addWidget(b)
        lay.addWidget(other)
        self.info = QLabel(objectName="hint")
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        go = QPushButton("Fit a model to the loaded spectrum  \u2192")
        go.clicked.connect(lambda: self.session.set_mode("fit"))
        blind = QPushButton("Blind analysis  \u2192")
        blind.clicked.connect(lambda: self.session.set_mode("blind"))
        lay.addWidget(go)
        lay.addWidget(blind)
        lay.addStretch(1)

    def _save(self):
        try:
            self.session.save_processed(label=self.name.text().strip(), out=self.out.text().strip() or None)
            self.window.show_page(self.window.jobs_panel)
        except Exception as exc:
            QMessageBox.warning(self, "Save", f"{type(exc).__name__}: {exc}")

    def refresh(self):
        d = self.session.data
        self.info.setText("loaded spectrum: none" if d is None else
                          f"loaded spectrum: {d['label']}, {len(d['freq'])} points, {len(d['ranges'])} fit ranges")
