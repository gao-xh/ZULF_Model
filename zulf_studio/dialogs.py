"""Import and export dialogs of ZULF Studio.

Import: inspect first (StudioSession.inspect_path: scans, sampling rate and its source, record, sequence, time span,
cached averages of other tools, problems, a preview of the FID start and its spectrum), then choose how to import
(average the scans, or use a cached average) and the processing (record, crop, frequency grid, excluded bands).
The last settings are kept. Nothing is written into the data folder.

Export: one folder to share (StudioSession.export_bundle): the spectrum (data, simulation, residual; real,
imaginary, magnitude; CSV with a header and/or NPZ with complex arrays), the FID, the parameters and session, the
applied fit, the plot, and information.json with the sources and their sha256.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QRadioButton,
                               QTableWidget, QTableWidgetItem, QVBoxLayout, QHeaderView)

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure


def _style_axes(fig, t):
    fig.set_facecolor(t["panel"])
    for ax in fig.axes:
        ax.set_facecolor(t["panel"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(t["line2"])
        ax.tick_params(colors=t["muted"], labelsize=8)
        ax.xaxis.label.set_color(t["muted"])
        ax.title.set_color(t["ink"])


class ImportDialog(QDialog):
    """Inspect a scan folder or an averaged FID, then import it with chosen processing."""

    def __init__(self, session, path, settings, theme, parent=None):
        super().__init__(parent)
        self.session, self.path, self.settings = session, str(path), settings
        self.setWindowTitle(f"Import {Path(path).name}")
        from PySide6.QtWidgets import QApplication
        screen = QApplication.primaryScreen()
        avail = screen.availableGeometry() if screen is not None else None
        self.resize(min(980, avail.width() - 40) if avail else 980, min(860, avail.height() - 30) if avail else 860)
        info = session.inspect_path(path)
        self.info = info
        lay = QVBoxLayout(self)
        head = QLabel(f"<b>{info['kind']}</b> &nbsp; {self.path}")
        head.setWordWrap(True)
        lay.addWidget(head)
        top = QHBoxLayout()
        table = QTableWidget(len(info["rows"]), 2)
        table.horizontalHeader().setVisible(False)
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        for i, (k, v) in enumerate(info["rows"]):
            table.setItem(i, 0, QTableWidgetItem(k))
            table.setItem(i, 1, QTableWidgetItem(v))
        table.resizeColumnToContents(0)
        table.setMinimumHeight(250)
        top.addWidget(table, 1)
        if info["preview"]:
            fig = Figure(figsize=(4.6, 3.6), layout="constrained")
            canvas = FigureCanvasQTAgg(fig)
            canvas.setMinimumSize(420, 270)
            a1, a2 = fig.subplots(2, 1)
            pv = info["preview"]
            a1.plot(pv["t_ms"], pv["fid"], color=theme["data"], lw=0.7)
            a1.set_title("FID start", fontsize=9)
            a1.set_xlabel("ms", fontsize=8)
            a2.plot(pv["f"], pv["mag"], color=theme["sim"], lw=0.7)
            a2.set_title("magnitude spectrum (crop 0.1 s, 0.3 1/s)", fontsize=9)
            a2.set_xlabel("Hz", fontsize=8)
            _style_axes(fig, theme)
            top.addWidget(canvas, 1)
        lay.addLayout(top, 1)
        for w in info["warnings"]:
            lab = QLabel("\u26a0 " + w, objectName="warn")
            lab.setWordWrap(True)
            lay.addWidget(lab)
        saved = json.loads(settings.value("import/recipe", "{}") or "{}")
        self.source = QButtonGroup(self)
        if info["kind"] == "scan folder":
            box = QGroupBox("Source")
            bl = QVBoxLayout(box)
            self.r_avg = QRadioButton("average the scans (average_scans.py; per-scan screening, even / odd halves)")
            self.r_cached = QRadioButton("use the cached average halp_compiled.npy (another tool; scans unknown)")
            self.r_avg.setChecked(True)
            self.r_cached.setEnabled(bool(info.get("cached")))
            self.source.addButton(self.r_avg)
            self.source.addButton(self.r_cached)
            row = QHBoxLayout()
            self.exclude_z = QDoubleSpinBox(decimals=1, minimum=0.0, maximum=50.0, value=saved.get("exclude_z", 5.0))
            self.exclude_z.setToolTip("leave out scans whose deviation exceeds this robust z (0 keeps every scan)")
            row.addWidget(QLabel("exclude scans with robust z above"))
            row.addWidget(self.exclude_z)
            row.addWidget(QLabel("(0 = keep all)"))
            row.addStretch(1)
            bl.addWidget(self.r_avg)
            bl.addLayout(row)
            bl.addWidget(self.r_cached)
            out_row = QHBoxLayout()
            self.avg_out = QLineEdit(str(session.workspace / "averages" / Path(path).name))
            pick = QPushButton("...", objectName="small")
            pick.clicked.connect(lambda: self._pick_dir(self.avg_out))
            out_row.addWidget(QLabel("average to"))
            out_row.addWidget(self.avg_out, 1)
            out_row.addWidget(pick)
            bl.addLayout(out_row)
            hint = QLabel("The run folder is read only: nothing is written there. Averages kept for other work belong "
                          "in ~/research/<project>/data/processed/<measurement>/.", objectName="hint")
            hint.setWordWrap(True)
            bl.addWidget(hint)
            lay.addWidget(box)
        proc = QGroupBox("Processing (make_series_entry.py)")
        form = QFormLayout(proc)
        self.record = QDoubleSpinBox(decimals=2, minimum=0.2, maximum=100.0, value=saved.get("record_s", 7.5))
        self.crop = QDoubleSpinBox(decimals=3, minimum=0.0, maximum=5.0, value=saved.get("crop_s", 0.1),
                                   singleStep=0.01)
        self.grid = QLineEdit(saved.get("grid", "20,380"))
        self.exclude = QLineEdit(saved.get("exclude", "81.5,86"))
        self.rate = QDoubleSpinBox(decimals=3, minimum=0.0, maximum=100000.0, value=0.0)
        self.rate.setSpecialValueText("from scans.json / .ini")
        self.label = QLineEdit(placeholderText="series name (default: folder and file name)")
        form.addRow("record (s)", self.record)
        form.addRow("crop: record start (s)", self.crop)
        form.addRow("frequency grid lo,hi (Hz)", self.grid)
        form.addRow("left out of the fit ranges", self.exclude)
        form.addRow("sampling rate (Hz)", self.rate)
        form.addRow("name", self.label)
        tip = QLabel("Fit ranges: the whole grid minus the power-line harmonics and instrument lines (+-0.4 Hz) and "
                     "the excluded bands (the 83 Hz NMRduino feature by default).", objectName="hint")
        tip.setWordWrap(True)
        form.addRow(tip)
        lay.addWidget(proc)
        buttons = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.go = buttons.addButton("Import", QDialogButtonBox.AcceptRole)
        self.go.setObjectName("primary")
        self.go.setEnabled(info["kind"] in ("scan folder", "FID"))
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        self.result_status = None

    def _pick_dir(self, edit):
        d = QFileDialog.getExistingDirectory(self, "Folder", edit.text())
        if d:
            edit.setText(d)

    def _accept(self):
        recipe = {"record_s": self.record.value(), "crop_s": self.crop.value(), "grid": self.grid.text().strip(),
                  "exclude": self.exclude.text().strip()}
        if self.info["kind"] == "scan folder":
            recipe["exclude_z"] = self.exclude_z.value()
        self.settings.setValue("import/recipe", json.dumps(recipe))
        opts = dict(record_s=recipe["record_s"], crop_s=recipe["crop_s"], grid=recipe["grid"],
                    exclude=recipe["exclude"], label=self.label.text().strip(), sampling_rate=self.rate.value())
        try:
            if self.info["kind"] == "FID":
                self.result_status = self.session.import_fid(self.path, **opts)
            elif self.r_cached.isChecked():
                opts["label"] = opts["label"] or f"{Path(self.path).name}-halp"
                self.result_status = self.session.import_fid(self.info["cached"], **opts)
            else:
                self.result_status = self.session.import_scans(self.path, out=self.avg_out.text().strip() or None,
                                                               exclude_z=self.exclude_z.value(), fid_options=opts)
        except Exception as exc:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Import", f"{type(exc).__name__}: {exc}")
            return
        self.accept()


class ExportDialog(QDialog):
    """Choose what goes into an export folder, then StudioSession.export_bundle (and the plot image)."""

    def __init__(self, session, settings, figure, parent=None):
        super().__init__(parent)
        self.session, self.settings, self.figure = session, settings, figure
        self.setWindowTitle("Export")
        self.resize(640, 560)
        saved = json.loads(settings.value("export/choice", "{}") or "{}")
        lay = QVBoxLayout(self)
        spec = QGroupBox("Spectrum (columns frequency_hz, then real, imaginary, magnitude of each)")
        sl = QVBoxLayout(spec)
        has_data = session.data is not None
        self.c_data = QCheckBox("data (as loaded, with the display phase)", checked=saved.get("data", True))
        self.c_sim = QCheckBox("simulation (quick look: Lorentzian lines, one rate)", checked=saved.get("sim", True))
        self.c_res = QCheckBox("residual (data - simulation)", checked=saved.get("res", False))
        for c in (self.c_data, self.c_res):
            c.setEnabled(has_data)
        row = QHBoxLayout()
        self.r_full = QRadioButton("whole spectrum")
        self.r_view = QRadioButton("current view only")
        (self.r_view if saved.get("view_only") else self.r_full).setChecked(True)
        self.f_csv = QCheckBox("CSV", checked=saved.get("csv", True))
        self.f_npz = QCheckBox("NPZ (complex arrays)", checked=saved.get("npz", False))
        for w in (self.r_full, self.r_view, QLabel("   format"), self.f_csv, self.f_npz):
            row.addWidget(w)
        row.addStretch(1)
        for w in (self.c_data, self.c_sim, self.c_res):
            sl.addWidget(w)
        sl.addLayout(row)
        lay.addWidget(spec)
        more = QGroupBox("Also")
        ml = QVBoxLayout(more)
        src = session.data.get("source_fid") if has_data else None
        self.c_fid = QCheckBox("FID of the spectrum (fid.csv: time_s, signal)" + ("" if src else " - no source FID"),
                               checked=bool(src) and saved.get("fid", False))
        self.c_fid.setEnabled(bool(src))
        self.c_par = QCheckBox("parameters: parameters.json, couplings.csv, lines.csv, session file",
                               checked=saved.get("par", True))
        applied = session.applied is not None
        self.c_fit = QCheckBox("applied fit: fit.json, J_table.csv" + ("" if applied else " - no fit applied"),
                               checked=applied and saved.get("fit", True))
        self.c_fit.setEnabled(applied)
        self.c_plot = QCheckBox("plot as shown:", checked=saved.get("plot", True))
        prow = QHBoxLayout()
        self.p_png = QCheckBox("PNG", checked=saved.get("png", True))
        self.p_pdf = QCheckBox("PDF", checked=saved.get("pdf", False))
        self.p_svg = QCheckBox("SVG", checked=saved.get("svg", True))
        for w in (self.c_plot, self.p_png, self.p_pdf, self.p_svg):
            prow.addWidget(w)
        prow.addStretch(1)
        for w in (self.c_fid, self.c_par, self.c_fit):
            ml.addWidget(w)
        ml.addLayout(prow)
        info = QLabel("information.json records the sources (paths, sha256), the time, the code commit and the "
                      "settings, so the export can be traced.", objectName="hint")
        info.setWordWrap(True)
        ml.addWidget(info)
        lay.addWidget(more)
        dest = QGroupBox("Folder")
        dl = QFormLayout(dest)
        label = (session.data or {}).get("label") or session.spec.get("compound") or "studio"
        base = settings.value("export/dir", str(session.workspace / "exports"))
        self.dir = QLineEdit(base)
        pick = QPushButton("...", objectName="small")
        pick.clicked.connect(self._pick)
        drow = QHBoxLayout()
        drow.addWidget(self.dir, 1)
        drow.addWidget(pick)
        self.name = QLineEdit(f"{label}_{time.strftime('%Y%m%d-%H%M%S')}")
        dl.addRow("in", drow)
        dl.addRow("name", self.name)
        lay.addWidget(dest)
        buttons = QDialogButtonBox(QDialogButtonBox.Cancel)
        go = buttons.addButton("Export", QDialogButtonBox.AcceptRole)
        go.setObjectName("primary")
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        lay.addWidget(buttons)
        self.result = None

    def _pick(self):
        d = QFileDialog.getExistingDirectory(self, "Export into", self.dir.text())
        if d:
            self.dir.setText(d)

    def _accept(self):
        choice = {"data": self.c_data.isChecked(), "sim": self.c_sim.isChecked(), "res": self.c_res.isChecked(),
                  "view_only": self.r_view.isChecked(), "csv": self.f_csv.isChecked(), "npz": self.f_npz.isChecked(),
                  "fid": self.c_fid.isChecked(), "par": self.c_par.isChecked(), "fit": self.c_fit.isChecked(),
                  "plot": self.c_plot.isChecked(), "png": self.p_png.isChecked(), "pdf": self.p_pdf.isChecked(),
                  "svg": self.p_svg.isChecked()}
        self.settings.setValue("export/choice", json.dumps(choice))
        self.settings.setValue("export/dir", self.dir.text().strip())
        out = Path(self.dir.text().strip()).expanduser() / self.name.text().strip()
        spectrum = tuple(n for n, on in (("data", choice["data"]), ("simulation", choice["sim"]),
                                         ("residual", choice["res"])) if on)
        formats = tuple(f for f, on in (("csv", choice["csv"]), ("npz", choice["npz"])) if on) or ("csv",)
        try:
            self.result = self.session.export_bundle(str(out), spectrum=spectrum, view_only=choice["view_only"],
                                                     formats=formats, fid=choice["fid"], parameters=choice["par"],
                                                     fit=choice["fit"])
            if choice["plot"]:
                for ext in ("png", "pdf", "svg"):
                    if choice[ext]:
                        self.figure.savefig(out / f"plot.{ext}", dpi=200)
                        self.result["files"].append(f"plot.{ext}")
        except Exception as exc:
            from PySide6.QtWidgets import QMessageBox
            QMessageBox.warning(self, "Export", f"{type(exc).__name__}: {exc}")
            return
        self.accept()
