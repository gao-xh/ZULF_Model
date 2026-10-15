"""Spin-system editor of ZULF Studio (PLAN 8c): components with isotopes, a J matrix (upper triangle) whose cells
are numbers or variable names, and weights. Apply builds the model; the variables and numeric couplings then have
sliders in the Couplings card. Text mode edits the matrix as JSON (paste from elsewhere)."""
from __future__ import annotations

import json

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                               QMessageBox, QPlainTextEdit, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout,
                               QWidget)

from zulf_hypothesis.spin_system import NAME


class SpinSystemEditor(QWidget):
    def __init__(self, session, window, parent=None):
        super().__init__(parent)
        self.session, self.window = session, window
        self.spec = None
        self.current = 0
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self.comp = QComboBox()
        self.comp.currentIndexChanged.connect(self._pick)
        add = QPushButton("+", objectName="small")
        add.setToolTip("add a component (another spin system with its own weight)")
        add.clicked.connect(self._add)
        rem = QPushButton("-", objectName="small")
        rem.setToolTip("remove this component")
        rem.clicked.connect(self._remove)
        row.addWidget(QLabel("component"))
        row.addWidget(self.comp, 1)
        row.addWidget(add)
        row.addWidget(rem)
        lay.addLayout(row)
        row2 = QHBoxLayout()
        self.name = QLineEdit()
        self.weight = QDoubleSpinBox(decimals=4, minimum=0.0, maximum=1e6, value=1.0)
        self.weight.setToolTip("weight of this component in the sum (a fit will fit amplitudes)")
        row2.addWidget(QLabel("name"))
        row2.addWidget(self.name, 1)
        row2.addWidget(QLabel("weight"))
        row2.addWidget(self.weight)
        lay.addLayout(row2)
        row3 = QHBoxLayout()
        self.isotopes = QLineEdit(placeholderText="e.g. 13C 1H 1H 1H")
        self.isotopes.setToolTip("isotopes in order, separated by spaces or commas (1H, 2H, 13C, 15N, 19F, 31P ...)")
        set_iso = QPushButton("Set", objectName="small")
        set_iso.setToolTip("resize the matrix to these isotopes (existing couplings are kept where they fit)")
        set_iso.clicked.connect(self._set_isotopes)
        row3.addWidget(QLabel("isotopes"))
        row3.addWidget(self.isotopes, 1)
        row3.addWidget(set_iso)
        lay.addLayout(row3)
        mode = QHBoxLayout()
        mode.addWidget(QLabel("J matrix (Hz or a variable name)", objectName="hint"))
        mode.addStretch(1)
        self.text_toggle = QPushButton("Text", objectName="small", checkable=True)
        self.text_toggle.setToolTip("edit the matrix as JSON text")
        self.text_toggle.toggled.connect(self._text_mode)
        mode.addWidget(self.text_toggle)
        lay.addLayout(mode)
        self.table = QTableWidget(0, 0)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)   # names readable
        self.table.horizontalHeader().setMinimumSectionSize(56)
        self.table.verticalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.setMinimumHeight(150)
        self.table.itemChanged.connect(self._cell_changed)
        lay.addWidget(self.table)
        self.text = QPlainTextEdit(objectName="mono")
        self.text.setMinimumHeight(150)
        self.text.hide()
        lay.addWidget(self.text)
        self.vars = QLabel(objectName="hint")
        self.vars.setWordWrap(True)
        lay.addWidget(self.vars)
        btn = QHBoxLayout()
        apply_ = QPushButton("Apply", objectName="primary")
        apply_.clicked.connect(self.apply)
        btn.addWidget(apply_)
        btn.addStretch(1)
        self.fit_weights = QCheckBox("fit weights", checked=True)
        self.fit_weights.setToolTip("a fit fits every component's amplitude; off: the weight ratios are held")
        btn.addWidget(self.fit_weights)
        lay.addLayout(btn)
        files = QHBoxLayout()
        for text, fn, tip in (("From structure", self._from_structure, "the current structure as an editable spin system"),
                              ("Load ...", self._load, "a model file (.json) or a ZULF_NMR_Suite molecule folder"),
                              ("Save ...", self._save, "save the model as JSON")):
            b = QPushButton(text, objectName="small")
            b.setToolTip(tip)
            b.clicked.connect(fn)
            files.addWidget(b)
        files.addStretch(1)
        lay.addLayout(files)

    # ---- state ---------------------------------------------------------------------------
    def refresh_values(self):
        """The variable values of the session (sliders, an applied fit) into the variables line; the cells being
        edited are left alone."""
        spec = self.session.spec
        if self.spec is None or "spin_system" not in spec or "spin_system" not in self.spec:
            return
        current = spec["spin_system"].get("variables", {})
        mine = self.spec["spin_system"].setdefault("variables", {})
        if any(mine.get(k) != v for k, v in current.items() if k in mine):
            mine.update({k: v for k, v in current.items() if k in mine})
            self._show_vars()

    def load_from_session(self):
        spec = self.session.spec
        if "spin_system" in spec:
            self.spec = json.loads(json.dumps(spec))
        elif self.spec is None:
            self.spec = {"compound": "custom", "spin_system": {"components": [
                {"name": "A", "isotopes": ["13C", "1H", "1H", "1H"],
                 "J": [[0, "J1", "J1", "J1"], ["J1", 0, 0, 0], ["J1", 0, 0, 0], ["J1", 0, 0, 0]], "weight": 1.0}],
                "variables": {"J1": 125.0}}}
        self.current = min(self.current, len(self._comps()) - 1)
        self.fit_weights.setChecked(not self.spec["spin_system"].get("fixed_weights", False))
        self._fill_components()

    def _comps(self):
        return self.spec["spin_system"]["components"]

    def _fill_components(self):
        self.comp.blockSignals(True)
        self.comp.clear()
        for c in self._comps():
            self.comp.addItem(f"{c.get('name', '')} ({len(c['isotopes'])} spins)")
        self.comp.setCurrentIndex(self.current)
        self.comp.blockSignals(False)
        self._show()

    def _store(self):
        """Name, weight and the matrix of the shown component into the specification."""
        if self.spec is None or not self._comps():
            return
        c = self._comps()[self.current]
        c["name"] = self.name.text().strip() or c.get("name", "A")
        c["weight"] = self.weight.value()
        if self.text_toggle.isChecked():
            c["J"] = json.loads(self.text.toPlainText())
        else:
            n = len(c["isotopes"])
            J = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    item = self.table.item(i, j)
                    t = item.text().strip() if item else ""
                    v = _value(t)
                    J[i][j] = J[j][i] = v
            c["J"] = J

    def _show(self):
        c = self._comps()[self.current]
        self.name.setText(str(c.get("name", "")))
        self.weight.setValue(float(c.get("weight", 1.0)))
        self.isotopes.setText(" ".join(c["isotopes"]))
        n = len(c["isotopes"])
        t = self.table
        t.blockSignals(True)
        t.clear()
        t.setRowCount(n)
        t.setColumnCount(n)
        heads = [f"{k + 1} {iso}" for k, iso in enumerate(c["isotopes"])]
        t.setHorizontalHeaderLabels(heads)
        t.setVerticalHeaderLabels(heads)
        J = c.get("J") or [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                v = J[min(i, j)][max(i, j)] if i != j else ""
                item = QTableWidgetItem("" if i == j else _text(v))
                item.setTextAlignment(Qt.AlignCenter)
                if j <= i:                              # lower triangle mirrors the upper one, diagonal unused
                    item.setFlags(Qt.ItemIsEnabled)
                    item.setForeground(QColor("#8a8a8f"))
                t.setItem(i, j, item)
        t.blockSignals(False)
        self.text.setPlainText(json.dumps(J))
        self._show_vars()

    def _show_vars(self):
        ss = self.spec["spin_system"]
        names = sorted({str(v) for c in ss["components"] for row in c["J"] for v in row
                        if isinstance(v, str) and NAME.match(str(v))})
        missing = [n for n in names if n not in ss.get("variables", {})]
        self.vars.setText(("variables: " + ", ".join(f"{n} = {ss['variables'][n]:g} Hz" for n in names
                                                     if n in ss.get("variables", {}))
                           if names else "no variables: type a name (e.g. J1) into cells that share one coupling")
                          + (f"<br>new (start at 0 Hz, set them with the sliders): {', '.join(missing)}"
                             if missing else ""))

    # ---- edits ---------------------------------------------------------------------------
    def _cell_changed(self, item):
        i, j = item.row(), item.column()
        if j > i:
            mirror = self.table.item(j, i)
            if mirror is not None:
                self.table.blockSignals(True)
                mirror.setText(item.text())
                self.table.blockSignals(False)
        self._store()
        self._show_vars()

    def _pick(self, index):
        if index < 0:
            return
        self._store()
        self.current = index
        self._show()

    def _add(self):
        self._store()
        self._comps().append({"name": f"S{len(self._comps()) + 1}", "isotopes": ["1H", "1H"], "J": [[0, 0], [0, 0]],
                              "weight": 1.0})
        self.current = len(self._comps()) - 1
        self._fill_components()

    def _remove(self):
        if len(self._comps()) <= 1:
            return
        self._comps().pop(self.current)
        self.current = max(self.current - 1, 0)
        self._fill_components()

    def _set_isotopes(self):
        self._store()
        iso = [x for x in self.isotopes.text().replace(",", " ").split() if x]
        if not iso:
            return
        c = self._comps()[self.current]
        old, n = c["J"], len(iso)
        c["J"] = [[old[i][j] if i < len(old) and j < len(old) else 0 for j in range(n)] for i in range(n)]
        c["isotopes"] = iso
        self._fill_components()

    def _text_mode(self, on):
        if on:
            self._store_table_only()
            self.text.setPlainText(json.dumps(self._comps()[self.current]["J"]))
        else:
            try:
                self._comps()[self.current]["J"] = json.loads(self.text.toPlainText())
            except json.JSONDecodeError as exc:
                QMessageBox.warning(self, "J matrix", f"not valid JSON: {exc}")
                self.text_toggle.blockSignals(True)
                self.text_toggle.setChecked(True)
                self.text_toggle.blockSignals(False)
                return
            self._show()
        self.table.setVisible(not on)
        self.text.setVisible(on)
        self._refit()

    def _refit(self):
        from PySide6.QtCore import QTimer
        from .app import fit_to_page
        stack = self.parentWidget()
        if stack is not None and hasattr(stack, "currentWidget"):
            QTimer.singleShot(0, lambda: fit_to_page(stack, True))

    def _store_table_only(self):
        was = self.text_toggle.isChecked()
        self.text_toggle.blockSignals(True)
        self.text_toggle.setChecked(False)
        self._store()
        self.text_toggle.setChecked(was)
        self.text_toggle.blockSignals(False)

    def apply(self):
        try:
            self._store()
        except json.JSONDecodeError as exc:
            QMessageBox.warning(self, "J matrix", f"not valid JSON: {exc}")
            return
        ss = self.spec["spin_system"]
        ss.setdefault("variables", {})
        for c in ss["components"]:
            for row in c["J"]:
                for v in row:
                    if isinstance(v, str) and NAME.match(v) and v not in ss["variables"]:
                        ss["variables"][v] = 0.0
        ss["fixed_weights"] = not self.fit_weights.isChecked()
        used = {v for c in ss["components"] for row in c["J"] for v in row if isinstance(v, str)}
        ss["variables"] = {k: v for k, v in ss["variables"].items() if k in used}
        self.window._guard(self.session.set_structure, json.loads(json.dumps(self.spec)))

    def _from_structure(self):
        if "spin_system" in self.session.spec:
            return
        self.window._guard(self.session.spin_system_from_structure)

    def _load(self):
        path, _ = QFileDialog.getOpenFileName(self, "Model (.json) or ZULF_NMR_Suite structure.csv",
                                              str(self.session.workspace), "Model (*.json *.csv)")
        if path:
            self.window._guard(self.session.load_model, path)

    def _save(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save model", str(self.session.workspace / "model.json"),
                                              "Model (*.json)")
        if path:
            self.window._guard(self.session.save_model, path)


def _value(text):
    if text == "":
        return 0
    try:
        return float(text)
    except ValueError:
        return text


def _text(v):
    if isinstance(v, (int, float)):
        return "" if float(v) == 0.0 else f"{float(v):g}"
    return str(v)
