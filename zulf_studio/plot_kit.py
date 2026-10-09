"""Plot kit: looking at the spectrum in Studio's main plot, tied to the session's view (so a redraw, e.g. while
following a fit, keeps it). Replaces matplotlib's navigation toolbar, whose zoom and pan only changed the axes
until the next redraw.

- left drag: a box (drawn by Qt over the plot, so redraws while following a fit do not wipe it): its frequency
  range becomes the view and, if it is tall enough on the signal axes, its signal range the y range (Auto y
  clears it); a horizontal stroke zooms the frequency only
- right drag, or left drag with Pan on: move the view
- wheel: zoom the frequency axis around the cursor; Shift+wheel: shift it; Ctrl+wheel: zoom the signal axis
- double click, Reset: the whole spectrum; Back / Forward: earlier views
- readout: frequency and signal under the cursor and the nearest transition (isotopologue, frequency, relative)
- Save: the plot as PNG, SVG or PDF; Copy: the plot image to the clipboard
"""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QRubberBand, QToolButton, QWidget

from .flow import FlowLayout, flow_policy


class PlotKit(QWidget):
    """Tool row of the main plot. window: the StudioWindow (session, fig, canvas, redraw, last lines)."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.w = window
        self.ylim = None                               # fixed y range of the signal axes (None: automatic)
        self.history, self.future = [], []
        self.mode = None                               # None (left drag zooms) or "pan" (left drag moves)
        self._drag = None
        self._band = QRubberBand(QRubberBand.Rectangle, self.w.canvas)
        flow_policy(self)
        lay = FlowLayout(self, spacing=6)              # wraps on narrow windows (never widens the window)

        def button(text, tip, slot, checkable=False):
            b = QToolButton(text=text, checkable=checkable, objectName="kit")
            b.setToolTip(tip)
            b.clicked.connect(slot)
            lay.addWidget(b)
            return b
        self.b_reset = button("Reset", "the whole spectrum (also: double click on the plot)", self.reset)
        self.b_back = button("Back", "the previous view", self.back)
        self.b_fwd = button("Fwd", "the next view", self.forward)
        self.b_pan = button("Pan", "left drag moves the view instead of drawing a zoom box (right drag always "
                            "moves it; Shift+wheel too)",
                            lambda on: self._set_mode("pan" if on else None), checkable=True)
        self.b_auto = button("Auto y", "automatic signal range (clears a box zoom's y range)", self.auto_y)
        self.b_save = button("Save", "save the plot as PNG, SVG or PDF", self.save)
        self.b_copy = button("Copy", "copy the plot image to the clipboard", self.copy)
        self.readout = QLabel(" ", objectName="hint")  # placed under the plot by the window
        self.readout.setMinimumHeight(16)
        c = self.w.canvas
        c.mpl_connect("scroll_event", self._scroll)
        c.mpl_connect("button_press_event", self._press)
        c.mpl_connect("motion_notify_event", self._motion)
        c.mpl_connect("button_release_event", self._release)
        self._update_buttons()

    # ---- views -----------------------------------------------------------------------------
    def _view(self):
        return list(self.w.session.view or self.w.session.full_view())

    def set_view(self, lo, hi, ylim="keep", remember=True):
        lo, hi = float(min(lo, hi)), float(max(lo, hi))
        if hi - lo < 1e-3:
            return
        if remember:
            self.history.append((self._view(), self.ylim))
            del self.history[:-50]
            self.future.clear()
        if ylim != "keep":
            self.ylim = ylim
        self.w._guard(self.w.session.set_view, lo, hi)
        self._update_buttons()

    def reset(self):
        self.set_view(*self.w.session.full_view(), ylim=None)

    def back(self):
        if self.history:
            self.future.append((self._view(), self.ylim))
            view, self.ylim = self.history.pop()
            self.set_view(*view, remember=False)

    def forward(self):
        if self.future:
            self.history.append((self._view(), self.ylim))
            view, self.ylim = self.future.pop()
            self.set_view(*view, remember=False)

    def auto_y(self):
        self.ylim = None
        self.w.schedule()

    def _update_buttons(self):
        self.b_back.setEnabled(bool(self.history))
        self.b_fwd.setEnabled(bool(self.future))

    def _set_mode(self, mode):
        self.mode = mode
        self.b_pan.setChecked(mode == "pan")
        self.w.canvas.setCursor(Qt.OpenHandCursor if mode == "pan" else Qt.ArrowCursor)

    # ---- mouse -----------------------------------------------------------------------------
    def _axes(self):
        """The axes that share the frequency axis (signal, residual, lines)."""
        return list(self.w.fig.axes)

    def _signal_axes(self):
        axes = self._axes()
        return axes[0] if axes else None

    @staticmethod
    def _qpos(e):
        """Widget position (logical pixels) of a matplotlib mouse event."""
        g = getattr(e, "guiEvent", None)
        if g is not None:
            pos = g.position() if hasattr(g, "position") else g.pos()
            return QPoint(int(pos.x()), int(pos.y()))
        return None

    def _scroll(self, e):
        if e.inaxes not in self._axes() or e.xdata is None:
            return
        lo, hi = self._view()
        up = e.button == "up"
        if e.key in ("control", "ctrl", "cmd", "super") and e.inaxes is self._signal_axes() and e.ydata is not None:
            y0, y1 = e.inaxes.get_ylim()                 # Ctrl+wheel: the signal axis around the cursor
            f = 0.8 if up else 1.25
            self.ylim = (e.ydata - (e.ydata - y0) * f, e.ydata + (y1 - e.ydata) * f)
            self.w.schedule()
            return
        if e.key == "shift":                           # Shift+wheel: move by 10 % of the width per step
            d = 0.1 * (hi - lo) * (1 if up else -1)
            self.set_view(lo + d, hi + d)
            return
        f = 0.8 if up else 1.25
        x = e.xdata
        self.set_view(x - (x - lo) * f, x + (hi - x) * f)

    def _press(self, e):
        if e.inaxes not in self._axes() or e.xdata is None:
            return
        if e.dblclick:
            self.reset()
            return
        pan = e.button == 3 or (e.button == 1 and self.mode == "pan")
        if e.button not in (1, 3):
            return
        self._drag = {"x": e.xdata, "y": e.ydata, "px": e.x, "py": e.y, "ax": e.inaxes, "view": self._view(),
                      "pan": pan, "q": self._qpos(e)}
        if pan:
            self.w.canvas.setCursor(Qt.ClosedHandCursor)
        elif self._drag["q"] is not None:
            self._band.setGeometry(QRect(self._drag["q"], QSize()))
            self._band.show()

    def _motion(self, e):
        self._show_readout(e)
        d = self._drag
        if d is None:
            return
        if d["pan"]:
            lo, hi = d["view"]
            width_px = max(d["ax"].bbox.width, 1.0)
            shift = (e.x - d["px"]) / width_px * (hi - lo)
            self.w._guard(self.w.session.set_view, lo - shift, hi - shift)
            return
        q = self._qpos(e)
        if q is not None and d["q"] is not None:
            box = QRect(d["q"], q).normalized()
            if abs(e.y - d["py"]) < 10 * self.w.canvas.devicePixelRatioF():   # a stroke: the full height
                ax = d["ax"]
                h = self.w.canvas.height()
                dpr = self.w.canvas.devicePixelRatioF()
                top, bottom = h - ax.bbox.y1 / dpr, h - ax.bbox.y0 / dpr
                box = QRect(QPoint(box.left(), int(top)), QPoint(box.right(), int(bottom)))
            self._band.setGeometry(box)

    def _release(self, e):
        d, self._drag = self._drag, None
        self._band.hide()
        if d is None:
            return
        if d["pan"]:
            self.history.append((d["view"], self.ylim))
            self.future.clear()
            self._update_buttons()
            self.w.canvas.setCursor(Qt.OpenHandCursor if self.mode == "pan" else Qt.ArrowCursor)
            return
        if e.x is None or abs(e.x - d["px"]) < 5:      # a click, not a drag
            return
        ax = d["ax"]
        x1, y1 = ax.transData.inverted().transform((e.x, e.y))
        ylim = "keep"
        if ax is self._signal_axes() and d["y"] is not None and abs(e.y - d["py"]) >= 10 * \
                self.w.canvas.devicePixelRatioF():
            ylim = tuple(sorted((d["y"], float(y1))))   # a real box: also the signal range
        self.set_view(d["x"], float(x1), ylim=ylim)

    def _show_readout(self, e):
        if e.inaxes not in self._axes() or e.xdata is None:
            self.readout.setText("")
            return
        text = f"{e.xdata:.3f} Hz"
        if e.inaxes is self._signal_axes() and e.ydata is not None:
            text += f"   signal {e.ydata:.4g}"
        lines = getattr(self.w, "_last_lines", None) or []
        if lines:
            near = min(lines, key=lambda r: abs(r["frequency_hz"] - e.xdata))
            lo, hi = self._view()
            if abs(near["frequency_hz"] - e.xdata) < 0.02 * (hi - lo):
                text += (f"   nearest line {near['component']} {near['frequency_hz']:.3f} Hz "
                         f"(relative {near['relative']:.3g})")
        self.readout.setText(text)

    # ---- output ----------------------------------------------------------------------------
    def save(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save the plot", "spectrum.png",
                                              "PNG (*.png);;SVG (*.svg);;PDF (*.pdf)")
        if path:
            self.w.fig.savefig(path, dpi=200, facecolor=self.w.fig.get_facecolor())
            self.w.statusBar().showMessage(f"plot saved: {path}", 8000)

    def copy(self):
        QApplication.clipboard().setPixmap(self.w.canvas.grab())
        self.w.statusBar().showMessage("plot copied to the clipboard", 5000)

    def apply_ylim(self, ax):
        """Called by the redraw on the signal axes."""
        if self.ylim is not None:
            ax.set_ylim(*self.ylim)
