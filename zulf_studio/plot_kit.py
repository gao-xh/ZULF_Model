"""Plot kit: looking at the spectrum in Studio's main plot, tied to the session's view (so a redraw, e.g. while
following a fit, keeps it). Replaces matplotlib's navigation toolbar, whose zoom and pan only changed the axes
until the next redraw.

- wheel: zoom the frequency axis around the cursor; Shift+wheel: shift it
- Zoom (box): drag a box: its frequency range becomes the view, its signal range the y range (Auto y clears it)
- Pan: drag to move the view
- double click, Reset: the whole spectrum; Back / Forward: earlier views
- readout: frequency and signal under the cursor and the nearest transition (isotopologue, frequency, relative)
- Save: the plot as PNG, SVG or PDF; Copy: the plot image to the clipboard
"""
from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QFileDialog, QLabel, QToolButton, QWidget

from .flow import FlowLayout, flow_policy


class PlotKit(QWidget):
    """Tool row of the main plot. window: the StudioWindow (session, fig, canvas, redraw, last lines)."""

    def __init__(self, window, parent=None):
        super().__init__(parent)
        self.w = window
        self.ylim = None                               # fixed y range of the signal axes (None: automatic)
        self.history, self.future = [], []
        self.mode = None                               # None (wheel and readout only), "zoom" or "pan"
        self._drag = None
        self._band = None
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
        self.b_zoom = button("Zoom", "drag a box to zoom in (frequency and signal); wheel: zoom around the cursor",
                             lambda on: self._set_mode("zoom" if on else None), checkable=True)
        self.b_pan = button("Pan", "drag to move the view; Shift+wheel also moves it",
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
        self.b_zoom.setChecked(mode == "zoom")
        self.b_pan.setChecked(mode == "pan")
        self.w.canvas.setCursor(Qt.CrossCursor if mode == "zoom" else Qt.OpenHandCursor if mode == "pan"
                                else Qt.ArrowCursor)

    # ---- mouse -----------------------------------------------------------------------------
    def _axes(self):
        """The axes that share the frequency axis (signal, residual, lines)."""
        return list(self.w.fig.axes)

    def _signal_axes(self):
        axes = self._axes()
        return axes[0] if axes else None

    def _scroll(self, e):
        if e.inaxes not in self._axes() or e.xdata is None:
            return
        lo, hi = self._view()
        if e.key == "shift":                           # Shift+wheel: move by 10 % of the width per step
            d = 0.1 * (hi - lo) * (1 if e.button == "up" else -1)
            self.set_view(lo + d, hi + d)
            return
        f = 0.8 if e.button == "up" else 1.25
        x = e.xdata
        self.set_view(x - (x - lo) * f, x + (hi - x) * f)

    def _press(self, e):
        if e.inaxes not in self._axes() or e.xdata is None:
            return
        if e.dblclick:
            self.reset()
            return
        if e.button == 1 and self.mode in ("zoom", "pan"):
            self._drag = {"x": e.xdata, "y": e.ydata, "px": e.x, "ax": e.inaxes, "view": self._view()}
            if self.mode == "pan":
                self.w.canvas.setCursor(Qt.ClosedHandCursor)

    def _motion(self, e):
        self._show_readout(e)
        d = self._drag
        if d is None or e.xdata is None:
            return
        if self.mode == "pan":
            lo, hi = d["view"]
            ax = d["ax"]
            width_px = max(ax.bbox.width, 1.0)
            shift = (e.x - d["px"]) / width_px * (hi - lo)
            self.w._guard(self.w.session.set_view, lo - shift, hi - shift)
        elif self.mode == "zoom":
            if self._band is not None:
                self._band.remove()
            t = self.w.t
            self._band = d["ax"].axvspan(min(d["x"], e.xdata), max(d["x"], e.xdata), color=t["accent"], alpha=0.15,
                                         lw=0)
            self.w.canvas.draw_idle()

    def _release(self, e):
        d, self._drag = self._drag, None
        if d is None:
            return
        if self._band is not None:
            self._band.remove()
            self._band = None
        if self.mode == "pan":
            self.history.append((d["view"], self.ylim))
            self.future.clear()
            self._update_buttons()
            self.w.canvas.setCursor(Qt.OpenHandCursor)
            return
        if e.xdata is None or abs(e.x - d["px"]) < 4:            # a click, not a drag
            self.w.canvas.draw_idle()
            return
        ylim = "keep"
        if d["ax"] is self._signal_axes() and e.ydata is not None and d["y"] is not None:
            y0, y1 = sorted((d["y"], e.ydata))
            span = np.diff(d["ax"].get_ylim())[0]
            if y1 - y0 > 0.02 * abs(span):            # a real box, not a horizontal stroke: also the y range
                ylim = (y0, y1)
        self.set_view(d["x"], e.xdata, ylim=ylim)

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
