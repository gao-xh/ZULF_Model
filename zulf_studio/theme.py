"""Look of ZULF Studio: colour tokens (light and dark), the Qt style sheet and the matching matplotlib style.

The plot colours follow the publication figures (scripts/figure_tools.py): ink for data, purple for the
simulation, the isotopologue palette for line sticks.
"""
from __future__ import annotations

ISOTOPOLOGUE = ("#c0469e", "#2f8f5b", "#2a78d6", "#eb6834", "#7a7a7a", "#8e5ac8")

LIGHT = {
    "bg": "#f4f4f1", "panel": "#ffffff", "panel2": "#fafaf8", "ink": "#1d1d1f", "muted": "#6b6b70",
    "line": "#e3e3df", "line2": "#d2d2cd", "accent": "#2a78d6", "accent_ink": "#ffffff", "accent_soft": "#e4eefb",
    "good": "#2a8a52", "bad": "#c0392b", "data": "#1d1d1f", "sim": "#7b2f8f", "trace": "#2f8f5b",
    "resid": "#8a8a8f", "grid": "#ecece8", "groove": "#e3e3df", "handle": "#ffffff",
}
DARK = {
    "bg": "#161618", "panel": "#1f1f22", "panel2": "#242428", "ink": "#ececee", "muted": "#9a9aa0",
    "line": "#323238", "line2": "#3d3d44", "accent": "#5b9cf0", "accent_ink": "#0d1420", "accent_soft": "#1c2a3d",
    "good": "#4cc27e", "bad": "#ff6b5e", "data": "#ececee", "sim": "#c88be0", "trace": "#4cc27e",
    "resid": "#8a8a92", "grid": "#2a2a30", "groove": "#3a3a42", "handle": "#ececee",
}


def _chevrons(color: str):
    """Small SVG chevrons (down, up) in `color` for combo and spin boxes; written once to the temp directory."""
    import tempfile
    from pathlib import Path
    d = Path(tempfile.gettempdir()) / "zulf_studio_theme"
    d.mkdir(exist_ok=True)
    paths = []
    for name, pts in (("down", "1,1.5 5,5.5 9,1.5"), ("up", "1,5.5 5,1.5 9,5.5")):
        f = d / f"chevron_{name}_{color.lstrip('#')}.svg"
        if not f.exists():
            f.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="10" height="7" viewBox="0 0 10 7">'
                         f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="1.6" '
                         f'stroke-linecap="round" stroke-linejoin="round"/></svg>')
        paths.append(f.as_posix())
    return paths


def stylesheet(t: dict) -> str:
    down, up = _chevrons(t["muted"])
    return f"""
* {{ font-size: 13px; }}
QMainWindow, QWidget#root {{ background: {t['bg']}; }}
QWidget {{ color: {t['ink']}; }}
QScrollArea, QScrollArea > QWidget > QWidget {{ background: transparent; border: none; }}
QLabel#title {{ font-size: 17px; font-weight: 600; }}
QLabel#subtitle {{ color: {t['muted']}; }}
QLabel#badge {{ background: {t['accent_soft']}; color: {t['accent']}; border-radius: 10px; padding: 3px 10px;
               font-weight: 600; }}
QLabel#hint {{ color: {t['muted']}; font-size: 11px; }}
QLabel#fine {{ color: {t['muted']}; font-size: 11px; }}
QWidget#header {{ background: {t['panel']}; border-bottom: 1px solid {t['line']}; }}
QGroupBox {{ background: {t['panel']}; border: 1px solid {t['line']}; border-radius: 10px; margin-top: 14px;
            padding: 12px 10px 8px 10px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 12px; top: 2px; padding: 0 4px; color: {t['muted']};
                   font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 1px; }}
QPushButton {{ background: {t['panel2']}; border: 1px solid {t['line2']}; border-radius: 7px; padding: 5px 12px; }}
QPushButton:hover {{ border-color: {t['accent']}; }}
QPushButton:pressed {{ background: {t['accent_soft']}; }}
QPushButton:disabled {{ color: {t['muted']}; border-color: {t['line']}; }}
QPushButton#primary {{ background: {t['accent']}; color: {t['accent_ink']}; border-color: {t['accent']};
                       font-weight: 600; }}
QPushButton#primary:hover {{ background: {t['accent']}; border-color: {t['ink']}; }}
QPushButton:checked {{ background: {t['accent_soft']}; border-color: {t['accent']}; color: {t['accent']}; }}
QLineEdit, QPlainTextEdit, QTextBrowser, QTableWidget {{
    background: {t['panel2']}; border: 1px solid {t['line2']}; border-radius: 7px; padding: 3px 6px;
    selection-background-color: {t['accent']}; selection-color: {t['accent_ink']}; }}
QDoubleSpinBox, QSpinBox, QComboBox {{ background: {t['panel2']}; border: 1px solid {t['line2']};
    border-radius: 7px; padding: 2px 4px 2px 6px; min-height: 20px; }}
QAbstractSpinBox::up-button, QAbstractSpinBox::down-button {{ subcontrol-origin: border; width: 16px;
    border: none; background: transparent; }}
QAbstractSpinBox::up-button {{ subcontrol-position: top right; }}
QAbstractSpinBox::down-button {{ subcontrol-position: bottom right; }}
QAbstractSpinBox::up-arrow {{ image: url({up}); width: 8px; height: 6px; }}
QAbstractSpinBox::down-arrow {{ image: url({down}); width: 8px; height: 6px; }}
QComboBox::drop-down {{ subcontrol-origin: padding; subcontrol-position: center right; width: 20px; border: none; }}
QComboBox::down-arrow {{ image: url({down}); width: 10px; height: 7px; }}
QLineEdit:focus, QPlainTextEdit:focus, QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {{
    border-color: {t['accent']}; }}
QPlainTextEdit#mono, QLineEdit#mono {{ font-family: "SF Mono", Menlo, monospace; font-size: 12px; }}
QComboBox QAbstractItemView {{ background: {t['panel']}; border: 1px solid {t['line2']};
                              selection-background-color: {t['accent_soft']}; selection-color: {t['ink']}; }}
QSlider::groove:horizontal {{ height: 4px; background: {t['groove']}; border-radius: 2px; }}
QSlider::sub-page:horizontal {{ background: {t['accent']}; border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {t['handle']}; border: 1px solid {t['line2']}; width: 14px; height: 14px;
                             margin: -6px 0; border-radius: 8px; }}
QSlider::handle:horizontal:hover {{ border-color: {t['accent']}; }}
QSlider#fineSlider::groove:horizontal {{ height: 2px; }}
QSlider#fineSlider::sub-page:horizontal {{ background: {t['groove']}; }}
QSlider#fineSlider::handle:horizontal {{ width: 10px; height: 10px; margin: -5px 0; border-radius: 6px; }}
QTabWidget::pane {{ background: {t['panel']}; border: 1px solid {t['line']}; border-radius: 10px; top: -1px; }}
QTabBar::tab {{ background: transparent; border: none; padding: 7px 14px; color: {t['muted']}; }}
QTabBar::tab:selected {{ color: {t['ink']}; border-bottom: 2px solid {t['accent']}; font-weight: 600; }}
QTabBar::tab:hover {{ color: {t['ink']}; }}
QHeaderView::section {{ background: {t['panel2']}; border: none; border-bottom: 1px solid {t['line']};
                       padding: 4px 6px; color: {t['muted']}; }}
QTableWidget {{ gridline-color: {t['line']}; alternate-background-color: {t['panel']}; }}
QTableWidget::item {{ padding: 0 6px; }}
QTableWidget::item:selected {{ background: {t['accent_soft']}; color: {t['ink']}; }}
QCheckBox {{ spacing: 6px; }}
QStatusBar {{ background: {t['panel']}; border-top: 1px solid {t['line']}; color: {t['muted']}; }}
QStatusBar::item {{ border: none; }}
QProgressBar {{ background: {t['groove']}; border: none; border-radius: 3px; max-height: 6px; min-width: 120px; }}
QProgressBar::chunk {{ background: {t['accent']}; border-radius: 3px; }}
QLabel#busy {{ color: {t['accent']}; font-weight: 600; }}
QLabel#idle {{ color: {t['muted']}; }}
QLabel#warn {{ color: {t['bad']}; font-weight: 600; }}
QLabel#ok {{ color: {t['good']}; }}
QLabel#pill {{ background: {t['accent_soft']}; color: {t['accent']}; border-radius: 10px; padding: 3px 10px;
              font-weight: 600; }}
QLabel#pillwarn {{ background: {t['panel2']}; color: {t['bad']}; border: 1px solid {t['bad']}; border-radius: 10px;
                  padding: 2px 9px; font-weight: 600; }}
QPushButton#small {{ padding: 2px 9px; border-radius: 6px; font-size: 12px; }}
QWidget#section {{ background: {t['panel']}; border: 1px solid {t['line']}; border-radius: 10px; }}
QWidget#sectionHead, QWidget#sectionBody {{ background: transparent; border: none; }}
QToolButton#sectionToggle {{ border: none; background: transparent; color: {t['ink']}; font-size: 12px;
                            font-weight: 700; letter-spacing: 0.5px; padding: 0; }}
QToolButton#kit {{ background: {t['panel2']}; border: 1px solid {t['line']}; border-radius: 6px; padding: 2px 7px;
                  color: {t['ink']}; }}
QToolButton#kit:hover {{ border-color: {t['line2']}; }}
QToolButton#kit:checked {{ background: {t['accent']}; color: {t['accent_ink']}; border-color: {t['accent']}; }}
QToolButton#kit:disabled {{ color: {t['muted']}; }}
QWidget#followBar {{ background: {t['accent_soft']}; border: 1px solid {t['accent']}; border-radius: 9px; }}
QWidget#plotCard {{ background: {t['panel']}; border: 1px solid {t['line']}; border-radius: 10px; }}
QLabel#dataInfo {{ color: {t['ink']}; }}
QWidget#modeBar {{ background: {t['panel2']}; border: 1px solid {t['line']}; border-radius: 9px; }}
QPushButton#mode {{ background: transparent; border: none; border-radius: 7px; padding: 4px 14px; color: {t['muted']};
                   font-weight: 600; }}
QPushButton#mode:hover {{ color: {t['ink']}; }}
QPushButton#mode:checked {{ background: {t['accent']}; color: {t['accent_ink']}; }}
QSplitter::handle {{ background: {t['bg']}; }}
QToolBar {{ background: transparent; border: none; }}
QMenuBar {{ background: {t['panel']}; }}
QScrollBar:vertical {{ background: transparent; width: 10px; }}
QScrollBar::handle:vertical {{ background: {t['line2']}; border-radius: 4px; min-height: 30px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; width: 0; }}
"""


def matplotlib_style(t: dict) -> dict:
    return {"figure.facecolor": t["panel"], "axes.facecolor": t["panel"], "axes.edgecolor": t["line2"],
            "axes.labelcolor": t["muted"], "axes.titlecolor": t["ink"], "xtick.color": t["muted"],
            "ytick.color": t["muted"], "text.color": t["ink"], "axes.spines.top": False,
            "axes.spines.right": False, "axes.grid": True, "grid.color": t["grid"], "grid.linewidth": 0.6,
            "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.frameon": False,
            "legend.fontsize": 8.5, "lines.solid_capstyle": "round"}


def pixel_margins(fig, canvas, nrows=1, left=80, right=14, top=12, bottom=40, gap=8):
    """Axes margins fixed in pixels (axis labels, titles), recomputed for the canvas size. Used instead of
    matplotlib's constrained layout, which leaves the axes collapsed once a canvas was very small (a narrow
    drawer). gap: pixels between stacked rows (titles and tick labels of the upper row)."""
    w, h = max(canvas.width(), 50), max(canvas.height(), 50)
    inner = max(h - top - bottom - gap * (nrows - 1), 1)
    fig.subplots_adjust(left=min(left / w, 0.45), right=max(1 - right / w, 0.55), top=max(1 - top / h, 0.55),
                        bottom=min(bottom / h, 0.4), hspace=min(gap * nrows / inner, 0.9) if nrows > 1 else 0.2)
