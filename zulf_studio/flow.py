"""A layout that places its widgets in a row and wraps them onto further rows when the width runs out (the Qt
"flow layout" example), so tool rows never set a large minimum window width."""
from __future__ import annotations

from PySide6.QtCore import QPoint, QRect, QSize, Qt
from PySide6.QtWidgets import QLayout, QSizePolicy


class FlowLayout(QLayout):
    def __init__(self, parent=None, spacing=6, right_align_last=False):
        super().__init__(parent)
        self._items = []
        self._spacing = spacing
        self._right_last = right_align_last           # the last item sits at the right end of its row
        self.setContentsMargins(0, 0, 0, 0)

    def addItem(self, item):
        self._items.append(item)

    def count(self):
        return len(self._items)

    def itemAt(self, index):
        return self._items[index] if 0 <= index < len(self._items) else None

    def takeAt(self, index):
        return self._items.pop(index) if 0 <= index < len(self._items) else None

    def expandingDirections(self):
        return Qt.Orientation(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self._place(QRect(0, 0, width, 0), apply=False)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self._place(rect, apply=True)

    def sizeHint(self):
        w = sum(i.sizeHint().width() for i in self._items) + self._spacing * max(len(self._items) - 1, 0)
        h = max((i.sizeHint().height() for i in self._items), default=0)
        m = self.contentsMargins()
        return QSize(w + m.left() + m.right(), h + m.top() + m.bottom())

    def minimumSize(self):
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        m = self.contentsMargins()
        return size + QSize(m.left() + m.right(), m.top() + m.bottom())

    def _place(self, rect, apply):
        m = self.contentsMargins()
        area = rect.adjusted(m.left(), m.top(), -m.right(), -m.bottom())
        x, y, row_h = area.x(), area.y(), 0
        visible = [i for i in self._items if not (i.widget() and i.widget().isHidden())]
        for n, item in enumerate(visible):
            hint = item.sizeHint()
            if x + hint.width() > area.right() + 1 and x > area.x():
                x, y, row_h = area.x(), y + row_h + self._spacing, 0
            pos_x = x
            if self._right_last and n == len(visible) - 1:
                pos_x = max(x, area.right() + 1 - hint.width())
            if apply:
                item.setGeometry(QRect(QPoint(pos_x, y), hint))
            x = pos_x + hint.width() + self._spacing
            row_h = max(row_h, hint.height())
        return y + row_h - rect.y() + m.bottom()


def flow_policy(widget):
    """Let a container with a FlowLayout grow in height when its row wraps."""
    pol = widget.sizePolicy()
    pol.setHeightForWidth(True)
    pol.setVerticalPolicy(QSizePolicy.Minimum)
    widget.setSizePolicy(pol)
    return widget
