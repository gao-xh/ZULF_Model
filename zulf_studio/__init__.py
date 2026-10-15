"""ZULF Studio: real-time zero-field spectrum simulator with fitting, a desktop window and an API for AI agents.

    python scripts/run_studio.py [--structure JSON] [--series runs/series/NAME/series.json] [--api-port 8766]

Layers (each usable without the next):
- `session.StudioSession`: the state (structure, couplings, field, line width, loaded spectrum, fit jobs, fit
  traces) and every operation on it; pure Python, no GUI.
- `api.StudioAPI` / `api.serve`: the same operations as JSON tools over local HTTP (discoverable at /api/tools),
  so an AI agent can drive the session the person sees.
- `app`: the PySide6 window (sliders, plot, line table, fit panel, log, terminal, Python console).
"""
from .session import StudioSession

__all__ = ["StudioSession"]
