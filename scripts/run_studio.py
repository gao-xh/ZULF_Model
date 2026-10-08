"""ZULF Studio: real-time zero-field spectrum simulator with fitting, a fit-progress slider, log, terminal, Python
console and a JSON API for AI agents.

    python scripts/run_studio.py                                   # acetonitrile, zero field
    python scripts/run_studio.py --series runs/series/acn/series.json --fit runs/processed/acn_field_fam2
    python scripts/run_studio.py --structure '{"motif": "N-ethyl (Et3N)", "one_bond": {"C1": 131, "C2": 125}}'
    python scripts/run_studio.py --no-gui --api-port 8766          # API only (agents)

Needs PySide6 (conda env zulf). The API listens on 127.0.0.1:8766 (/api/tools); see zulf_studio/api.py.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

if __name__ == "__main__":
    from zulf_studio.app import main
    main()
