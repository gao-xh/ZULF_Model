"""Per-run log written next to a script's outputs (RUN_LOG.md): command, code version, times, key results.

    from run_log import RunLog
    log = RunLog(out_dir, "fit_staged")          # at the start of main(), after parsing arguments
    ...
    log.finish({"residual": 0.171, "k_ex": 78.0}, figures=["fit.png"])

Every analysis of a spectrum also gets its own analysis log under docs/analysis/ (AGENTS.md); the run log is the
machine record of one execution that the analysis log points to.
"""
from __future__ import annotations

import json
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]


def _git(*args) -> str:
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, timeout=10).stdout.strip()
    except Exception:
        return ""


class RunLog:
    def __init__(self, out_dir, script: str, argv: Optional[Sequence[str]] = None):
        self.out = Path(out_dir)
        self.out.mkdir(parents=True, exist_ok=True)
        self.script = script
        self.argv = list(sys.argv if argv is None else argv)
        self.start = time.time()
        self.commit = _git("rev-parse", "--short", "HEAD")
        self.dirty = bool(_git("status", "--porcelain", "--untracked-files=no"))
        self._write(status="running")

    def _write(self, status: str, results: Optional[dict] = None, figures: Sequence[str] = (), notes: str = ""):
        lines = [f"# Run log: {self.script}", "",
                 f"- status: {status}",
                 f"- started: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime(self.start))}",
                 f"- code: {self.commit or 'unknown'}{' (uncommitted changes)' if self.dirty else ''}",
                 f"- output directory: {self.out}", "",
                 "## Command", "", "```", " ".join(shlex.quote(a) for a in self.argv), "```"]
        if status != "running":
            lines[4:4] = [f"- finished: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}",
                          f"- wall time: {round(time.time() - self.start)} s"]
        if results:
            lines += ["", "## Results", "", "```json", json.dumps(results, indent=1, default=str), "```"]
        if figures:
            lines += ["", "## Figures", ""] + [f"- {f}" for f in figures]
        if notes:
            lines += ["", "## Notes", "", notes]
        (self.out / "RUN_LOG.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    def finish(self, results: dict, figures: Sequence[str] = (), notes: str = ""):
        self._write("finished", results, figures, notes)
