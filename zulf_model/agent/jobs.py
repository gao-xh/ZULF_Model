"""Background jobs for long-running tools (training, generation, refinement, benchmarks).

A job is a child process running `python -m zulf_model.agent.worker JOB_ID`.
Requests, status and logs live in `<workspace>/jobs/<job_id>/` and survive
client reconnects, so an agent can poll `get_job` instead of blocking. Cancel
writes a flag and terminates the process if it does not stop by itself.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path

from .registry import workspace


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def job_dir(job_id: str) -> Path:
    if not job_id or any(c not in "0123456789abcdef" for c in job_id):
        raise ValueError("Invalid job id.")
    return workspace() / "jobs" / job_id


def _write(path: Path, value: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(value, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def submit(tool: str, arguments: dict) -> dict:
    job_id = uuid.uuid4().hex
    path = job_dir(job_id)
    path.mkdir(parents=True)
    _write(path / "request.json", {"tool": tool, "arguments": arguments})
    _write(path / "status.json", {"job_id": job_id, "tool": tool, "status": "queued", "created_utc": _now()})
    env = dict(os.environ)
    env.setdefault("ZULF_MODEL_WORKSPACE", str(workspace()))
    with (path / "worker.log").open("wb") as log:
        process = subprocess.Popen([sys.executable, "-m", "zulf_model.agent.worker", job_id], stdout=log, stderr=log,
                                   stdin=subprocess.DEVNULL, env=env,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    status = json.loads((path / "status.json").read_text(encoding="utf-8"))
    status["pid"] = process.pid
    _write(path / "status.json", status)
    return {"job_id": job_id, "status": "submitted", "pid": process.pid, "log_path": str(path / "worker.log"),
            "next": "Poll get_job with this job_id."}


def get(job_id: str) -> dict:
    path = job_dir(job_id)
    status = json.loads((path / "status.json").read_text(encoding="utf-8"))
    status["cancel_requested"] = (path / "cancel").exists()
    if (path / "result.json").exists():
        status["result"] = json.loads((path / "result.json").read_text(encoding="utf-8"))
    return status


def cancel(job_id: str, grace_s: float = 5.0) -> dict:
    path = job_dir(job_id)
    (path / "cancel").touch()
    status = json.loads((path / "status.json").read_text(encoding="utf-8"))
    pid = status.get("pid")
    if status.get("status") in ("queued", "running") and pid:
        deadline = time.time() + grace_s
        while time.time() < deadline:
            if get(job_id)["status"] not in ("queued", "running"):
                return get(job_id)
            time.sleep(0.2)
        try:
            os.kill(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            pass
        status.update(status="cancelled", finished_utc=_now())
        _write(path / "status.json", status)
    return get(job_id)


def list_jobs(limit: int = 20) -> dict:
    root = workspace() / "jobs"
    if not root.exists():
        return {"jobs": []}
    rows = []
    for path in sorted(root.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)[:limit]:
        try:
            rows.append(json.loads((path / "status.json").read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError):
            continue
    return {"jobs": rows}


def run_worker(job_id: str) -> None:
    from .tools import REGISTRY
    path = job_dir(job_id)
    request = json.loads((path / "request.json").read_text(encoding="utf-8"))
    status = json.loads((path / "status.json").read_text(encoding="utf-8"))
    status.update(status="running", started_utc=_now(), pid=os.getpid())
    _write(path / "status.json", status)
    try:
        arguments = dict(request["arguments"], _job_dir=str(path))
        result = REGISTRY.call(request["tool"], arguments)
        _write(path / "result.json", result)
        status.update(status="cancelled" if (path / "cancel").exists() else "complete", finished_utc=_now())
    except Exception as exc:  # recorded for the agent, never swallowed silently
        status.update(status="failed", error=f"{type(exc).__name__}: {exc}", finished_utc=_now())
    _write(path / "status.json", status)
