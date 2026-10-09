"""Live monitor of fit_joint_series runs: watch the objective of every start, its stage and its best couplings while
the fit runs, without changing the fit.

Writer side (used by fit_joint_series, on by default, `--monitor off` to disable): every residual evaluation of a
start appends its objective and its parameter vector ("x") to RUN/monitor/start_NNN.jsonl (buffered, flushed at most
every 0.5 s), so the spectrum at any evaluation can be drawn afterwards; the couplings and vector of the best point
so far are added when it improves (at most once a second, and at the end). RUN/monitor/status.json
holds the command, the phase, the number of starts and the scores of the finished starts; RUN/monitor/console.log
copies the console output. The writer only reads the objective value the fit has already computed; the optimizer,
its steps and its results are unchanged.

Viewer side:

    python scripts/fit_monitor.py RUN_OR_PARENT_DIR [--port 8770]     # browser page, http://127.0.0.1:8770
    python scripts/fit_monitor.py RUN_DIR --text [--interval 2]          # terminal table

The page lists every run with a monitor below the given directory (newest first) and shows the phase, the objective
of every start against its evaluations (best so far), the stage of each start, the couplings of the best point of a
start against its starting values, and the console output. "Show spectrum" renders data and model of the selected
start's best point (the problem is rebuilt from the recorded command in the viewer process, on demand only).
"""
import argparse
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
UI = Path(__file__).with_name("fit_monitor_ui.html")


# ------------------------------------------------------------------ writer side
class MonitorWriter:
    """Called by JointSeries.residual as monitor(label, cost, z) for one start (or stage) in its own process."""

    def __init__(self, directory, name, joint, keys, flush_s=0.5, snapshot_s=1.0):
        self.path = Path(directory) / f"{name}.jsonl"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(self.path, "a")
        self.joint, self.keys = joint, keys
        self.flush_s, self.snapshot_s = flush_s, snapshot_s
        self.n, self.best, self.best_z = 0, None, None
        self.t0 = self.last_flush = time.time()
        self.last_snapshot = 0.0
        self.snapshot_pending = False
        self.buffer = []
        self._write({"event": "start", "t": self.t0})

    def _couplings(self, z):
        return {key: self.joint.coupling_values(z, k).tolist() for k, key in enumerate(self.keys)}

    def _write(self, record):
        self.buffer.append(json.dumps(record))

    def __call__(self, label, cost, z):
        self.n += 1
        now = time.time()
        improved = self.best is None or cost < self.best
        if improved:
            self.best, self.best_z = cost, np.array(z, float)
            self.snapshot_pending = True
        rec = {"n": self.n, "t": round(now - self.t0, 3), "cost": cost, "best": self.best, "label": label,
               "x": [float(v) for v in z]}            # this evaluation's vector (Studio: spectrum at any point)
        if self.n == 1:
            rec["J0"] = self._couplings(z)                   # this start's own starting couplings
        if self.snapshot_pending and now - self.last_snapshot >= self.snapshot_s:
            rec["J"] = self._couplings(self.best_z)
            rec["z"] = self.best_z.tolist()
            self.last_snapshot, self.snapshot_pending = now, False
        self._write(rec)
        if now - self.last_flush >= self.flush_s:
            self.flush()

    def flush(self):
        if self.buffer:
            self.fh.write("\n".join(self.buffer) + "\n")
            self.fh.flush()
            self.buffer = []
        self.last_flush = time.time()

    def close(self, score=None):
        end = {"event": "end", "n": self.n, "t": round(time.time() - self.t0, 3), "best": self.best, "score": score}
        if self.best_z is not None:
            end["J"] = self._couplings(self.best_z)
            end["z"] = self.best_z.tolist()
        self._write(end)
        self.flush()
        self.fh.close()


class RunStatus:
    """RUN/monitor/status.json, rewritten atomically by the main process at every phase change."""

    def __init__(self, directory, argv, starts, keys, nodes, start_couplings):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        for old in self.dir.glob("*.jsonl"):            # a rerun into the same directory starts a fresh record
            old.unlink()
        self.data = {"argv": list(argv), "cwd": os.getcwd(), "started": time.time(), "phase": "starting",
                     "starts": int(starts), "keys": list(keys), "nodes": list(nodes), "finished": {},
                     "start_couplings": start_couplings, "pid": os.getpid()}
        self.write()

    def set(self, **kw):
        self.data.update(kw)
        self.write()

    def finished(self, k, score):
        self.data["finished"][str(k)] = score
        self.write()

    def write(self):
        self.data["updated"] = time.time()
        tmp = self.dir / "status.json.tmp"
        with open(tmp, "w") as fh:
            json.dump(self.data, fh)
        os.replace(tmp, self.dir / "status.json")


class Tee:
    """Console output also into RUN/monitor/console.log."""

    def __init__(self, stream, path):
        self.stream, self.fh = stream, open(path, "a")

    def write(self, s):
        self.stream.write(s)
        self.fh.write(s)
        self.fh.flush()

    def flush(self):
        self.stream.flush()
        self.fh.flush()

    def __getattr__(self, name):
        return getattr(self.stream, name)


# ------------------------------------------------------------------ reader side
def read_run(run, max_points=400):
    """Status and per-start histories (decimated to max_points, best so far) of one run directory."""
    mon = Path(run) / "monitor"
    with open(mon / "status.json") as fh:
        status = json.load(fh)
    starts = {}
    for p in sorted(mon.glob("*.jsonl")):
        recs = []
        with open(p) as fh:
            for line in fh:
                try:
                    recs.append(json.loads(line))
                except json.JSONDecodeError:        # a line being written
                    pass
        evals = [r for r in recs if "cost" in r]
        snaps = [r for r in recs if "J" in r]
        first = next((r for r in recs if "J0" in r), None)
        end = next((r for r in recs if r.get("event") == "end"), None)
        step = max(len(evals) // max_points, 1)
        hist = evals[::step]
        if evals and hist[-1] is not evals[-1]:
            hist.append(evals[-1])
        starts[p.stem] = {"evaluations": len(evals), "running": end is None,
                          "label": evals[-1]["label"] if evals else "", "best": evals[-1]["best"] if evals else None,
                          "cost": evals[-1]["cost"] if evals else None, "seconds": evals[-1]["t"] if evals else 0.0,
                          "score": end.get("score") if end else None,
                          "history": [[r["n"], r["cost"], r["best"]] for r in hist],
                          "J": snaps[-1]["J"] if snaps else None, "z": snaps[-1]["z"] if snaps else None,
                          "J0": first["J0"] if first else None}
    alive = None
    pid = status.get("pid")
    if pid:
        try:
            os.kill(int(pid), 0)
            alive = True
        except (OSError, ValueError):
            alive = False
    return {"status": status, "starts": starts, "alive": alive, "now": time.time()}


def find_runs(root):
    root = Path(root)
    found = [p.parent.parent for p in root.glob("**/monitor/status.json")]
    if (root / "monitor" / "status.json").exists() and root not in found:
        found.append(root)
    return sorted(set(found), key=lambda p: -(p / "monitor" / "status.json").stat().st_mtime)


def text_view(run, interval=2.0):
    try:
        while True:
            d = read_run(run)
            st = d["status"]
            os.system("clear" if os.name != "nt" else "cls")
            el = d["now"] - st["started"]
            print(f"{run}  phase: {st['phase']}  elapsed {el:.0f} s  finished {len(st['finished'])}/{st['starts']}"
                  f"  process {'alive' if d['alive'] else 'gone' if d['alive'] is False else '?'}")
            print(f"{'start':>20} {'evals':>6} {'stage':>22} {'current':>11} {'best':>11} {'finished':>11} {'s':>7}")
            best_name, best_val = None, None
            for name, s in d["starts"].items():
                fin = st["finished"].get(name.replace("start_", "").lstrip("0") or "0")
                print(f"{name:>20} {s['evaluations']:6d} {s['label'][:22]:>22} {s['cost'] or 0:11.5f} "
                      f"{s['best'] or 0:11.5f} {fin if fin is not None else float('nan'):11.5f} {s['seconds']:7.0f}")
                if s["best"] is not None and (best_val is None or s["best"] < best_val):
                    best_name, best_val = name, s["best"]
            if best_name and d["starts"][best_name]["J"]:
                print(f"\nbest so far: {best_name} {best_val:.5f}")
                for key, v in d["starts"][best_name]["J"].items():
                    v0 = (d["starts"][best_name]["J0"] or st["start_couplings"]).get(key)
                    dv = f"  ({v[0] - v0[0]:+.3f} from start)" if v0 else ""
                    print(f"  {key:>14} {', '.join(f'{x:.3f}' for x in v)}{dv}")
            if st["phase"] == "finished":
                break
            time.sleep(interval)
    except KeyboardInterrupt:
        pass


def point_result(prob, z, args=None) -> dict:
    """A fit.json-like result at the parameter vector z of a run's problem: couplings (J_at_x, by the fit's keys),
    spectrum_parameters (field, decay-rate families, phase delay, ...), family_edges_hz, z_final, the objective
    as the fit ranks it (hard missing-peak rows) and, with the command's args, structure and exchange. Built the
    way fit_joint_series writes fit.json, so Studio can apply any point of a run (running or finished)."""
    joint = prob.joint
    z = np.asarray(z, float)
    out = {"couplings": {prob.key_of.get(n, n): {"J_at_x": joint.coupling_values(z, k).tolist()}
                         for k, n in enumerate(joint.coupling)}}
    shared = {n: float(z[joint.ntheta + i]) for i, n in enumerate(joint.shared)}
    out["spectrum_parameters"] = {prob.series[si]["id"]: {**shared, **{n: float(z[joint.nt + si * joint.nl + i])
                                                                        for i, n in enumerate(joint.local)}}
                                  for si in range(joint.ns)}
    out["family_edges_hz"] = [float(v) for v in joint.params[0].policy.family_edges_hz]
    out["x"] = [float(v) for v in joint.nodes]          # as fit.json: the concentration nodes of J_at_x
    out["z_final"] = [float(v) for v in z]
    smooth = getattr(joint, "peak_smooth", 0.0)
    try:
        if getattr(joint, "peaks", None) is not None:
            joint.peak_smooth = 0.0
        out["scores"] = [float(np.sum(joint.residual(z) ** 2))]
    finally:
        joint.peak_smooth = smooth
    if args is not None:
        out["structure"] = json.loads(args.structure)
        out["exchange"] = args.exchange
    return out


def read_point(run, start, n):
    """Parameter vector at evaluation n of one start: that evaluation's own vector ("x", recorded since
    2026-10-08) or, for older records, the best point recorded up to n. Returns {n, cost, best, z, kind}."""
    path = Path(run) / "monitor" / f"{start}.jsonl"
    own, best = None, None
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "n" not in r:
                continue
            if r["n"] > n:
                break
            if "z" in r:
                best = r
            if r["n"] == n and "cost" in r:
                own = r
    if own is not None and "x" in own:
        return {"n": n, "cost": own["cost"], "best": own["best"], "z": own["x"], "kind": "evaluation"}
    if best is not None:
        return {"n": best["n"], "cost": best.get("cost", best["best"]), "best": best["best"], "z": best["z"], "kind": "best up to"}
    raise ValueError(f"no parameter vector recorded up to evaluation {n} of {start}")


class SpectrumCache:
    """build_problem from a run's recorded command, once per run, for the on-demand spectrum view."""

    def __init__(self):
        self.lock = threading.Lock()
        self.problems = {}
        self.args = {}

    def get(self, run):
        run = str(run)
        with self.lock:
            if run not in self.problems:
                sys.path.insert(0, str(ROOT / "scripts"))
                import fit_joint_series as fj
                st = json.load(open(Path(run) / "monitor" / "status.json"))
                args = fj.make_parser().parse_args(st["argv"][1:])
                cwd = os.getcwd()
                try:
                    os.chdir(st.get("cwd", cwd))
                    self.problems[run] = fj.build_problem(args)
                    self.args[run] = args
                finally:
                    os.chdir(cwd)
            return self.problems[run]

    def spectrum(self, run, start, s=0, view=None, max_points=4000, z=None):
        """Data and model of spectrum s at the best point so far of one start (its full parameter vector), or at
        the vector z."""
        prob = self.get(run)
        joint = prob.joint
        if z is None:
            z = read_run(run)["starts"][start]["z"]
        if z is None:
            raise ValueError("no best point recorded yet")
        with self.lock:
            z = np.asarray(z, float)
            f = joint.forwards[s]
            pred = f.predict(joint.spectrum_vector(z, s))
            lo, hi = view if view else (float(f.f.min()), float(f.f.max()))
            sel = np.flatnonzero((f.f >= lo) & (f.f <= hi))
            sel = sel[::max(len(sel) // max_points, 1)]
            ranges = [list(map(float, r)) for r in prob.series[s].get("ranges", [(prob.lo, prob.hi)])]
            return {"f": f.f[sel].tolist(), "data_re": f.y[sel].real.tolist(), "data_im": f.y[sel].imag.tolist(),
                    "model_re": pred.model[sel].real.tolist(), "model_im": pred.model[sel].imag.tolist(),
                    "ranges": ranges, "extent": [float(f.f.min()), float(f.f.max())],
                    "objective": float(np.sum(joint.residual(z) ** 2))}


def make_handler(root, cache):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body, ctype="application/json"):
            data = body if isinstance(body, bytes) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            u = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(u.query).items()}
            try:
                if u.path in ("/", "/index.html"):
                    return self._send(200, UI.read_bytes(), "text/html; charset=utf-8")
                if u.path == "/favicon.ico":
                    return self._send(204, b"", "image/x-icon")
                if u.path == "/api/runs":
                    return self._send(200, {"runs": [str(p) for p in find_runs(root)]})
                if u.path == "/api/run":
                    d = read_run(q["run"])
                    log = Path(q["run"]) / "monitor" / "console.log"
                    d["console"] = open(log).read()[-6000:] if log.exists() else ""
                    return self._send(200, d)
                if u.path == "/api/spectrum":
                    view = [float(q["lo"]), float(q["hi"])] if "lo" in q else None
                    return self._send(200, cache.spectrum(q["run"], q["start"], view=view))
                self._send(404, {"error": "not found"})
            except Exception as exc:
                self._send(400, {"error": f"{type(exc).__name__}: {exc}"})
    return Handler


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path", nargs="?", default="runs", help="a run directory or a directory containing runs")
    ap.add_argument("--port", type=int, default=8770)
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--text", action="store_true", help="terminal table instead of the browser page")
    ap.add_argument("--interval", type=float, default=2.0)
    args = ap.parse_args()
    if args.text:
        runs = find_runs(args.path)
        if not runs:
            raise SystemExit(f"no run with a monitor below {args.path}")
        return text_view(runs[0], args.interval)
    server = ThreadingHTTPServer((args.host, args.port), make_handler(args.path, SpectrumCache()))
    print(f"fit monitor on http://{args.host}:{server.server_address[1]} (runs below {args.path}; Ctrl+C to quit)",
          flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
