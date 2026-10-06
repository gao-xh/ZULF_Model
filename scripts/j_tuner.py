"""Interactive J tuner: adjust couplings by hand (coarse and fine sliders), see the objective, the residual and a live
analysis while dragging, and start an automatic refinement from the current values.

    python scripts/j_tuner.py <the fit_joint_series options of the fit> --fit RUN/fit.json [--spectrum 0]
        [--port 8765] [--host 127.0.0.1] [--out RUN/tuner]

Then open http://127.0.0.1:8765 in a browser. The problem (data, processing, weights, missing-peak rows, rate
families, bounds) is built by fit_joint_series.build_problem from the same command line as the fit, so the numbers
shown are the fitter's own objective. --fit loads the couplings and spectrum parameters of a fit (keys missing from it
keep their structure values).

Shown on every change: the objective (with the missing-peak rows, as ranked by the fitter), the plain weighted
residual, the relative residual on the data cores, the cost inside the current view, their change against the
baseline (the loaded fit or the last accepted state), data / model / residual in the view, the model lines (the
analytic transitions) and the localized residual peaks (JointSeries.find_residual_peaks). On demand: the sources of a
residual peak (which couplings shift or split the lines under it, JointSeries.peak_sources) and the gradient of the
objective with a one-coupling Gauss-Newton step for every coupling.

Refinement: least squares on the full objective with only the chosen couplings (and optionally the decay rates of
the families in view, or every spectrum parameter) free, everything else held; progress (objective, couplings) is
streamed to the page and the fit can be stopped (the best point so far is kept). "Accept" makes the current state
the baseline; "Save" writes OUT/tuned.json, which fit_joint_series reads with --from-joint.
"""
import json
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import fit_joint_series as fj  # noqa: E402

UI = Path(__file__).with_name("j_tuner_ui.html")


class StopRefine(Exception):
    pass


def load_fit(prob, fit):
    """Couplings (J at the nearest node) and spectrum parameters (same spectrum id) of a fit_joint_series result
    into prob.z0; keys the fit does not have keep their start values. Returns the keys that were missing."""
    joint, key_of = prob.joint, prob.key_of
    rows = [int(np.argmin(np.abs(np.asarray(fit["x"]) - x))) for x in joint.nodes]
    table = np.array([[joint.coupling_values(prob.z0, k)[i] for k in range(joint.nc)] for i in range(joint.nn)])
    missing = []
    for k, n in enumerate(joint.coupling):
        c = fit["couplings"].get(key_of.get(n, n))
        if c is None:
            missing.append(key_of.get(n, n))
            continue
        table[:, k] = [c["J_at_x"][r] for r in rows]
    xs = [x.copy() for x in prob.xs0]
    from fit_joint_series import remap_family_rates
    edges_now = joint.params[0].policy.family_edges_hz
    for e, x in zip(prob.series, xs):
        spectrum = remap_family_rates(fit.get("spectrum_parameters", {}).get(e["id"], {}), fit.get("family_edges_hz"),
                                      edges_now)
        for n, v in spectrum.items():
            if n in joint.col:
                x[joint.col[n]] = v
    prob.z0 = np.clip(joint.pack(table, xs), prob.lower + 1e-9, prob.upper - 1e-9)
    return missing


class TuningSession:
    """State of one tuning session on spectrum s of a build_problem result: the current vector z, the baseline it is
    compared with, an undo stack, and at most one refinement running in a background thread."""

    def __init__(self, prob, s=0):
        self.prob, self.joint, self.s = prob, prob.joint, int(s)
        self.z = prob.z0.copy()
        self.initial = self.z.copy()
        self.baseline = self.z.copy()
        self.undo_stack = []
        self.lock = threading.Lock()
        self.job = None
        self.keys = [prob.key_of.get(n, n) for n in self.joint.coupling]
        self.leader = dict(zip(self.keys, self.joint.coupling))
        self._baseline_numbers = self.scores(self.baseline)

    # ---- couplings ----
    def _index(self, k):
        j = self.joint
        return k * j.m + (j.node_of[self.s] if j.shape == "free" else 0)

    def coupling_value(self, z, k):
        return float(self.joint.coupling_values(z, k)[self.joint.node_of[self.s]])

    def couplings(self):
        out = []
        for k, key in enumerate(self.keys):
            i = self._index(k)
            out.append({"key": key, "value": self.coupling_value(self.z, k),
                        "baseline": self.coupling_value(self.baseline, k),
                        "initial": self.coupling_value(self.initial, k),
                        "lower": float(self.prob.lower[i]), "upper": float(self.prob.upper[i])})
        return out

    def set_couplings(self, values, z=None):
        """values {key: Hz}: J of this spectrum's node (free shape) or the level of a monotone coupling moved so
        that its value at this node is the given one; clipped to the bounds."""
        z = self.z if z is None else z
        for key, v in values.items():
            k = self.keys.index(key)
            i = self._index(k)
            target = float(v) - self.coupling_value(z, k) + z[i]
            z[i] = float(np.clip(target, self.prob.lower[i] + 1e-9, self.prob.upper[i] - 1e-9))
        return z

    def push_undo(self):
        self.undo_stack.append(self.z.copy())
        del self.undo_stack[:-200]

    def undo(self):
        if self.undo_stack:
            self.z = self.undo_stack.pop()

    def reset(self, to="baseline"):
        self.push_undo()
        self.z = (self.initial if to == "initial" else self.baseline).copy()

    def accept(self):
        self.baseline = self.z.copy()
        self._baseline_numbers = self.scores(self.baseline)

    # ---- numbers ----
    def scores(self, z, view=None):
        """objective (fitter's ranking objective: all spectra, missing-peak rows, priors), plain (weighted residual
        of this spectrum), relative (residual on the data cores), window (plain cost inside view)."""
        j, s = self.joint, self.s
        f = j.forwards[s]
        objective = float(np.sum(j.residual(z) ** 2))
        pred = f.predict(j.spectrum_vector(z, s))
        m = f.data_signal_mask if f.data_signal_mask is not None else np.ones(len(f.y), bool)
        rel = float(np.linalg.norm(f.mismatch(pred.model, m)) / max(np.linalg.norm(f.mismatch(np.zeros_like(f.y), m)),
                                                                      1e-30))
        out = {"objective": objective, "plain": float(pred.residual @ pred.residual), "relative": rel}
        if view is not None:
            out["window"] = self._window_cost(pred.residual, view)
        return out, pred

    def _window_cost(self, r, view):
        f = self.joint.forwards[self.s]
        sel = np.flatnonzero((f.f >= view[0]) & (f.f <= view[1]))
        nf = len(f.f)
        idx = np.r_[sel, sel + nf] if len(r) == 2 * nf else sel
        return float(np.sum(r[idx] ** 2))

    def lines(self, z, view, min_relative=0.02):
        j, s = self.joint, self.s
        f, P = j.forwards[s], j.params[s]
        values = P.values(j.spectrum_vector(z, s))
        out = []
        for c, system in enumerate(P.systems(values)):
            tl = f.transitions(values, c, system)
            a = np.abs(np.asarray(tl.amplitudes))
            if not len(a):
                continue
            fr = np.asarray(tl.frequencies_hz)
            keep = (fr >= view[0]) & (fr <= view[1]) & (a >= min_relative * a.max())
            out += [{"component": c, "label": j.model.component_labels[c], "frequency_hz": float(v),
                     "relative_amplitude": float(w)} for v, w in zip(fr[keep], a[keep] / a.max())]
        return out

    def evaluate(self, view, analysis=True, max_points=4000):
        t0 = time.time()
        with self.lock:
            z = self.z.copy()
            numbers, pred = self.scores(z, view)
            base = self._baseline_numbers[0]
            f = self.joint.forwards[self.s]
            sel = np.flatnonzero((f.f >= view[0]) & (f.f <= view[1]))
            step = max(len(sel) // max_points, 1)
            sel = sel[::step]
            y, model = f.y[sel], pred.model[sel]
            key = (tuple(view), self.baseline.tobytes())
            if getattr(self, "_bw", (None,))[0] != key:            # baseline cost of this view, cached
                self._bw = (key, self._window_cost(f.predict(self.joint.spectrum_vector(self.baseline, self.s)).residual,
                                                   view))
            out = {"scores": numbers, "baseline": base, "baseline_window": self._bw[1],
                   "couplings": self.couplings(),
                   "spectrum": {"f": f.f[sel].tolist(), "data_re": y.real.tolist(), "data_im": y.imag.tolist(),
                                "model_re": model.real.tolist(), "model_im": model.imag.tolist(),
                                "weight": (f.weight[sel] / max(float(f.weight.max()), 1e-30)).tolist()},
                   "lines": self.lines(z, view)}
            if analysis:
                _, report = self.joint.find_residual_peaks(z)
                out["residual_peaks"] = [p for p in report if p["spectrum"] == self.s]
        out["seconds"] = round(time.time() - t0, 3)
        return out

    def sources(self, frequency_hz, assign_hz=0.6):
        with self.lock:
            src = self.joint.peak_sources(self.z, self.s, frequency_hz, assign_hz=assign_hz)
        key = {n: k for k, n in self.leader.items()}
        couplings = sorted(({"key": key.get(n, n), **v} for n, v in src["couplings"].items()),
                           key=lambda d: -d["local"])
        lines = [{**l, "df_dJ": {key.get(n, n): v for n, v in l["df_dJ"].items()}} for l in src["lines"]]
        return {"frequency_hz": frequency_hz, "couplings": couplings, "lines": lines}

    def gradient(self, view=None):
        """d objective / dJ for every coupling at this node and the one-coupling Gauss-Newton step
        -g_k / (2 sum_i J_ik^2) (what each coupling alone would move to lower the objective, to first order)."""
        j = self.joint
        with self.lock:
            z = self.z.copy()
            for f in j.forwards:
                f.jacobian_only = set(j.coupling)
            try:
                r = j.residual(z)
                jac = j.jacobian(z)
            finally:
                for f in j.forwards:
                    f.jacobian_only = None
        out = []
        for k, key in enumerate(self.keys):
            col = jac[:, self._index(k)]
            g = float(2 * col @ r)
            h = float(2 * col @ col)
            out.append({"key": key, "gradient": g, "step_hz": -g / h if h > 0 else 0.0})
        return out

    # ---- refinement ----
    def spectrum_names(self, which, view):
        """Spectrum parameters to free with the couplings: "none", "rates_in_view" (decay-rate families with a line
        in the view) or "all"."""
        j = self.joint
        if which == "all":
            return list(j.local) + list(j.shared)
        if which != "rates_in_view":
            return []
        edges = np.asarray(j.params[self.s].policy.family_edges_hz, float)
        fams = set()
        for line in self.lines(self.z, view, min_relative=0.0):
            fams.add((line["component"], int(np.searchsorted(edges, line["frequency_hz"], side="right"))))
        names = [f"c{c}.log_rate{fam}" for c, fam in sorted(fams)]
        names += [f"c{c}.log_rate" for c, _ in sorted(fams)]
        return [n for n in names if n in j.local or n in j.shared]

    def refine(self, keys, spectrum="rates_in_view", view=None, max_nfev=50):
        if self.job is not None and self.job["running"]:
            raise RuntimeError("a refinement is running")
        j = self.joint
        names = self.spectrum_names(spectrum, view) if view is not None or spectrum == "all" else []
        idx = [self._index(self.keys.index(k)) for k in keys]
        for n in names:
            idx.append(j.ntheta + j.shared.index(n) if n in j.shared else j.nt + self.s * j.nl + j.local.index(n))
        idx = np.asarray(sorted(set(idx)), int)
        if not len(idx):
            raise ValueError("nothing to refine: choose couplings or spectrum parameters")
        self.push_undo()
        job = {"running": True, "stop": False, "evaluations": 0, "history": [], "best": None, "best_z": None,
               "message": "", "free": [self.keys[k] for k in range(j.nc) if self._index(k) in idx] + names,
               "started": time.time(), "snapshot": None}
        self.job = job
        z_start = self.z.copy()
        lo, hi = self.prob.lower[idx], self.prob.upper[idx]
        wanted = {self.leader[k] for k in keys} | set(names)

        def full(p):
            zz = z_start.copy()
            zz[idx] = p
            return zz

        def res(p):
            if job["stop"]:
                raise StopRefine()
            zz = full(p)
            r = j.residual(zz)
            cost = float(r @ r)
            job["evaluations"] += 1
            job["history"].append(cost)
            if job["best"] is None or cost < job["best"]:
                job["best"], job["best_z"] = cost, zz
                now = time.time()
                if job["snapshot"] is None or now - job["snapshot"][0] > 1.0:
                    job["snapshot"] = (now, {"couplings": {self.keys[k]: self.coupling_value(zz, k)
                                                           for k in range(j.nc)}})
            return r

        def jac(p):
            if job["stop"]:
                raise StopRefine()
            return j.jacobian(full(p))[:, idx]

        def run():
            with self.lock:
                for f in j.forwards:
                    f.jacobian_only = wanted
                try:
                    p0 = np.clip(z_start[idx], lo + 1e-9, hi - 1e-9)
                    least_squares(res, p0, jac=jac, bounds=(lo, hi), x_scale="jac", max_nfev=max_nfev,
                                  ftol=1e-10, xtol=1e-10, gtol=1e-10)
                    job["message"] = "converged or reached max_nfev"
                except StopRefine:
                    job["message"] = "stopped"
                except Exception as exc:                        # report, keep the best point
                    job["message"] = f"error: {exc}"
                finally:
                    for f in j.forwards:
                        f.jacobian_only = None
                    if job["best_z"] is not None and job["best"] <= float(np.sum(j.residual(z_start) ** 2)):
                        self.z = job["best_z"].copy()
                    job["running"] = False
                    job["seconds"] = round(time.time() - job["started"], 1)

        threading.Thread(target=run, daemon=True).start()
        return {"free": job["free"]}

    def progress(self):
        job = self.job
        if job is None:
            return {"running": False}
        out = {k: job[k] for k in ("running", "evaluations", "best", "message", "free")}
        out["history"] = job["history"][-400:]
        out["seconds"] = job.get("seconds", round(time.time() - job["started"], 1))
        if job["snapshot"] is not None:
            out["couplings"] = job["snapshot"][1]["couplings"]
        return out

    def stop(self):
        if self.job is not None:
            self.job["stop"] = True

    # ---- output ----
    def result(self):
        j, z = self.joint, self.z
        numbers, _ = self.scores(z)
        return {"source": "scripts/j_tuner.py", "x": j.nodes.tolist(),
                "spectra": [{"id": e["id"], "x": e["x"]} for e in self.prob.series],
                "scores": [numbers["objective"]],
                "data_region_residuals": dict(zip([e["id"] for e in self.prob.series], j.data_residuals(z))),
                "couplings": {self.keys[k]: {"J_at_x": j.coupling_values(z, k).tolist()} for k in range(j.nc)},
                "shared": j.shared,
                "spectrum_parameters": {self.prob.series[si]["id"]: {
                    **{n: float(z[j.ntheta + i]) for i, n in enumerate(j.shared)},
                    **{n: float(z[j.nt + si * j.nl + i]) for i, n in enumerate(j.local)}} for si in range(j.ns)}}

    def save(self, path):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with self.lock:
            json.dump(self.result(), open(path, "w"), indent=1)
        return str(path)

    def state(self):
        f = self.joint.forwards[self.s]
        e = self.prob.series[self.s]
        ranges = [list(map(float, r)) for r in e.get("ranges", [(self.prob.lo, self.prob.hi)])]
        return {"spectrum": e["id"], "spectra": [x["id"] for x in self.prob.series], "ranges": ranges,
                "extent": [float(f.f.min()), float(f.f.max())], "complex": not f.real_only,
                "couplings": self.couplings(), "components": list(self.joint.model.component_labels),
                "spectrum_parameters": {n: float(self.z[self.joint.nt + self.s * self.joint.nl + i])
                                        for i, n in enumerate(self.joint.local)}}


def make_handler(session, out_dir):
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
            if self.path in ("/", "/index.html"):
                return self._send(200, UI.read_bytes(), "text/html; charset=utf-8")
            if self.path == "/api/state":
                return self._send(200, session.state())
            if self.path == "/api/progress":
                return self._send(200, session.progress())
            if self.path == "/favicon.ico":
                return self._send(204, b"", "image/x-icon")
            self._send(404, {"error": "not found"})

        def do_POST(self):
            try:
                n = int(self.headers.get("Content-Length", 0))
                req = json.loads(self.rfile.read(n) or b"{}")
                view = req.get("view")
                if self.path == "/api/evaluate":
                    if session.job is not None and session.job["running"]:
                        return self._send(409, {"error": "refinement running"})
                    if req.get("commit"):
                        session.push_undo()
                    if req.get("values"):
                        with session.lock:
                            session.set_couplings(req["values"])
                    return self._send(200, session.evaluate(view, analysis=req.get("analysis", True)))
                if self.path == "/api/sources":
                    return self._send(200, session.sources(float(req["frequency_hz"]), req.get("assign_hz", 0.6)))
                if self.path == "/api/gradient":
                    return self._send(200, {"couplings": session.gradient()})
                if self.path == "/api/refine":
                    return self._send(200, session.refine(req.get("couplings", []), req.get("spectrum", "rates_in_view"),
                                                          view, int(req.get("max_nfev", 50))))
                if self.path == "/api/stop":
                    session.stop()
                    return self._send(200, {"ok": True})
                if self.path in ("/api/accept", "/api/reset", "/api/undo"):
                    with session.lock:
                        if self.path == "/api/accept":
                            session.accept()
                        elif self.path == "/api/reset":
                            session.reset(req.get("to", "baseline"))
                        else:
                            session.undo()
                    return self._send(200, session.evaluate(view, analysis=req.get("analysis", True)))
                if self.path == "/api/save":
                    return self._send(200, {"path": session.save(Path(out_dir) / req.get("name", "tuned.json"))})
                self._send(404, {"error": "not found"})
            except Exception as exc:                               # the page shows the message
                self._send(400, {"error": f"{type(exc).__name__}: {exc}"})
    return Handler


def serve(session, host="127.0.0.1", port=8765, out_dir="runs/tuner"):
    server = ThreadingHTTPServer((host, port), make_handler(session, out_dir))
    return server


def main():
    ap = fj.make_parser()
    ap.add_argument("--fit", default="", help="fit.json (fit_joint_series or j_tuner) to start from")
    ap.add_argument("--spectrum", type=int, default=0, help="index of the spectrum of the series to tune")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    args = ap.parse_args()
    prob = fj.build_problem(args)
    if args.fit:
        missing = load_fit(prob, json.load(open(args.fit)))
        if missing:
            print("not in the fit (structure values kept):", ", ".join(missing))
    session = TuningSession(prob, args.spectrum)
    out_dir = args.out if args.out != "runs/processed/joint" else str(Path(args.fit).parent / "tuner" if args.fit
                                                                    else "runs/tuner")
    server = serve(session, args.host, args.port, out_dir)
    print(f"J tuner on http://{args.host}:{server.server_address[1]}  (objective {session._baseline_numbers[0]['objective']:.5f};"
          f" saves to {out_dir}/tuned.json; Ctrl+C to quit)", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
