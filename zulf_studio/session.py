"""Studio session: the simulator state and every operation on it, without a GUI.

The spin model comes from a structure specification (zulf_hypothesis.structure_spec): natural-abundance 13C
isotopologues with fixed abundance ratios, exchangeable protons decoupled ("fast") or kept ("slow"). Every
change rebuilds the model and recomputes the transitions exactly (zulf_core.physics.compute_transitions with
the static field of the session, D47), so the spectrum follows the sliders in real time for small molecules.

The displayed spectrum is a sum of complex Lorentzian lines (one decay rate for every line). It is a quick
look, not the fit model: a fit (fit_joint_series, started from the session) renders the model through the
same processing as the data, and its trace frames are shown as they were computed.

Units: couplings in Hz, field components in nT (B transverse to the detection axis, and B along it), decay
rate in 1/s, data phase in degrees and delay in ms.
"""
from __future__ import annotations

import csv
import json
import math
import re
import subprocess
import sys
import threading
import time
from collections import deque
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np

from zulf_core.physics.protocol import Protocol
from zulf_core.physics.transitions import compute_transitions
from zulf_hypothesis import exchange_variants
from zulf_hypothesis.builder import build_model
from zulf_hypothesis.motifs import MOTIFS
from zulf_hypothesis.structure_spec import fragment_from_spec, override_couplings, parse_coupling_key

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_SPEC = {"compound": "acetonitrile", "chain": {"groups": [["C1", "C", 3], ["C2", "C", 0]],
                                                       "bonds": [["C1", "C2"]]},
                "one_bond": {"C1": 136.3}}

# fit_joint_series options of the amine and acetonitrile analyses (docs/WORKFLOW.md)
FIT_DEFAULTS = {"starts": 8, "workers": 4, "max_nfev": 300, "trace": 40, "fit_field": True, "family_edges": "",
                "rate_bounds": "0.2,15", "signal_threshold": 2.5, "signal_taper": 4.0, "peak_penalty": 5.0,
                "peak_smooth": 0.03, "peak_min_sigma": 3.0, "component_search": "both", "extra_args": []}


def _resolve(path) -> Path:
    p = Path(path).expanduser()
    if p.exists() or p.is_absolute():
        return p
    return ROOT / p


class FitJob:
    """One fit_joint_series run in a subprocess; its output lines go to the session log."""

    def __init__(self, argv: List[str], out_dir: Path, log: Callable[[str, str], None],
                 done: Callable[["FitJob"], None]):
        self.argv, self.out_dir = argv, out_dir
        self.started = time.time()
        self.finished: Optional[float] = None
        self.returncode: Optional[int] = None
        self.best: Optional[float] = None
        self.starts_finished = 0
        self.last_line = ""
        self._log, self._done = log, done
        env = dict(**__import__("os").environ, OMP_NUM_THREADS="1")
        self.proc = subprocess.Popen(argv, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, env=env)
        self.thread = threading.Thread(target=self._read, daemon=True)
        self.thread.start()

    def _read(self):
        for line in self.proc.stdout:
            line = line.rstrip()
            if not line:
                continue
            self.last_line = line
            m = re.search(r"objective ([0-9.eE+-]+)", line)
            if m:
                v = float(m.group(1))
                self.best = v if self.best is None else min(self.best, v)
            if re.match(r"start \d+ finished", line):
                self.starts_finished += 1
            self._log(line, "fit")
        self.returncode = self.proc.wait()
        self.finished = time.time()
        self._log(f"fit finished with code {self.returncode} ({self.finished - self.started:.0f} s): {self.out_dir}",
                  "fit")
        self._done(self)

    @property
    def running(self) -> bool:
        return self.returncode is None

    def stop(self):
        if self.running:
            self.proc.terminate()

    def status(self) -> dict:
        return {"running": self.running, "returncode": self.returncode, "out_dir": str(self.out_dir),
                "seconds": round((self.finished or time.time()) - self.started, 1), "best_objective": self.best,
                "starts_finished": self.starts_finished, "last_line": self.last_line,
                "has_result": (self.out_dir / "fit.json").exists()}


class StudioSession:
    def __init__(self, spec: Optional[dict] = None, workspace="runs/studio", log_file: bool = True):
        self.lock = threading.RLock()
        self.workspace = _resolve(workspace)
        self.workspace.mkdir(parents=True, exist_ok=True)
        self.log_path = self.workspace / "studio.log" if log_file else None
        self.log_entries: deque = deque(maxlen=5000)
        self.listeners: List[Callable[[str], None]] = []
        self.log_listeners: List[Callable[[dict], None]] = []
        self.version = 0
        self.spec = dict(spec or DEFAULT_SPEC)
        self.overrides: Dict[str, float] = {}
        self.exchange = "fast"
        self.field_nt = [0.0, 0.0]            # (transverse, z)
        self.rate_per_s = 1.0
        self.view: Optional[List[float]] = None
        self.data: Optional[dict] = None
        self.data_phase_deg, self.data_delay_ms = 0.0, 0.0
        self.display = "re"
        self.fit_job: Optional[FitJob] = None
        self.trace: Optional[dict] = None
        self.trace_index = -1
        self._model_cache = None
        self._build()
        self.log(f"session started; structure {self.spec.get('compound', self.spec.get('motif', ''))}")

    # ---- events and log ----------------------------------------------------------------
    def log(self, message: str, source: str = "session"):
        entry = {"time": time.strftime("%H:%M:%S"), "source": source, "message": str(message)}
        self.log_entries.append(entry)
        if self.log_path is not None:
            with open(self.log_path, "a") as fh:
                fh.write(f"{time.strftime('%Y-%m-%d')} {entry['time']} [{source}] {message}\n")
        for cb in list(self.log_listeners):
            try:
                cb(entry)
            except Exception:
                pass

    def read_log(self, n: int = 200, source: Optional[str] = None) -> List[dict]:
        entries = [e for e in self.log_entries if source is None or e["source"] == source]
        return entries[-int(n):]

    def _changed(self, event: str):
        self.version += 1
        for cb in list(self.listeners):
            try:
                cb(event)
            except Exception as exc:
                self.log(f"listener error: {exc}", "session")

    # ---- model ---------------------------------------------------------------------------
    def _build(self):
        fragment = fragment_from_spec(self.spec)
        if self.overrides:
            fragment = override_couplings(fragment, self.overrides)
        fragment = exchange_variants(fragment, self.exchange)[0]
        model = build_model(fragment)
        self._model_cache = (fragment, model)
        if self.view is None:
            js = [abs(v["value"]) for v in self.couplings() if abs(v["value"]) >= 50] or [10.0]
            self.view = [max(0.0, 0.85 * min(js)), 2.15 * max(js)]

    def model(self):
        return self._model_cache

    def couplings(self) -> List[dict]:
        fragment, model = self._model_cache
        from zulf_hypothesis.fragment import pair
        out, seen = [], set()
        for key in list(model.coupling_names) + [k for k in self.overrides if k not in model.coupling_names]:
            if key in seen:
                continue
            seen.add(key)
            a, b = parse_coupling_key(key)
            value = float(fragment.couplings.get(pair(a, b), 0.0))
            out.append({"key": key, "value": value, "one_bond": abs(value) >= 50.0,
                        "unspecified": key in model.unspecified and key not in self.overrides,
                        "overridden": key in self.overrides})
        return out

    def components(self) -> List[dict]:
        _, model = self._model_cache
        return [{"label": lab, "contribution": float(c.contribution), "spins": len(c.system.isotopes)}
                for lab, c in zip(model.component_labels, model.interpretation.components)]

    def protocol(self) -> Protocol:
        bt, bz = self.field_nt
        return Protocol(field_ut=(1e-3 * bt, 0.0, 1e-3 * bz))

    # ---- setters (API tools) -------------------------------------------------------------
    def set_structure(self, spec: dict, keep_couplings: bool = False):
        with self.lock:
            old = (self.spec, self.overrides, self.view)
            self.spec = dict(spec)
            if not keep_couplings:
                self.overrides = {}
            self.view = None
            try:
                self._build()
            except Exception:
                self.spec, self.overrides, self.view = old
                self._build()
                raise
        self.log(f"structure: {json.dumps(spec)}")
        self._changed("structure")
        return self.state()

    def set_couplings(self, values: Dict[str, float]):
        with self.lock:
            old = dict(self.overrides)
            for key, v in values.items():
                parse_coupling_key(key)
                self.overrides[key] = float(v)
            try:
                self._build()
            except Exception:
                self.overrides = old
                self._build()
                raise
        self._changed("couplings")
        return self.couplings()

    def remove_coupling(self, key: str):
        with self.lock:
            self.overrides.pop(key, None)
            self._build()
        self._changed("couplings")
        return self.couplings()

    def set_exchange(self, mode: str):
        if mode not in ("fast", "slow"):
            raise ValueError("exchange must be 'fast' or 'slow'")
        with self.lock:
            self.exchange = mode
            self._build()
        self._changed("structure")
        return self.state()

    def set_field(self, transverse_nt: Optional[float] = None, z_nt: Optional[float] = None):
        with self.lock:
            if transverse_nt is not None:
                self.field_nt[0] = max(0.0, float(transverse_nt))
            if z_nt is not None:
                self.field_nt[1] = max(0.0, float(z_nt))
        self._changed("field")
        return {"transverse_nt": self.field_nt[0], "z_nt": self.field_nt[1],
                "magnitude_nt": math.hypot(*self.field_nt)}

    def set_linewidth(self, rate_per_s: float):
        with self.lock:
            self.rate_per_s = float(np.clip(rate_per_s, 1e-3, 1e3))
        self._changed("linewidth")
        return {"rate_per_s": self.rate_per_s, "fwhm_hz": self.rate_per_s / math.pi}

    def set_view(self, lo: float, hi: float):
        if not hi > lo:
            raise ValueError("view needs lo < hi")
        with self.lock:
            self.view = [float(lo), float(hi)]
        self._changed("view")
        return self.view

    def set_display(self, part: str = "re", phase_deg: Optional[float] = None, delay_ms: Optional[float] = None):
        if part not in ("re", "im", "abs"):
            raise ValueError("part must be re, im or abs")
        with self.lock:
            self.display = part
            if phase_deg is not None:
                self.data_phase_deg = float(phase_deg)
            if delay_ms is not None:
                self.data_delay_ms = float(delay_ms)
        self._changed("display")
        return {"part": self.display, "phase_deg": self.data_phase_deg, "delay_ms": self.data_delay_ms}

    # ---- simulation ----------------------------------------------------------------------
    def lines(self, min_relative: float = 0.0, view: Optional[List[float]] = None) -> List[dict]:
        """Every transition of every isotopologue: frequency, amplitude (times the abundance) and the amplitude
        relative to the strongest line."""
        with self.lock:
            _, model = self._model_cache
            prot = self.protocol()
            rows = []
            for lab, comp in zip(model.component_labels, model.interpretation.components):
                tl = compute_transitions(comp.system, prot)
                f = np.asarray(tl.frequencies_hz, float)
                a = np.real(np.asarray(tl.amplitudes)) * float(comp.contribution)
                rows += [{"component": lab, "frequency_hz": float(x), "amplitude": float(y)} for x, y in zip(f, a)]
        top = max((abs(r["amplitude"]) for r in rows), default=1.0) or 1.0
        lo, hi = view if view is not None else (-np.inf, np.inf)
        out = []
        for r in sorted(rows, key=lambda r: r["frequency_hz"]):
            r["relative"] = abs(r["amplitude"]) / top
            if r["relative"] >= min_relative and lo <= r["frequency_hz"] <= hi:
                out.append(r)
        return out

    def _data_display(self, f: np.ndarray, values: np.ndarray) -> np.ndarray:
        phase = np.radians(self.data_phase_deg) + 2 * np.pi * f * 1e-3 * self.data_delay_ms
        return values * np.exp(1j * phase)

    @staticmethod
    def _part(z: np.ndarray, part: str) -> np.ndarray:
        return z.real if part == "re" else z.imag if part == "im" else np.abs(z)

    def simulate(self, points: int = 3000, view: Optional[List[float]] = None) -> dict:
        """Simulated spectrum in the view (complex Lorentzians, one decay rate), the loaded data there (with the
        display phase) and the scale that matches the simulation to the data (least squares on the shown part)."""
        t0 = time.time()
        with self.lock:
            lo, hi = view or self.view
            lines = self.lines(view=(lo - 50, hi + 50))
            gamma = self.rate_per_s / (2 * np.pi)
            if self.data is not None:
                fd = self.data["freq"]
                sel = np.flatnonzero((fd >= lo) & (fd <= hi))
                step = max(len(sel) // points, 1)
                f = fd[sel[::step]]
                data = self._data_display(f, self.data["values"][sel[::step]])
            else:
                f = np.linspace(lo, hi, points)
                data = None
            sim = np.zeros(len(f), complex)
            for r in lines:
                if abs(r["frequency_hz"] - 0.5 * (lo + hi)) <= 0.5 * (hi - lo) + 40 * gamma + 1:
                    sim += r["amplitude"] * gamma / (gamma - 1j * (f - r["frequency_hz"]))
            part = self.display
            scale = 1.0
            if data is not None:
                m, d = self._part(sim, part), self._part(data, part)
                mm = float(m @ m)
                scale = float(m @ d) / mm if mm > 0 else 1.0
                if part == "abs" or scale <= 0:
                    scale = float(np.max(np.abs(d)) / max(np.max(np.abs(m)), 1e-30))
            out = {"f": f.tolist(), "sim_re": (scale * sim.real).tolist(), "sim_im": (scale * sim.imag).tolist(),
                   "part": part, "scale": scale, "view": [lo, hi],
                   "lines": [r for r in lines if lo <= r["frequency_hz"] <= hi and r["relative"] >= 1e-3]}
            if data is not None:
                out.update({"data_re": data.real.tolist(), "data_im": data.imag.tolist(),
                            "residual_rms": float(np.sqrt(np.mean((self._part(data, part) - scale * self._part(sim, part)) ** 2)))})
        out["seconds"] = round(time.time() - t0, 4)
        return out

    # ---- experimental spectrum -----------------------------------------------------------
    def load_spectrum(self, series: Optional[str] = None, index: int = 0, freq: Optional[str] = None,
                      values: Optional[str] = None, label: str = ""):
        """A processed spectrum: entry `index` of a series file (scripts/make_series_entry.py), or a frequency
        array and a (complex) value array (.npy)."""
        with self.lock:
            if series:
                spath = _resolve(series)
                entry = json.loads(spath.read_text())[int(index)]
                fr, va = np.load(_resolve(entry["freq"])), np.load(_resolve(entry["values"]))
                self.data = {"freq": np.asarray(fr, float), "values": np.asarray(va, complex), "label": entry["id"],
                             "series": str(spath.resolve()), "index": int(index),
                             "ranges": [list(map(float, r)) for r in entry.get("ranges", [])]}
            else:
                fr, va = np.load(_resolve(freq)), np.load(_resolve(values))
                self.data = {"freq": np.asarray(fr, float), "values": np.asarray(va, complex),
                             "label": label or Path(values).stem, "series": None, "index": 0, "ranges": []}
            if self.data["ranges"]:
                self.view = [min(r[0] for r in self.data["ranges"]), max(r[1] for r in self.data["ranges"])]
        self.log(f"spectrum loaded: {self.data['label']} ({len(self.data['freq'])} points)")
        self._changed("data")
        return {"label": self.data["label"], "points": len(self.data["freq"]), "ranges": self.data["ranges"]}

    # ---- fitting -------------------------------------------------------------------------
    def fit_command(self, **options) -> List[str]:
        """fit_joint_series command for the loaded series, starting from the current couplings and field."""
        if self.data is None or not self.data.get("series"):
            raise ValueError("load a series file (load_spectrum series=...) before fitting")
        o = dict(FIT_DEFAULTS, **options)
        ranges = self.data["ranges"] or [self.view]
        lo, hi = options.get("range") or [min(r[0] for r in ranges), max(r[1] for r in ranges)]
        out = Path(options.get("out") or self.workspace / "fits" /
                   f"{time.strftime('%Y%m%d-%H%M%S')}_{re.sub(r'[^A-Za-z0-9]+', '-', self.data['label'])}")
        current = {c["key"]: c["value"] for c in self.couplings()}
        argv = [sys.executable, str(ROOT / "scripts" / "fit_joint_series.py"), "--series", self.data["series"],
                "--real-only", "false", "--shape", "free", "--exchange", self.exchange, "--range", f"{lo},{hi}",
                "--structure", json.dumps(self.spec), "--couplings", json.dumps(current),
                "--rate-bounds", str(o["rate_bounds"]), "--signal-threshold", str(o["signal_threshold"]),
                "--signal-taper", str(o["signal_taper"]), "--peak-penalty", str(o["peak_penalty"]),
                "--peak-smooth", str(o["peak_smooth"]), "--peak-min-sigma", str(o["peak_min_sigma"]),
                "--starts", str(int(o["starts"])), "--workers", str(int(o["workers"])),
                "--max-nfev", str(int(o["max_nfev"])), "--trace", str(int(o["trace"])),
                "--component-search", str(o["component_search"]), "--out", str(out)]
        if o["fit_field"]:
            bt, bz = (v if v > 0 else 20.0 for v in self.field_nt)      # zero field has zero gradient (D47)
            argv += ["--fit-field", "--field-start", f"{1e-3 * bt},{1e-3 * bz}"]
        if o["family_edges"]:
            argv += ["--family-edges", str(o["family_edges"])]
        return argv + [str(a) for a in o["extra_args"]]

    def start_fit(self, **options) -> dict:
        if self.fit_job is not None and self.fit_job.running:
            raise RuntimeError("a fit is already running")
        argv = self.fit_command(**options)
        out = Path(argv[argv.index("--out") + 1])
        self.log("fit: " + " ".join(argv[1:]), "fit")
        self.fit_job = FitJob(argv, out, self.log, lambda job: self._changed("fit_done"))
        self._changed("fit_started")
        return self.fit_job.status()

    def fit_status(self) -> dict:
        return self.fit_job.status() if self.fit_job is not None else {"running": False, "out_dir": None}

    def stop_fit(self) -> dict:
        if self.fit_job is not None:
            self.fit_job.stop()
            self.log("fit stop requested", "fit")
        return self.fit_status()

    def apply_fit(self, run_dir: Optional[str] = None, spectrum: Optional[str] = None) -> dict:
        """Couplings, field and decay rate of a finished fit (fit.json) into the session; loads its trace."""
        run = _resolve(run_dir) if run_dir else (self.fit_job.out_dir if self.fit_job else None)
        if run is None:
            raise ValueError("no run directory")
        fit = json.loads((run / "fit.json").read_text())
        couplings = {k: float(v["J_at_x"][0]) for k, v in fit["couplings"].items()}
        sp_all = fit.get("spectrum_parameters", {})
        sp = sp_all.get(spectrum) if spectrum else next(iter(sp_all.values()), {})
        with self.lock:
            self.overrides.update(couplings)
            self._build()
            bt = sp.get("field_transverse_ut", sp.get("field_perp_ut"))
            if bt is not None or "field_z_ut" in sp:
                self.field_nt = [1e3 * float(bt or 0.0), 1e3 * float(sp.get("field_z_ut", 0.0))]
            rates = [math.exp(v) for k, v in sp.items() if ".log_rate" in k]
            if rates:
                self.rate_per_s = float(np.median(rates))
        self.log(f"applied fit {run}: score {fit.get('scores', [None])[0]}, field {self.field_nt} nT, "
                 f"rate {self.rate_per_s:.3g} 1/s (median of {len(rates)} families)")
        if (run / "trace.json").exists():
            self.load_trace(str(run))
        self._changed("couplings")
        return {"couplings": couplings, "field_nt": self.field_nt, "rate_per_s": self.rate_per_s}

    # ---- fit trace (progress slider) -----------------------------------------------------
    def load_trace(self, run_dir: str) -> dict:
        run = _resolve(run_dir)
        meta = json.loads((run / "trace.json").read_text())
        arrays = np.load(run / "trace.npz")
        with self.lock:
            self.trace = {"run": str(run), "meta": meta, "f": np.asarray(arrays["f0"]), "y": np.asarray(arrays["y0"]),
                          "models": np.asarray(arrays["model0"])}
            self.trace_index = len(meta["frames"]) - 1
        self.log(f"trace loaded: {run} ({len(meta['frames'])} frames of {meta['evaluations']} evaluations)")
        self._changed("trace")
        return {"frames": len(meta["frames"]), "evaluations": meta["evaluations"]}

    def trace_frame(self, index: int, apply: bool = False) -> dict:
        if self.trace is None:
            raise ValueError("no trace loaded")
        frames = self.trace["meta"]["frames"]
        index = int(np.clip(index, 0, len(frames) - 1))
        fr = frames[index]
        with self.lock:
            self.trace_index = index
        if apply:
            self.set_couplings({k: float(v[0]) for k, v in fr["J"].items()})
        else:
            self._changed("trace_frame")
        return {"index": index, "frames": len(frames), "evaluation": fr["evaluation"], "stage": fr["stage"],
                "objective": fr["objective"], "couplings": {k: float(v[0]) for k, v in fr["J"].items()}}

    def clear_trace(self):
        with self.lock:
            self.trace, self.trace_index = None, -1
        self._changed("trace")
        return {"ok": True}

    # ---- state and export ----------------------------------------------------------------
    def state(self) -> dict:
        with self.lock:
            return {"version": self.version, "structure": self.spec, "exchange": self.exchange,
                    "components": self.components(), "couplings": self.couplings(),
                    "field_nt": {"transverse": self.field_nt[0], "z": self.field_nt[1],
                                 "magnitude": math.hypot(*self.field_nt)},
                    "rate_per_s": self.rate_per_s, "view": self.view,
                    "display": {"part": self.display, "phase_deg": self.data_phase_deg,
                                "delay_ms": self.data_delay_ms},
                    "data": None if self.data is None else {k: self.data[k] for k in ("label", "series", "index",
                                                                                       "ranges")},
                    "fit": self.fit_status(),
                    "trace": None if self.trace is None else {"run": self.trace["run"], "index": self.trace_index,
                                                              "frames": len(self.trace["meta"]["frames"])}}

    def list_motifs(self) -> List[dict]:
        return [{"name": n, "description": m.description,
                 "one_bond": {s.site: round(0.5 * sum(s.j_range_hz), 1) for s in m.one_bond}}
                for n, m in MOTIFS.items()]

    def export(self, directory: Optional[str] = None) -> dict:
        """parameters.json (structure, couplings, field, rate), lines.csv and spectrum.csv into a directory."""
        out = _resolve(directory) if directory else self.workspace / "exports" / time.strftime("%Y%m%d-%H%M%S")
        out.mkdir(parents=True, exist_ok=True)
        st = self.state()
        (out / "parameters.json").write_text(json.dumps(
            {k: st[k] for k in ("structure", "exchange", "components", "couplings", "field_nt", "rate_per_s",
                                "view", "data")}, indent=1))
        with open(out / "lines.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["component", "frequency_hz", "amplitude", "relative"])
            for r in self.lines(min_relative=1e-4):
                w.writerow([r["component"], f"{r['frequency_hz']:.6f}", f"{r['amplitude']:.6g}", f"{r['relative']:.6g}"])
        sim = self.simulate()
        with open(out / "spectrum.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            cols = ["frequency_hz", "sim_re", "sim_im"] + (["data_re", "data_im"] if "data_re" in sim else [])
            w.writerow(cols)
            for row in zip(*(sim[c] if c != "frequency_hz" else sim["f"] for c in cols)):
                w.writerow([f"{v:.8g}" for v in row])
        self.log(f"exported to {out}")
        return {"directory": str(out), "files": ["parameters.json", "lines.csv", "spectrum.csv"]}
