"""Studio session: the simulator state and every operation on it, without a GUI.

The spin model comes from a structure specification (zulf_hypothesis.structure_spec): natural-abundance 13C
isotopologues with fixed abundance ratios, exchangeable protons decoupled ("fast") or kept ("slow"). Every
change rebuilds the model and recomputes the transitions exactly (zulf_core.physics.compute_transitions with
the static field of the session, D47), so the spectrum follows the sliders in real time for small molecules.

The displayed spectrum is a sum of complex Lorentzian lines a gamma / (gamma + i (f - f_k)), gamma = rate / 2 pi
(the sign convention of the processed spectra, numpy FFT exp(-2 pi i f t); one decay rate). It is a quick look, not the fit model: a fit (fit_joint_series, started from the session) renders the model through the
same processing as the data, and its trace frames are shown as they were computed.

Units: couplings in Hz, field components in nT (B transverse to the detection axis, and B along it), decay
rate in 1/s, data phase in degrees and delay in ms.
"""
from __future__ import annotations

import csv
import json
import math
import os
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
                "peak_smooth": 0.03, "peak_min_sigma": 3.0, "component_search": "both", "precision": 0.01,
                "extra_args": []}


def _resolve(path) -> Path:
    p = Path(path).expanduser()
    if p.exists() or p.is_absolute():
        return p
    return ROOT / p


class FitJob:
    """One script run in a subprocess (a fit, or a figure); its output lines go to the session log under
    `source`."""

    def __init__(self, argv: List[str], out_dir: Path, log: Callable[[str, str], None],
                 done: Callable[["FitJob"], None], source: str = "fit", kind: str = "", title: str = ""):
        self.argv, self.out_dir, self.source = argv, out_dir, source
        self.started = time.time()
        self.finished: Optional[float] = None
        self.returncode: Optional[int] = None
        self.best: Optional[float] = None
        self.starts_finished = 0
        self.last_line = ""
        self._log, self._done = log, done
        env = {**os.environ, **ONE_THREAD}          # one BLAS thread per worker process (inherited values replaced)
        self.kind, self.title = kind or source, title or Path(out_dir).name
        self.total = int(argv[argv.index("--starts") + 1]) if "--starts" in argv else None
        self.workers = int(argv[argv.index("--workers") + 1]) if "--workers" in argv else 1
        self.phase = "starting"
        self.proc = subprocess.Popen(argv, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, bufsize=1, env=env)
        self.pid = self.proc.pid
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
            self.phase = _phase_of(line, self.phase)
            self._log(line, self.source)
        self.returncode = self.proc.wait()
        self.proc.stdout.close()
        self.finished = time.time()
        self._log(f"{self.source} finished with code {self.returncode} ({self.finished - self.started:.0f} s): "
                  f"{self.out_dir}", self.source)
        self._done(self)

    @property
    def running(self) -> bool:
        return self.returncode is None

    def stop(self):
        if self.running:
            self.proc.terminate()

    def status(self) -> dict:
        return {"kind": self.kind, "title": self.title, "pid": self.pid, "workers": self.workers,
                "running": self.running, "returncode": self.returncode, "out_dir": str(self.out_dir),
                "seconds": round((self.finished or time.time()) - self.started, 1), "best_objective": self.best,
                "starts_finished": self.starts_finished, "starts_total": self.total, "phase": self.phase,
                "last_line": self.last_line, "has_result": (self.out_dir / "fit.json").exists()}


SESSION_FORMAT = "zulf-studio-session"
SESSION_SUFFIX = ".zulfstudio"
ONE_THREAD = {"OMP_NUM_THREADS": "1", "VECLIB_MAXIMUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
              "MKL_NUM_THREADS": "1"}
# scripts whose processes count as analysis load (machine_status)
ANALYSIS_SCRIPTS = {"fit_joint_series.py": 4, "analyze_sample.py": 4, "paper_figure.py": 1, "average_scans.py": 1,
                    "make_series_entry.py": 1, "reliability_series.py": 4, "regression_confirmed.py": 4}


def _phase_of(line: str, previous: str) -> str:
    """Short stage of a running fit or analysis from one output line (for the progress display)."""
    rules = (("component search (start)", "component search before the starts"),
             ("component search (end)", "component search after the starts"), ("global refit", "global refit"),
             ("family-edges auto", "choosing rate families"), ("guard rows", "setting up"),
             ("start ", "multi-start fit"), ("crop ", "processing"), ("figure", "figure"),
             ("averag", "averaging scans"), ('"score"', "writing results"))
    for key, phase in rules:
        if key in line:
            return phase
    return previous


def analysis_processes(ps_text: Optional[str] = None) -> List[dict]:
    """Analysis scripts running on this machine (any user session, any program): pid, script, workers, output.
    Read from `ps -Ao pid=,ppid=,command=` (or `ps_text`) by parsing each command into its arguments; a process
    whose parent is a counted run is one of its forked workers and is counted there (only counted, never
    signalled)."""
    text = ps_text
    if text is None:
        try:
            text = subprocess.run(["ps", "-Ao", "pid=,ppid=,command="], capture_output=True, text=True,
                                  timeout=5).stdout
        except (OSError, subprocess.SubprocessError):
            return []
    rows = []
    for row in text.splitlines():
        parts = row.split()
        if len(parts) < 4 or "python" not in Path(parts[2]).name:
            continue
        script = next((Path(a).name for a in parts[3:5] if Path(a).name in ANALYSIS_SCRIPTS), None)
        if script is not None:
            rows.append((int(parts[0]), int(parts[1]), script, parts[3:]))
    pids = {r[0] for r in rows}
    out = []
    for pid, ppid, script, args in rows:
        if ppid in pids:                  # a worker forked by a counted run (same command line): counted there
            continue
        workers = int(args[args.index("--workers") + 1]) if "--workers" in args[:-1] and \
            args[args.index("--workers") + 1].isdigit() else ANALYSIS_SCRIPTS[script]
        if "--starts" in args[:-1] and args[args.index("--starts") + 1].isdigit():
            workers = min(workers, int(args[args.index("--starts") + 1]))
        outdir = args[args.index("--out") + 1] if "--out" in args[:-1] else ""
        out.append({"pid": pid, "script": script, "workers": workers, "out": outdir})
    return out


def fit_structure(run: Path, fit: dict):
    """The structure specification of a fit: fit.json "structure" (fits since 2026-10-08), else the --structure
    of the command in its RUN_LOG.md; None when neither has one (a list of structures is returned as written)."""
    if fit.get("structure"):
        return fit["structure"]
    log = run / "RUN_LOG.md"
    if log.exists():
        import shlex
        for line in log.read_text().splitlines():
            if "fit_joint_series.py" in line and "--structure" in line:
                try:
                    args = shlex.split(line.strip().strip("`"))
                    return json.loads(args[args.index("--structure") + 1])
                except (ValueError, IndexError):
                    return None
    return None


_MONITOR_CACHE: Dict[str, tuple] = {}


def monitor_progress(run: Path) -> dict:
    """Progress of a fit_joint_series run from RUN/monitor (written while every start runs): elapsed seconds,
    best objective so far, starts finished and total, phase with the number of evaluations. Empty without a
    monitor; cached for 2 s per run."""
    key = str(run)
    hit = _MONITOR_CACHE.get(key)
    if hit and time.time() - hit[0] < 2.0:
        return hit[1]
    out = {}
    if (run / "monitor" / "status.json").exists():
        try:
            sys.path.insert(0, str(ROOT / "scripts"))
            from fit_monitor import read_run
            r = read_run(run, max_points=2)
            starts = list(r["starts"].values())
            bests = [v["best"] for v in starts if v["best"] is not None]
            evals = sum(v["evaluations"] for v in starts)
            out = {"seconds": max([v["seconds"] for v in starts] or [0.0]) or None,
                   "best_objective": min(bests) if bests else None,
                   "starts_finished": sum(not v["running"] for v in starts) if starts else None,
                   "starts_total": r["status"].get("starts") or None,
                   "phase": f"{r['status'].get('phase', 'running')}, {evals} evaluations" if starts else None}
        except Exception:                               # a monitor being written: no progress this time
            out = {}
    _MONITOR_CACHE[key] = (time.time(), out)
    return out


def cpu_cores() -> dict:
    """Logical cores, and on Apple silicon the performance and efficiency cores."""
    info = {"logical": os.cpu_count() or 1, "performance": None, "efficiency": None}
    if sys.platform == "darwin":
        for key, name in (("performance", "hw.perflevel0.physicalcpu"), ("efficiency", "hw.perflevel1.physicalcpu")):
            try:
                info[key] = int(subprocess.run(["sysctl", "-n", name], capture_output=True, text=True,
                                               timeout=2).stdout.strip())
            except (OSError, ValueError, subprocess.SubprocessError):
                pass
    return info


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
        self.scale_lock: Optional[float] = None      # fixed display scale of the simulation (None: automatic)
        self.fit_job: Optional[FitJob] = None
        self.figure_job: Optional[FitJob] = None
        self.blind_job: Optional[FitJob] = None
        self.jobs: List[FitJob] = []                 # every background job of the session (fits, blind, figures)
        self.session_file: Optional[str] = None      # last saved or opened session file
        self.mode = "simulate"                       # task mode of the window (D58): simulate, process, fit, blind
        self.process: Optional[dict] = None          # Process mode: source FID, its scans record, recipe, preview
        self.figure: Optional[dict] = None           # last figure: directory, files, manual flag
        self.applied: Optional[dict] = None          # run directory of the applied fit and the state it gave
        self.fit_curve: Optional[dict] = None        # exact model of the applied fit (fit_model_curve)
        self._fit_problems = None                    # fit_monitor.SpectrumCache: the problems of fit runs
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
        if "spin_system" in self.spec:                # PLAN 8c: isotopes and a J matrix typed in
            from .spin_model import CustomModel
            self._model_cache = (None, CustomModel(self.spec))
            if self.view is None:
                js = [abs(v["value"]) for v in self.couplings() if abs(v["value"]) >= 50] or [10.0]
                self.view = [max(0.0, 0.85 * min(js)), 2.15 * max(js)]
            return
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
        if fragment is None:
            return model.couplings()
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
    def structure_from_molecule(self, molecule: str, compound: Optional[str] = None,
                                one_bond: Optional[Dict[str, float]] = None) -> dict:
        """Build the model from a molecule (SMILES or mol-file text; zulf_hypothesis.molecule): every heavy atom a
        site, protons on O and S left out (exchangeable), the symmetry found from the molecule, 1J starting guesses
        from the hybridisation (one_bond overrides them by site label C1, C2, ...). The interface for drawing
        molecules later. Returns the specification and its notes (warnings)."""
        from zulf_hypothesis.molecule import structure_from_molecule
        spec = structure_from_molecule(molecule, one_bond=one_bond, compound=compound)
        notes = spec.pop("notes", [])
        self.set_structure(spec)
        for note in notes:
            self.log(f"molecule: {note}", "session")
        return {"structure": spec, "notes": notes, "components": [c["label"] for c in self.components()]}

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
        if "spin_system" in self.spec:
            from .spin_model import set_values
            with self.lock:
                old = self.spec
                self.spec = set_values(self.spec, values)
                try:
                    self._build()
                except Exception:
                    self.spec = old
                    self._build()
                    raise
            self._changed("couplings")
            return self.couplings()
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

    def auto_phase(self, method: str = "model", fit_delay: bool = True, delay_range_ms=None,
                   view: Optional[List[float]] = None) -> dict:
        """Display phase (degrees) and delay (ms) of the loaded data, on top of its own phasing.

        method "model": the phase0 and delay that best match the current simulation in the view (data x
        exp(i(phase0 + 2 pi f delay)) closest to a positive multiple of the simulation: delay on a 0.005 ms grid,
        phase0 in closed form); needs the lines roughly in place. method "data": model-free,
        zulf_core.render.phasing.estimate_phase on the data in the view (phase of the strongest peaks against
        frequency, modulo pi; the overall sign follows the peak weight). fit_delay=False keeps the delay. Default
        delay range: +-3 ms (model), +-0.5 ms (data: with lines at J and 2J, delays about 1 / (2 J) apart are
        nearly equivalent for a model-free criterion, so a wide range lands on an alias)."""
        if self.data is None:
            raise ValueError("no data loaded")
        if method not in ("model", "data"):
            raise ValueError("method must be 'model' or 'data'")
        if delay_range_ms is None:
            delay_range_ms = (-3.0, 3.0) if method == "model" else (-0.5, 0.5)
        with self.lock:
            lo, hi = view or self.view
            fd = self.data["freq"]
            sel = (fd >= lo) & (fd <= hi)
            ranges = self.data.get("ranges") or []
            if ranges:
                sel &= np.any([(fd >= r[0]) & (fd <= r[1]) for r in ranges], axis=0)
            f, d = fd[sel], self.data["values"][sel]
            if method == "data":
                from zulf_core.render.phasing import estimate_phase
                rng = (-1e-3 * delay_range_ms[1], -1e-3 * delay_range_ms[0]) if fit_delay else \
                    (-1e-3 * self.data_delay_ms, -1e-3 * self.data_delay_ms)
                est = estimate_phase(d, f, delay_range_s=rng, delay_step_s=5e-6 if fit_delay else 1.0)
                phase, delay = -np.degrees(est["phase0_rad"]), -1e3 * est["delay_s"]
                quality = 1.0 - est["misfit"]
            else:
                gamma = self.rate_per_s / (2 * np.pi)
                sim = np.zeros(len(f), complex)
                for r in self.lines(view=(lo - 20, hi + 20)):
                    sim += r["amplitude"] * gamma / (gamma + 1j * (f - r["frequency_hz"]))
                w = np.conj(sim) * d
                taus = (np.arange(delay_range_ms[0], delay_range_ms[1] + 1e-9, 0.005) if fit_delay
                        else np.array([self.data_delay_ms]))
                best = None
                for chunk in np.array_split(taus, max(1, len(taus) // 200)):
                    a = np.exp(2j * np.pi * np.outer(1e-3 * chunk, f)) @ w
                    k = int(np.argmax(np.abs(a)))
                    if best is None or abs(a[k]) > abs(best[1]):
                        best = (float(chunk[k]), a[k])
                delay, amp = best
                phase = -np.degrees(np.angle(amp))
                quality = float(abs(amp) / max(np.linalg.norm(sim) * np.linalg.norm(d), 1e-300))
            phase = float((phase + 180.0) % 360.0 - 180.0)
            self.data_phase_deg, self.data_delay_ms = phase, float(delay)
        self.log(f"auto phase ({method}{'' if fit_delay else ', delay held'}): phase {phase:.1f} deg, "
                 f"delay {delay:.3f} ms, match {quality:.3f}")
        self._changed("display")
        return {"phase_deg": phase, "delay_ms": float(delay), "method": method, "match": quality}

    def lock_scale(self, lock: bool = True, value: Optional[float] = None):
        """Freeze the display scale of the simulation at `value` (default: the current automatic scale), or
        release it (lock=False)."""
        with self.lock:
            if not lock:
                self.scale_lock = None
            else:
                if value is None:
                    self.scale_lock = None
                    value = self.simulate(points=1500)["scale"]
                self.scale_lock = float(value)
        self._changed("display")
        return {"locked": self.scale_lock is not None, "scale": self.scale_lock}

    # ---- simulation ----------------------------------------------------------------------
    def lines(self, min_relative: float = 0.0, view: Optional[List[float]] = None) -> List[dict]:
        """Every transition of every isotopologue: frequency, amplitude (times the abundance), the amplitude
        relative to the strongest line, and its decay rate in 1/s (decay_per_s): with an applied, unchanged fit the
        fitted rate of the line's family (the family is found from the fit's family edges as the forward model
        does), otherwise the one decay rate of the quick look."""
        with self.lock:
            _, model = self._model_cache
            prot = self.protocol()
            rows = []
            fitted = self._fitted_rates()
            for ci, (lab, comp) in enumerate(zip(model.component_labels, model.interpretation.components)):
                tl = compute_transitions(comp.system, prot)
                f = np.asarray(tl.frequencies_hz, float)
                a = np.real(np.asarray(tl.amplitudes)) * float(comp.contribution)
                rates = np.full(len(f), float(self.rate_per_s))
                if fitted is not None and ci in fitted["rates"]:
                    fam = np.searchsorted(fitted["edges"], f, side="right") if len(fitted["edges"]) else \
                        np.zeros(len(f), int)
                    fam = np.minimum(fam, len(fitted["rates"][ci]) - 1)
                    rates = fitted["rates"][ci][fam]
                rows += [{"component": lab, "frequency_hz": float(x), "amplitude": float(y), "decay_per_s": float(r)}
                         for x, y, r in zip(f, a, rates)]
        top = max((abs(r["amplitude"]) for r in rows), default=1.0) or 1.0
        lo, hi = view if view is not None else (-np.inf, np.inf)
        out = []
        for r in sorted(rows, key=lambda r: r["frequency_hz"]):
            r["relative"] = abs(r["amplitude"]) / top
            if r["relative"] >= min_relative and lo <= r["frequency_hz"] <= hi:
                out.append(r)
        return out

    def _fitted_rates(self) -> Optional[dict]:
        """Decay rates of the applied fit per component index (family order) and its family edges in Hz; None
        without an applied fit or once the parameters were changed."""
        if self.applied is None or self.applied["state"] != self._snapshot():
            return None
        fit = self.applied["fit"]
        sp_all = fit.get("spectrum_parameters", {})
        sp = sp_all.get(self.applied.get("spectrum")) or next(iter(sp_all.values()), {})
        rates = {}
        for key, v in sp.items():
            m = re.fullmatch(r"c(\d+)\.log_rate(\d+)", key)
            if m:
                rates.setdefault(int(m.group(1)), {})[int(m.group(2))] = math.exp(v)
        if not rates:
            return None
        return {"edges": np.asarray(fit.get("family_edges_hz") or [], float),
                "rates": {c: np.array([fams[k] for k in sorted(fams)]) for c, fams in rates.items()}}

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
                    sim += r["amplitude"] * gamma / (gamma + 1j * (f - r["frequency_hz"]))
            part = self.display
            scale = 1.0
            if self.scale_lock is not None:
                scale = self.scale_lock
            elif data is not None:
                # least squares on the magnitudes: positive and continuous in the parameters, the same for every
                # shown part (a fit on the real part fades to zero and changes sign when a line moves through a
                # dispersive data feature)
                m, d = np.abs(sim), np.abs(data)
                mm = float(m @ m)
                scale = float(m @ d) / mm if mm > 0 else 1.0
            out = {"f": f.tolist(), "sim_re": (scale * sim.real).tolist(), "sim_im": (scale * sim.imag).tolist(),
                   "part": part, "scale": scale, "view": [lo, hi],
                   "lines": [r for r in lines if lo <= r["frequency_hz"] <= hi and r["relative"] >= 1e-3]}
            curve = self._fit_curve_shown()
            if curve is not None:                     # the applied fit's own model, same display phase as the data
                fc = curve["f"]
                inside = (f >= fc[0]) & (f <= fc[-1])
                fm = np.full(len(f), np.nan, complex)
                fm[inside] = (np.interp(f[inside], fc, curve["model"].real)
                              + 1j * np.interp(f[inside], fc, curve["model"].imag))
                fm = self._data_display(f, fm)
                out.update({"fit_re": fm.real.tolist(), "fit_im": fm.imag.tolist(), "fit_source": curve["source"],
                            "fit_objective": curve["objective"]})
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
                             "source_fid": entry.get("source_fid"),
                             "ranges": [list(map(float, r)) for r in entry.get("ranges", [])]}
            else:
                fr, va = np.load(_resolve(freq)), np.load(_resolve(values))
                self.data = {"freq": np.asarray(fr, float), "values": np.asarray(va, complex),
                             "label": label or Path(values).stem, "series": None, "index": 0, "ranges": [],
                             "source_fid": None}
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
                "--component-search", str(o["component_search"]), "--precision", str(float(o["precision"])),
                "--out", str(out)]
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
        self.fit_job = FitJob(argv, out, self.log, lambda job: self._changed("fit_done"), kind="fit",
                              title=f"fit {self.data['label']}")
        self.jobs.append(self.fit_job)
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
        structure = fit_structure(run, fit)
        bare = {k: v for k, v in self.spec.items() if k != "molecule"}
        if isinstance(structure, dict) and {k: v for k, v in structure.items() if k != "molecule"} != bare:
            if ("molecule" in self.spec and "molecule" not in structure
                    and structure.get("compound") == self.spec.get("compound")):
                structure = dict(structure, molecule=self.spec["molecule"])     # same compound: keep the drawing
            self.set_structure(structure)
            if fit.get("exchange") in ("fast", "slow"):
                self.set_exchange(fit["exchange"])
        sp_all = fit.get("spectrum_parameters", {})
        sp = sp_all.get(spectrum) if spectrum else next(iter(sp_all.values()), {})
        with self.lock:
            if "spin_system" in self.spec:            # typed-in spin system: the values live in the J matrices
                from .spin_model import set_values
                self.spec = set_values(self.spec, couplings)
            else:
                self.overrides.update(couplings)
            self._build()
            bt = sp.get("field_transverse_ut", sp.get("field_perp_ut"))
            if bt is not None or "field_z_ut" in sp:
                self.field_nt = [1e3 * float(bt or 0.0), 1e3 * float(sp.get("field_z_ut", 0.0))]
            rates = [math.exp(v) for k, v in sp.items() if ".log_rate" in k]
            if rates:
                self.rate_per_s = float(np.median(rates))
        self.applied = {"run": str(run), "fit": fit, "spectrum": spectrum or next(iter(sp_all), None),
                        "state": self._snapshot()}
        self.log(f"applied fit {run}: score {fit.get('scores', [None])[0]}, field {self.field_nt} nT, "
                 f"rate {self.rate_per_s:.3g} 1/s (median of {len(rates)} families)")
        if (run / "trace.json").exists():
            self.load_trace(str(run))
        self._changed("couplings")
        return {"couplings": couplings, "field_nt": self.field_nt, "rate_per_s": self.rate_per_s}

    def fit_model_curve(self) -> Optional[dict]:
        """Exact model of the applied fit for the loaded spectrum: the fit's own forward model (phase, delay, decay
        rates of every family, amplitudes) at its final parameter vector (fit.json z_final; older fits: the best
        point of the monitor record), built from the run's recorded command (scripts/fit_monitor.py). It is the
        curve the fit compared with the data, on the data's grid and in the data's units. Slow (the problem is
        built once per run, then cached): call it off the GUI thread. None when it cannot be made."""
        if self.applied is None or self.data is None:
            return None
        run = _resolve(self.applied["run"]).resolve()
        if self.fit_curve is not None and self.fit_curve["run"] == str(run):
            return self.fit_curve
        if not (run / "monitor" / "status.json").exists():
            raise ValueError(f"{run} has no monitor record (fit run with --monitor off)")
        sys.path.insert(0, str(ROOT / "scripts"))
        from fit_monitor import SpectrumCache, read_run
        if self._fit_problems is None:
            self._fit_problems = SpectrumCache()
        prob = self._fit_problems.get(str(run))
        ids = [e["id"] for e in prob.series]
        label = self.applied.get("spectrum") if self.applied.get("spectrum") in ids else self.data["label"]
        if label not in ids:
            raise ValueError(f"the fit has no spectrum {label!r} (it has {ids})")
        k = ids.index(label)
        z = self.applied["fit"].get("z_final")
        source = "final vector (fit.json)"
        if z is None:
            starts = [v for v in read_run(run)["starts"].values() if v["z"] is not None and v["best"] is not None]
            if not starts:
                raise ValueError("no best point in the monitor record")
            best = min(starts, key=lambda v: v["best"])
            z, source = best["z"], f"best point of the monitor record (objective {best['best']:.5g})"
        z = np.asarray(z, float)
        joint = prob.joint
        fw = joint.forwards[k]
        model = fw.predict(joint.spectrum_vector(z, k)).model
        smooth = getattr(joint, "peak_smooth", 0.0)
        try:                                          # the objective as the fit ranks and reports it (hard
            if getattr(joint, "peaks", None) is not None:   # missing-peak rows): equals fit.json scores[0]
                joint.peak_smooth = 0.0
            objective = float(np.sum(joint.residual(z) ** 2))
        finally:
            joint.peak_smooth = smooth
        curve = {"run": str(run), "spectrum": label, "f": np.asarray(fw.f, float), "model": np.asarray(model, complex),
                 "objective": objective, "source": source}
        with self.lock:
            self.fit_curve = curve
        self.log(f"fit model of {run.name} ({label}): {source}", "session")
        return curve

    def _fit_curve_shown(self) -> Optional[dict]:
        """The applied fit's exact model when the parameters are still the fit's and it was computed."""
        c = self.fit_curve
        if c is None or self.applied is None or self.data is None or not self.state_is_applied_fit():
            return None
        if c["run"] != str(_resolve(self.applied["run"]).resolve()) or c["spectrum"] != self.data["label"]:
            return None
        return c

    def _snapshot(self):
        model_spec = {k: v for k, v in self.spec.items() if k != "molecule"}      # drawing only (D61)
        return (json.dumps(model_spec, sort_keys=True), self.exchange,
                tuple(sorted((c["key"], round(c["value"], 9)) for c in self.couplings())),
                tuple(round(v, 9) for v in self.field_nt), round(self.rate_per_s, 9))

    def state_is_applied_fit(self) -> bool:
        return self.applied is not None and self.applied["state"] == self._snapshot()

    # ---- figures (scripts/paper_figure.py) -----------------------------------------------
    def _manual_fit_file(self, out: Path) -> Path:
        """A fit.json-compatible file of the current sliders (the format load_fit / --from-joint read): couplings,
        field, decay rates (the applied fit's families when the rate slider was not moved) and its delay."""
        base = self.applied["fit"] if self.applied else {}
        sid = self.data["label"]
        old = dict(next(iter(base.get("spectrum_parameters", {}).values()), {})) if base else {}
        rates = {k: v for k, v in old.items() if ".log_rate" in k}
        moved_rate = not self.applied or abs(self.applied["state"][4] - round(self.rate_per_s, 9)) > 1e-12
        if moved_rate or not rates:
            _, model = self._model_cache
            edges = base.get("family_edges_hz") or []
            rates = {f"c{c}.log_rate{f}": math.log(self.rate_per_s)
                     for c in range(len(model.component_labels)) for f in range(len(edges) + 1)}
        sp = {k: v for k, v in old.items() if k == "phase_delay"}
        sp.update(rates)
        sp.update({"field_transverse_ut": 1e-3 * self.field_nt[0], "field_z_ut": 1e-3 * self.field_nt[1]})
        fit = {"x": [1.0], "manual_parameters": True, "family_edges_hz": base.get("family_edges_hz", []),
               "couplings": {c["key"]: {"J_at_x": [c["value"]]} for c in self.couplings()},
               "spectrum_parameters": {sid: sp}}
        path = out / "parameters_fit.json"
        path.write_text(json.dumps(fit, indent=1))
        return path

    def figure_command(self, **options) -> dict:
        """paper_figure.py command for the current state: the applied fit's run unchanged, or (sliders moved) a
        parameter file of the current state, labelled as manual parameters. Options: title, wide, segments, gains,
        display_window, colors (JSON), insets (path), formats, dpi, fid, out."""
        if self.data is None or not self.data.get("series"):
            raise ValueError("load a series file before making a figure")
        fid = options.get("fid") or self.data.get("source_fid")
        if not fid:
            raise ValueError("the series entry has no source_fid: give fid (the averaged FID .npy)")
        label = re.sub(r"[^A-Za-z0-9]+", "-", self.data["label"])
        out = Path(options.get("out") or self.workspace / "figures" / f"{time.strftime('%Y%m%d-%H%M%S')}_{label}")
        out.mkdir(parents=True, exist_ok=True)
        manual = not self.state_is_applied_fit()
        if self.applied:                          # the problem (options) of the applied fit's run
            if str(ROOT / "scripts") not in sys.path:
                sys.path.insert(0, str(ROOT / "scripts"))
            from run_problem import problem_argv, run_argv
            base = problem_argv(run_argv(self.applied["run"]))
        else:                                     # the problem a fit from the current state would solve
            full = self.fit_command(fit_field=any(self.field_nt))
            base = full[2:full.index("--out")]
        fit_path = self._manual_fit_file(out) if manual else Path(self.applied["run"]) / "fit.json"
        base = self._strip(base, ("--couplings", "--trace", "--starts", "--workers", "--max-nfev"))
        argv = [sys.executable, str(ROOT / "scripts" / "paper_figure.py"), *base, "--fit", str(fit_path),
                "--fid", str(fid), "--figure", str(out / "figure.png"),
                "--formats", str(options.get("formats", "png,pdf,svg")), "--dpi", str(int(options.get("dpi", 300)))]
        for key, flag in (("title", "--title"), ("wide", "--wide"), ("segments", "--segments"), ("gains", "--gains"),
                          ("display_window", "--display-window"), ("colors", "--colors"), ("insets", "--insets")):
            if options.get(key) not in (None, ""):
                argv += [flag, str(options[key])]
        if manual:
            argv.append("--manual")
        return {"argv": argv, "out_dir": str(out), "manual": manual, "fit_file": str(fit_path)}

    @staticmethod
    def _strip(argv, flags):
        out, skip = [], False
        for a in argv:
            if skip:
                skip = False
                continue
            if a in flags:
                skip = True
                continue
            out.append(a)
        return out

    def make_figure(self, **options) -> dict:
        """Run paper_figure.py in the background (PNG, PDF, SVG and caption.txt); poll figure_status."""
        if self.figure_job is not None and self.figure_job.running:
            raise RuntimeError("a figure is already being made")
        cmd = self.figure_command(**options)
        out = Path(cmd["out_dir"])
        self.figure = {"directory": str(out), "manual": cmd["manual"], "files": []}
        self.log(("figure (manual parameters, not a fit): " if cmd["manual"] else "figure of the applied fit: ")
                 + " ".join(cmd["argv"][1:]), "figure")

        def done(job):
            self.figure["files"] = sorted(p.name for p in out.glob("figure*"))
            self._changed("figure_done")
        self.figure_job = FitJob(cmd["argv"], out, self.log, done, source="figure", kind="figure",
                                 title=f"figure {out.name}")
        self.jobs.append(self.figure_job)
        return self.figure_status()

    def figure_status(self) -> dict:
        if self.figure_job is None:
            return {"running": False, "directory": None}
        st = self.figure_job.status()
        out = Path(st["out_dir"])
        return {"running": st["running"], "returncode": st["returncode"], "directory": str(out),
                "manual": self.figure["manual"], "seconds": st["seconds"],
                "files": sorted(p.name for p in out.glob("figure*")),
                "png": str(out / "figure.png") if (out / "figure.png").exists() else None}

    def export_figure(self, directory: str) -> dict:
        """Copy the last figure (PNG, PDF, SVG, caption, parameter file) to a directory."""
        import shutil
        st = self.figure_status()
        if not st.get("directory") or st["running"] or not st["files"]:
            raise ValueError("no finished figure to export")
        dest = _resolve(directory)
        dest.mkdir(parents=True, exist_ok=True)
        src = Path(st["directory"])
        copied = []
        for p in list(src.glob("figure*")) + list(src.glob("parameters_fit.json")):
            shutil.copy2(p, dest / p.name)
            copied.append(p.name)
        self.log(f"figure exported to {dest}: {', '.join(sorted(copied))}", "figure")
        return {"directory": str(dest), "files": sorted(copied)}

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

    # ---- custom spin systems (PLAN 8c) ------------------------------------------------------
    def spin_system_from_structure(self) -> dict:
        """Turn the current structure into an editable spin system (every isotopologue a component with its
        abundance weight and its numeric J matrix)."""
        from .spin_model import from_structure_model
        _, model = self._model_cache
        name = self.spec.get("compound") or self.spec.get("motif") or "custom"
        spec = from_structure_model(model, name)
        try:                                          # keep the molecule for drawing (D61)
            from zulf_hypothesis.molecule import molecule_record
            spec["molecule"] = molecule_record(self.spec)
        except Exception as exc:
            self.log(f"no molecule kept for the spin system: {exc}", "session")
        return self.set_structure(spec)

    def attach_molecule(self, molecule: str) -> dict:
        """Attach a molecule (SMILES or mol-file text) to the current model for drawing only (D61): the model, its
        couplings and an applied fit are unchanged. Sites are labelled C1, C2, O1, ... in atom order, so
        isotopologue names such as 13C@C1 mark their atom."""
        from zulf_hypothesis.molecule import molecule_record
        rec = molecule_record(molecule)
        with self.lock:
            self.spec = dict(self.spec, molecule=rec)
        self.log(f"molecule attached for drawing: {rec['smiles']}", "session")
        self._changed("molecule")
        return rec

    def load_model(self, path: str) -> dict:
        """A model file (a structure or spin-system specification as JSON) or a ZULF_NMR_Suite molecule folder."""
        p = _resolve(path)
        if p.is_dir() or p.name == "structure.csv":
            from .spin_model import read_suite_molecule
            spec = read_suite_molecule(p)
        else:
            spec = json.loads(p.read_text())
            spec = spec.get("structure", spec)
        return self.set_structure(spec)

    def save_model(self, path: str) -> dict:
        p = _resolve(path)
        if p.suffix != ".json":
            p = p.with_name(p.name + ".json")
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.spec, indent=1))
        self.log(f"model saved: {p}")
        return {"path": str(p)}

    # ---- task mode (D58) -----------------------------------------------------------------
    MODES = ("simulate", "process", "fit", "blind")

    def set_mode(self, mode: str) -> dict:
        """Task mode of the window: simulate (model only), process (scans to spectrum), fit (spectrum and a known
        model), blind (spectrum, unknown model). The model is the same in every mode."""
        if mode not in self.MODES:
            raise ValueError(f"mode must be one of {self.MODES}")
        if mode != self.mode:
            self.mode = mode
            self.log(f"mode: {mode}")
            self._changed("mode")
        return {"mode": self.mode}

    # ---- jobs and machine load -----------------------------------------------------------
    def job_list(self) -> List[dict]:
        """Every background job of the session, newest first, with its index for stop_job; then the fits running
        on this machine outside Studio (command line, other sessions) with their progress from RUN/monitor
        (index None: Studio does not stop them)."""
        own = [dict(job.status(), index=i, external=False) for i, job in reversed(list(enumerate(self.jobs)))]
        for st in own:                             # live evaluations and best objective while a start runs
            if st["kind"] == "fit" and st["running"]:
                st.update({k: v for k, v in monitor_progress(Path(st["out_dir"])).items() if v is not None})
        return own + self.external_runs()

    def external_runs(self) -> List[dict]:
        """fit_joint_series runs of other programs, read from their monitor files (cached for 3 s)."""
        cache = getattr(self, "_external_cache", None)
        if cache and time.time() - cache[0] < 3.0:
            return cache[1]
        own = {job.pid for job in self.jobs}
        out = []
        for p in analysis_processes():
            if p["pid"] in own or p["script"] != "fit_joint_series.py" or not p["out"]:
                continue
            run = _resolve(p["out"])
            st = {"kind": "fit", "title": f"{run.name} (outside Studio)", "pid": p["pid"], "workers": p["workers"],
                  "running": True, "returncode": None, "out_dir": str(run), "seconds": 0.0, "best_objective": None,
                  "starts_finished": 0, "starts_total": None, "phase": "running", "last_line": "",
                  "has_result": (run / "fit.json").exists(), "index": None, "external": True}
            st.update({k: v for k, v in monitor_progress(run).items() if v is not None})
            out.append(st)
        self._external_cache = (time.time(), out)
        return out

    def stop_job(self, index: int) -> dict:
        job = self.jobs[int(index)]
        job.stop()
        self.log(f"stop requested: {job.title}", job.source)
        return job.status()

    def monitor_url(self, run=None) -> dict:
        """URL of the live fit monitor page (scripts/fit_monitor.py: objective per start, current simulated
        spectrum, console). Studio serves it itself on 127.0.0.1 (a free port, started on first use, stopped with
        Studio); with run, the page opens on that run directory."""
        if getattr(self, "_monitor", None) is None:
            import threading
            from http.server import ThreadingHTTPServer
            sys.path.insert(0, str(ROOT / "scripts"))
            from fit_monitor import SpectrumCache, make_handler
            server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(str(ROOT / "runs"), SpectrumCache()))
            threading.Thread(target=server.serve_forever, daemon=True, name="fit-monitor").start()
            self._monitor = server
            self.log(f"fit monitor on http://127.0.0.1:{server.server_address[1]}", "studio")
        url = f"http://127.0.0.1:{self._monitor.server_address[1]}/"
        has_monitor = False
        if run:
            path = _resolve(run).resolve()
            has_monitor = (path / "monitor" / "status.json").exists()
            from urllib.parse import quote
            url += "?run=" + quote(str(path), safe="")
        return {"url": url, "run": None if not run else str(_resolve(run).resolve()), "has_monitor": has_monitor}

    def machine_status(self) -> dict:
        """Cores, load average and the analysis processes running on this machine (any program), with the
        workers they use and the free worker slots (cores minus busy workers). Fits use one BLAS thread per worker
        (ONE_THREAD); more workers than free cores slow every run down."""
        cores = cpu_cores()
        procs = analysis_processes()
        own = {job.pid for job in self.jobs if job.running}
        busy = sum(p["workers"] for p in procs)
        try:
            load = [round(v, 2) for v in os.getloadavg()]
        except OSError:
            load = None
        return {"cores": cores, "load_average": load, "processes": [dict(p, own=p["pid"] in own) for p in procs],
                "busy_workers": busy, "free_workers": max(cores["logical"] - busy, 0),
                "suggested_workers": max(1, min(8, cores["logical"] - busy))}

    # ---- blind analysis (scripts/analyze_sample.py) ----------------------------------------
    def blind_command(self, fid: Optional[str] = None, workers: int = 4, structure: Optional[dict] = None,
                      labeling: str = "", out: Optional[str] = None) -> List[str]:
        """analyze_sample command for an averaged FID (default: the FID the loaded series came from)."""
        fid = fid or (self.data or {}).get("source_fid")
        if not fid:
            raise ValueError("no FID: load a series made by make_series_entry (it records source_fid) or pass fid")
        path = _resolve(fid)
        if not path.exists():
            raise ValueError(f"FID not found: {path}")
        label = re.sub(r"[^A-Za-z0-9]+", "-", (self.data or {}).get("label") or path.parent.name).strip("-")
        outdir = Path(out) if out else self.workspace / "blind" / f"{time.strftime('%Y%m%d-%H%M%S')}_{label}"
        argv = [sys.executable, str(ROOT / "scripts" / "analyze_sample.py"), str(path), "--id", label,
                "--out", str(outdir), "--workers", str(int(workers))]
        if structure:
            argv += ["--structure", json.dumps(structure)]
            if labeling:
                argv += ["--labeling", labeling]
        return argv

    def start_blind(self, **options) -> dict:
        """Blind analysis (or a known-structure analysis with structure=...) of an averaged FID in the background:
        processing, hypotheses, search, report (OUT/blind.md or OUT/structure.md)."""
        if self.blind_job is not None and self.blind_job.running:
            raise RuntimeError("a blind analysis is already running")
        argv = self.blind_command(**options)
        out = Path(argv[argv.index("--out") + 1])
        self.log("blind analysis: " + " ".join(argv[1:]), "blind")

        def done(job):
            report = next((p for p in (out / "blind.md", out / "structure.md") if p.exists()), None)
            self.log(f"blind analysis report: {report}" if report else "blind analysis wrote no report", "blind")
            self._changed("blind_done")
        self.blind_job = FitJob(argv, out, self.log, done, source="blind", kind="blind",
                                title=f"blind {Path(argv[2]).parent.name}")
        self.jobs.append(self.blind_job)
        self._changed("blind_started")
        return self.blind_status()

    def blind_status(self) -> dict:
        if self.blind_job is None:
            return {"running": False, "out_dir": None}
        st = self.blind_job.status()
        out = Path(st["out_dir"])
        st["report"] = next((str(p) for p in (out / "blind.md", out / "structure.md") if p.exists()), None)
        return st

    # ---- import ------------------------------------------------------------------------------
    def import_fid(self, fid: str, record_s: float = 7.5, crop_s: float = 0.1, grid: str = "20,380",
                   exclude: str = "", label: str = "", sampling_rate: float = 0.0) -> dict:
        """An averaged FID (.npy) -> processed series (scripts/make_series_entry.py; sampling rate from scans.json
        or the .ini next to it, whole-grid fit ranges) in the workspace, loaded when it is done."""
        path = _resolve(fid)
        if not path.exists():
            raise ValueError(f"FID not found: {path}")
        label = re.sub(r"[^A-Za-z0-9]+", "-", label or f"{path.parent.name}-{path.stem}").strip("-")
        out = self.workspace / "series" / label
        argv = [sys.executable, str(ROOT / "scripts" / "make_series_entry.py"), "--fid", str(path), "--id", label,
                "--out", str(out), "--record", str(record_s), "--crop", str(crop_s), "--grid", grid]
        if exclude:
            argv += ["--exclude", exclude]
        if sampling_rate:
            argv += ["--sampling-rate", str(sampling_rate)]
        self.log("import FID: " + " ".join(argv[1:]), "import")

        def done(job):
            if job.returncode == 0 and (out / "series.json").exists():
                self.load_spectrum(series=str(out / "series.json"))
            self._changed("import_done")
        job = FitJob(argv, out, self.log, done, source="import", kind="import", title=f"import {path.name}")
        self.jobs.append(job)
        self._changed("import_started")
        return job.status()

    def import_scans(self, run_folder: str, out: Optional[str] = None, exclude_z: float = 0.0,
                     fid_options: Optional[dict] = None, then: str = "series") -> dict:
        """An instrument run folder (<n>.dat, <n>.ini; read-only) -> averaged FID (scripts/average_scans.py), then
        import_fid of the average. The average goes to `out` (default the workspace; averages kept for other work
        belong in ~/research/<project>/data/processed/<measurement>/)."""
        run = _resolve(run_folder)
        if not any(run.glob("*.dat")):
            raise ValueError(f"no <n>.dat scans in {run}")
        dest = _resolve(out) if out else self.workspace / "averages" / run.name
        argv = [sys.executable, str(ROOT / "scripts" / "average_scans.py"), str(run), str(dest)]
        if exclude_z:
            argv += ["--exclude-z", str(exclude_z)]
        self.log("average scans: " + " ".join(argv[1:]), "import")

        def done(job):
            if job.returncode == 0 and (dest / "average_fid.npy").exists() and then == "process":
                self.set_process_source(str(dest / "average_fid.npy"))
            elif job.returncode == 0 and (dest / "average_fid.npy").exists():
                opts = dict(fid_options or {})
                opts["label"] = opts.get("label") or run.name
                self.import_fid(str(dest / "average_fid.npy"), **opts)
            self._changed("import_done")
        job = FitJob(argv, dest, self.log, done, source="import", kind="import", title=f"average {run.name}")
        self.jobs.append(job)
        self._changed("import_started")
        return job.status()

    def open_path(self, path: str) -> dict:
        """Open whatever a path holds: a session file, a series file, a fit run (directory or its fit.json), a
        directory with series.json, an averaged FID (.npy: imported), or an instrument run folder (averaged)."""
        p = _resolve(path)
        if p.is_dir():
            if (p / "fit.json").exists():
                return {"opened": "fit run", **self.apply_fit(str(p))}
            if (p / "trace.json").exists():
                return {"opened": "trace", **self.load_trace(str(p))}
            if (p / "series.json").exists():
                return {"opened": "series", **self.load_spectrum(series=str(p / "series.json"))}
            if (p / "average_fid.npy").exists():
                return {"opened": "FID", **self.import_fid(str(p / "average_fid.npy"))}
            if any(p.glob("*.dat")):
                return {"opened": "scan folder", **self.import_scans(str(p))}
            raise ValueError(f"nothing to open in {p}")
        if p.suffix == ".npy":
            return {"opened": "FID", **self.import_fid(str(p))}
        if p.suffix == ".json" or p.suffix == SESSION_SUFFIX:
            content = json.loads(p.read_text())
            if isinstance(content, dict) and content.get("format") == SESSION_FORMAT:
                return {"opened": "session", **self.open_session(str(p))}
            if isinstance(content, dict) and "couplings" in content and "spectrum_parameters" in content:
                return {"opened": "fit run", **self.apply_fit(str(p.parent))}
            if isinstance(content, list) and content and "freq" in content[0]:
                return {"opened": "series", **self.load_spectrum(series=str(p))}
        raise ValueError(f"unknown file type: {p.name}")

    # ---- Process mode: an averaged FID and its recipe (PLAN 8b) ----------------------------------
    RECIPE = {"crop_s": 0.1, "record_s": 7.5, "apodization_per_s": 0.3, "zero_fill": 3, "grid": [20.0, 380.0],
              "sg_window_s": 0.0, "phase0_deg": None, "delay_ms": None, "exclude": "81.5,86"}

    def set_process_source(self, fid: str) -> dict:
        """The averaged FID to process (and, when average_scans.py made it, its per-scan record for the scan
        table). The recipe keeps its values; the sampling rate comes from scans.json or the .ini next to the FID."""
        from zulf_processing import find_sampling_rate
        path = _resolve(fid)
        y = np.asarray(np.load(path), float)
        fs, fs_source = find_sampling_rate(path, None)
        record = path.parent / "scans.json"
        scans = json.loads(record.read_text()) if record.exists() else None
        recipe = dict(self.RECIPE, **(self.process or {}).get("recipe", {}))
        recipe["record_s"] = min(recipe["record_s"], len(y) / fs - recipe["crop_s"])
        self.process = {"fid": str(path), "y": y, "fs": fs, "fs_source": fs_source, "scans": scans,
                        "run": (scans or {}).get("run"), "recipe": recipe, "preview": None}
        self.log(f"process source: {path} ({len(y)} points, {fs:g} Hz from {Path(fs_source).name})", "process")
        self.process_preview()
        self._changed("process_source")
        return self.process_status()

    def set_recipe(self, **values) -> dict:
        if self.process is None:
            raise ValueError("no FID to process (set_process_source)")
        for k, v in values.items():
            if k not in self.RECIPE:
                raise ValueError(f"unknown recipe key {k!r}")
            self.process["recipe"][k] = v
        self.process_preview()
        self._changed("process_recipe")
        return self.process_status()

    def process_preview(self) -> dict:
        """The processed spectrum of the current recipe (zulf_processing.series_spectrum, the same steps as a
        saved series) and the FID start, for the live preview."""
        from zulf_processing.series_spectrum import series_spectrum
        p = self.process
        r = p["recipe"]
        sp = series_spectrum(p["y"], p["fs"], r["crop_s"], r["record_s"], r["apodization_per_s"], r["zero_fill"],
                             r["grid"], sg_window_s=r["sg_window_s"] or None, phase0_deg=r["phase0_deg"],
                             delay_ms=r["delay_ms"])
        n0 = min(int(0.15 * p["fs"]), len(p["y"]))
        p["preview"] = {"f": sp["f"], "spectrum": sp["spectrum"], "t_ms": np.arange(n0) / p["fs"] * 1e3,
                        "fid": p["y"][:n0], "edge_ms": 1e3 * sp["edge_s"],
                        "phase0_deg": float(np.degrees(sp["phasing"]["phase0_rad"])),
                        "delay_ms": 1e3 * sp["phasing"]["delay_s"]}
        return p["preview"]

    def process_status(self) -> dict:
        if self.process is None:
            return {"source": None}
        p = self.process
        pv = p["preview"] or {}
        return {"source": p["fid"], "sampling_rate_hz": p["fs"], "points": len(p["y"]), "recipe": dict(p["recipe"]),
                "edge_ms": pv.get("edge_ms"), "phase0_deg": pv.get("phase0_deg"), "delay_ms": pv.get("delay_ms"),
                "scans": None if p["scans"] is None else {"found": p["scans"].get("scans_found"),
                                                          "kept": p["scans"].get("scans_kept"), "run": p["run"]}}

    def average_selection(self, keep: List[int], exclude_z: float = 0.0, out: Optional[str] = None) -> dict:
        """Average the chosen scans of the source run again (average_scans.py --keep; the run stays read only),
        then make the new average the process source."""
        if self.process is None or not self.process.get("run"):
            raise ValueError("the process source has no scan record (average a scan folder first)")
        run = Path(self.process["run"])
        dest = _resolve(out) if out else self.workspace / "averages" / f"{run.name}_{time.strftime('%Y%m%d-%H%M%S')}"
        dest.mkdir(parents=True, exist_ok=True)
        keep_file = dest / "keep.json"
        keep_file.write_text(json.dumps(sorted(int(k) for k in keep)))
        argv = [sys.executable, str(ROOT / "scripts" / "average_scans.py"), str(run), str(dest), "--keep",
                str(keep_file)] + (["--exclude-z", str(exclude_z)] if exclude_z else [])
        self.log(f"average {len(keep)} chosen scans: " + " ".join(argv[1:]), "process")

        def done(job):
            if job.returncode == 0 and (dest / "average_fid.npy").exists():
                self.set_process_source(str(dest / "average_fid.npy"))
            self._changed("import_done")
        job = FitJob(argv, dest, self.log, done, source="process", kind="import", title=f"average {len(keep)} scans")
        self.jobs.append(job)
        self._changed("import_started")
        return job.status()

    def save_processed(self, label: str = "", out: Optional[str] = None, load: bool = True) -> dict:
        """The processed spectrum of the current recipe as a series (make_series_entry.py with the same recipe,
        whole-grid fit ranges minus the excluded bands) and the recipe as recipe.json beside it; loaded when done."""
        if self.process is None:
            raise ValueError("no FID to process")
        p, r = self.process, self.process["recipe"]
        src = Path(p["fid"])
        label = re.sub(r"[^A-Za-z0-9]+", "-", label or f"{src.parent.name}-{src.stem}").strip("-")
        dest = _resolve(out) if out else self.workspace / "series" / label
        argv = [sys.executable, str(ROOT / "scripts" / "make_series_entry.py"), "--fid", str(src), "--id", label,
                "--out", str(dest), "--record", str(r["record_s"]), "--crop", str(r["crop_s"]),
                "--apodization", str(r["apodization_per_s"]), "--zero-fill", str(int(r["zero_fill"])),
                "--grid", f"{r['grid'][0]},{r['grid'][1]}", "--sampling-rate", str(p["fs"])]
        if r.get("exclude"):
            argv += ["--exclude", r["exclude"]]
        if r["sg_window_s"]:
            argv += ["--sg-window", str(r["sg_window_s"])]
        if r["phase0_deg"] is not None:
            argv += ["--phase0-deg", str(r["phase0_deg"])]
        if r["delay_ms"] is not None:
            argv += ["--delay-ms", str(r["delay_ms"])]
        dest.mkdir(parents=True, exist_ok=True)
        (dest / "recipe.json").write_text(json.dumps({"source_fid": str(src), "sampling_rate_hz": p["fs"],
                                                      "recipe": r}, indent=1))
        self.log("save processed spectrum: " + " ".join(argv[1:]), "process")

        def done(job):
            if load and job.returncode == 0 and (dest / "series.json").exists():
                self.load_spectrum(series=str(dest / "series.json"))
            self._changed("import_done")
        job = FitJob(argv, dest, self.log, done, source="process", kind="import", title=f"save {label}")
        self.jobs.append(job)
        self._changed("import_started")
        return job.status()

    # ---- inspection before an import ----------------------------------------------------------
    def inspect_path(self, path: str) -> dict:
        """What a path holds and what an import would use, without loading or writing anything: kind (scan folder,
        FID, series, session, fit run), sampling rate and where it comes from, points, duration, scan count and
        time span, cached averages of other tools (halp_compiled.npy), problems; plus a preview (the FID start and
        its magnitude spectrum) as arrays."""
        from zulf_processing import find_sampling_rate, read_settings
        p = _resolve(path)
        info = {"path": str(p), "kind": "unknown", "rows": [], "warnings": [], "preview": None, "cached": None}
        rows, warn = info["rows"], info["warnings"]
        fid = None
        if p.is_dir() and any(p.glob("*.dat")):
            info["kind"] = "scan folder"
            scans = sorted(p.glob("*.dat"), key=lambda q: int(q.stem) if q.stem.isdigit() else 10 ** 9)
            numbered = [q for q in scans if q.stem.isdigit()]
            inis = {q.stem for q in p.glob("*.ini")}
            missing = [q.stem for q in numbered if q.stem not in inis]
            first = read_settings(p / f"{numbered[0].stem}.ini") if numbered and numbered[0].stem in inis else {}
            last = read_settings(p / f"{numbered[-1].stem}.ini") if numbered and numbered[-1].stem in inis else {}
            fs = first.get("sampling_rate_hz")
            sizes = {q.stat().st_size for q in numbered[:50]}
            rows += [("scans", f"{len(numbered)} (numbers {numbered[0].stem}-{numbered[-1].stem})" if numbered else "0"),
                     ("sampling rate", f"{fs:g} Hz (0.ini)" if fs else "not in the .ini"),
                     ("points per scan", str(first.get("points", "?"))),
                     ("record", f"{first['points'] / fs:.2f} s" if fs and first.get("points") else "?"),
                     ("sequence", Path(str(first.get("pulse_sequence", "?")).replace("\\", "/")).name),
                     ("first / last scan", f"{first.get('date_time', '?')} / {last.get('date_time', '?')}"),
                     ("size", f"{sum(q.stat().st_size for q in numbered) / 1e9:.2f} GB")]
            if missing:
                warn.append(f"{len(missing)} scans without .ini (e.g. {', '.join(missing[:5])})")
            if len(sizes) > 1:
                warn.append(f"scan files differ in size ({sorted(sizes)[:4]} bytes): different record lengths?")
            cached = p / "halp_compiled.npy"
            if cached.exists():
                arr = np.load(cached, mmap_mode="r")
                info["cached"] = str(cached)
                rows.append(("cached average", f"halp_compiled.npy {arr.shape}, made by another tool (which scans "
                                               "and how is not recorded)"))
                fid = np.asarray(arr, float)
            else:
                from zulf_core.io import decode_dat
                fid = np.asarray(decode_dat(numbered[0]), float) if numbered else None
                rows.append(("preview", f"scan {numbered[0].stem} (one scan)" if numbered else "none"))
        elif p.suffix == ".npy":
            info["kind"] = "FID"
            fid = np.asarray(np.load(p), float)
            fs, source = find_sampling_rate(p, None) if (p.parent / "scans.json").exists() or \
                any(p.parent.glob("*.ini")) else (None, "")
            rows += [("points", str(len(fid))),
                     ("sampling rate", f"{fs:g} Hz ({Path(source).name})" if fs else "not found next to the file: "
                      "give it in the import (a wrong rate scales every frequency)"),
                     ("record", f"{len(fid) / fs:.2f} s" if fs else "?")]
            rec = p.parent / "scans.json"
            if rec.exists():
                meta = json.loads(rec.read_text())
                rows.append(("average of", f"{meta.get('scans_kept')} of {meta.get('scans_found')} scans "
                                           f"(exclude z {meta.get('exclude_z')}); {meta.get('run', '')}"))
            if not fs:
                warn.append("no sampling rate next to the FID (scans.json or .ini)")
        elif p.is_dir() or p.suffix in (".json", SESSION_SUFFIX):
            info["kind"] = "series, session or fit run (opened directly)"
            return info
        else:
            warn.append("unknown file type")
            return info
        if fid is not None and len(fid) > 16:
            fs = next((float(r[1].split()[0]) for r in rows if r[0] == "sampling rate" and r[1][0].isdigit()), None)
            if fs:
                from zulf_processing.diagnostics import switching_edge
                try:
                    edge = switching_edge(fid, fs)["edge_time_s"]
                    rows.append(("switching edge", f"{1e3 * edge:.2f} ms"))
                except ValueError:
                    warn.append("no switching edge found in the FID start")
                start = int(0.1 * fs)
                seg = fid[start:start + int(min(8.0, (len(fid) - start) / fs) * fs)]
                from scipy.signal import savgol_filter          # drift out, as in the processing (SG low-pass)
                win = max(int(0.05 * fs) // 2 * 2 + 1, 5)
                seg = seg - savgol_filter(seg, win, 2) if len(seg) > win else seg - seg.mean()
                spec = np.abs(np.fft.rfft(seg * np.exp(-0.3 * np.arange(len(seg)) / fs), 4 * len(seg)))
                freq = np.fft.rfftfreq(4 * len(seg), 1 / fs)
                keep = (freq >= 20.0) & (freq <= min(400.0, fs / 2))      # above the drift region
                n0 = int(0.1 * fs)
                info["preview"] = {"t_ms": (np.arange(n0) / fs * 1e3).tolist(), "fid": fid[:n0].tolist(),
                                   "f": freq[keep][::4].tolist(), "mag": spec[keep][::4].tolist()}
        return info

    # ---- export bundles ---------------------------------------------------------------------------
    def export_bundle(self, directory: str, spectrum=("data", "simulation", "residual"), view_only: bool = False,
                      formats=("csv",), fid: bool = False, parameters: bool = True, fit: bool = True) -> dict:
        """One folder that can be shared: spectrum.csv (frequency_hz and, per chosen source, real, imaginary,
        magnitude; header line), spectrum.npz (complex, with csv or alone), fid.csv (time_s, signal of the source
        FID), parameters (parameters.json, couplings.csv, lines.csv, the session file), the applied fit
        (fit.json, J_table.csv), and information.json: sources with sha256, time, settings and the code commit."""
        import hashlib
        import shutil
        out = _resolve(directory)
        out.mkdir(parents=True, exist_ok=True)
        files, sources = [], {}

        def sha(path):
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for block in iter(lambda: fh.read(1 << 20), b""):
                    h.update(block)
            return h.hexdigest()
        if parameters:                       # first: export() writes its own spectrum.csv, replaced below
            r = self.export(str(out))
            (out / "spectrum.csv").unlink()
            files += [f for f in r["files"] if f != "spectrum.csv" and (fit or not f.startswith("applied_"))]
            if not fit:
                for name in ("applied_fit.json", "applied_J_table.csv"):
                    if (out / name).exists():
                        (out / name).unlink()
        fq = None
        if spectrum:
            if self.data is not None:
                f = self.data["freq"]
                view = list(self.view) if view_only else [float(f.min()), float(f.max())]
                sim = self.simulate(points=len(f), view=view)
            else:
                sim = self.simulate(points=4000, view=list(self.view))
            fq = np.asarray(sim["f"])
            cols = {"frequency_hz": fq}
            arrays = {"frequency_hz": fq}
            series = {"simulation": np.asarray(sim["sim_re"]) + 1j * np.asarray(sim["sim_im"])}
            if "data_re" in sim:
                series["data"] = np.asarray(sim["data_re"]) + 1j * np.asarray(sim["data_im"])
                series["residual"] = series["data"] - series["simulation"]
            for name in spectrum:
                if name in series:
                    z = series[name]
                    cols.update({f"{name}_real": z.real, f"{name}_imaginary": z.imag, f"{name}_magnitude": np.abs(z)})
                    arrays[name] = z
            if "csv" in formats:
                with open(out / "spectrum.csv", "w", newline="") as fh:
                    w = csv.writer(fh)
                    w.writerow(list(cols))
                    for row in zip(*cols.values()):
                        w.writerow([f"{v:.10g}" for v in row])
                files.append("spectrum.csv")
            if "npz" in formats:
                np.savez(out / "spectrum.npz", **arrays)
                files.append("spectrum.npz")
        if fid and self.data is not None and self.data.get("source_fid"):
            from zulf_processing import find_sampling_rate
            src = _resolve(self.data["source_fid"])
            y = np.load(src)
            fs, fs_source = find_sampling_rate(src, None)
            with open(out / "fid.csv", "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["time_s", "signal"])
                for k, v in enumerate(y):
                    w.writerow([f"{k / fs:.9g}", f"{float(v):.10g}"])
            files.append("fid.csv")
            sources["fid"] = {"path": str(src), "sha256": sha(src), "sampling_rate_hz": fs,
                              "sampling_rate_from": fs_source}
        if self.data is not None and self.data.get("series"):
            sp = Path(self.data["series"])
            sources["series"] = {"path": str(sp), "sha256": sha(sp)}
        if self.applied:
            sources["applied_fit"] = {"path": self.applied["run"]}
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=str(ROOT), capture_output=True,
                                text=True).stdout.strip()
        info = {"exported": time.strftime("%Y-%m-%dT%H:%M:%S"), "program": "ZULF Studio", "commit": commit,
                "range_hz": None if fq is None else [float(fq.min()), float(fq.max())],
                "spectrum_sources": list(spectrum), "simulation": "quick look: Lorentzian lines, one decay rate "
                "(not the processed fit model; the fit's own model is in the fit run)", "sources": sources,
                "state": self.session_dict(), "files": sorted(set(files))}
        (out / "information.json").write_text(json.dumps(info, indent=1))
        files.append("information.json")
        self.log(f"export bundle {out}: {', '.join(sorted(set(files)))}")
        return {"directory": str(out), "files": sorted(set(files))}

    # ---- session files -------------------------------------------------------------------------
    def session_dict(self) -> dict:
        with self.lock:
            return {"format": SESSION_FORMAT, "version": 1, "saved": time.strftime("%Y-%m-%dT%H:%M:%S"),
                    "mode": self.mode,
                    "structure": self.spec, "overrides": dict(self.overrides), "exchange": self.exchange,
                    "field_nt": list(self.field_nt), "rate_per_s": self.rate_per_s, "view": self.view,
                    "display": {"part": self.display, "phase_deg": self.data_phase_deg,
                                "delay_ms": self.data_delay_ms, "scale_lock": self.scale_lock},
                    "data": None if self.data is None else {"series": self.data.get("series"),
                                                            "index": self.data.get("index", 0)},
                    "applied_fit": self.applied["run"] if self.applied else None}

    def save_session(self, path: str) -> dict:
        """Structure, couplings, field, line width, view, display and the paths of the loaded series and applied
        fit as JSON (data and fits stay where they are; the file refers to them)."""
        p = _resolve(path)
        if p.suffix != SESSION_SUFFIX:
            p = p.with_name(p.name + SESSION_SUFFIX)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.session_dict(), indent=1))
        self.session_file = str(p)
        self.log(f"session saved: {p}")
        self._changed("session_saved")
        return {"path": str(p)}

    def open_session(self, path: str) -> dict:
        """Restore a saved session; unknown keys are ignored, missing ones keep their defaults, and a series or
        fit that no longer exists is reported, not fatal."""
        p = _resolve(path)
        d = json.loads(p.read_text())
        if d.get("format") != SESSION_FORMAT:
            raise ValueError(f"not a ZULF Studio session: {p}")
        missing = []
        with self.lock:
            self.spec = d.get("structure") or self.spec
            self.overrides = {k: float(v) for k, v in (d.get("overrides") or {}).items()}
            self.exchange = d.get("exchange", self.exchange)
            self.field_nt = [float(v) for v in d.get("field_nt", self.field_nt)]
            self.rate_per_s = float(d.get("rate_per_s", self.rate_per_s))
            disp = d.get("display") or {}
            self.display = disp.get("part", self.display)
            self.data_phase_deg = float(disp.get("phase_deg", self.data_phase_deg))
            self.data_delay_ms = float(disp.get("delay_ms", self.data_delay_ms))
            self.scale_lock = disp.get("scale_lock")
            self._build()
        data = d.get("data") or {}
        if data.get("series"):
            if _resolve(data["series"]).exists():
                self.load_spectrum(series=data["series"], index=int(data.get("index", 0)))
            else:
                missing.append(data["series"])
        if d.get("view"):
            self.view = [float(v) for v in d["view"]]
        run = d.get("applied_fit")
        if run and (_resolve(run) / "fit.json").exists():
            self.applied = {"run": run, "fit": json.loads((_resolve(run) / "fit.json").read_text()),
                            "spectrum": None, "state": self._snapshot()}
        elif run:
            missing.append(run)
        self.session_file = str(p)
        if d.get("mode") in self.MODES:
            self.mode = d["mode"]
        self.log(f"session opened: {p}" + (f"; missing: {', '.join(missing)}" if missing else ""))
        self._changed("structure")
        self._changed("session_opened")
        return {"path": str(p), "missing": missing}

    # ---- state and export ----------------------------------------------------------------
    def state(self) -> dict:
        with self.lock:
            return {"version": self.version, "structure": self.spec, "exchange": self.exchange,
                    "components": self.components(), "couplings": self.couplings(),
                    "field_nt": {"transverse": self.field_nt[0], "z": self.field_nt[1],
                                 "magnitude": math.hypot(*self.field_nt)},
                    "rate_per_s": self.rate_per_s, "view": self.view,
                    "display": {"part": self.display, "phase_deg": self.data_phase_deg,
                                "delay_ms": self.data_delay_ms, "scale_lock": self.scale_lock},
                    "data": None if self.data is None else {k: self.data[k] for k in ("label", "series", "index",
                                                                                       "ranges", "source_fid")},
                    "fit": self.fit_status(), "figure": self.figure_status(), "blind": self.blind_status(),
                    "jobs_running": sum(job.running for job in self.jobs), "session_file": self.session_file,
                    "mode": self.mode,
                    "parameters_are_applied_fit": self.state_is_applied_fit(),
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
            w.writerow(["component", "frequency_hz", "amplitude", "relative", "decay_per_s", "decay_time_s"])
            for r in self.lines(min_relative=1e-4):
                w.writerow([r["component"], f"{r['frequency_hz']:.6f}", f"{r['amplitude']:.6g}", f"{r['relative']:.6g}",
                            f"{r['decay_per_s']:.6g}", f"{1.0 / r['decay_per_s']:.6g}" if r["decay_per_s"] > 0 else ""])
        sim = self.simulate()
        with open(out / "spectrum.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            cols = ["frequency_hz", "sim_re", "sim_im"] + (["data_re", "data_im"] if "data_re" in sim else [])
            w.writerow(cols)
            for row in zip(*(sim[c] if c != "frequency_hz" else sim["f"] for c in cols)):
                w.writerow([f"{v:.8g}" for v in row])
        with open(out / "couplings.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["coupling", "J_hz", "set_by_hand"])
            for c in self.couplings():
                w.writerow([c["key"], f"{c['value']:.6f}", int(c["overridden"])])
        (out / ("session" + SESSION_SUFFIX)).write_text(json.dumps(self.session_dict(), indent=1))
        files = ["parameters.json", "lines.csv", "spectrum.csv", "couplings.csv", "session" + SESSION_SUFFIX]
        if self.applied and (Path(self.applied["run"]) / "fit.json").exists():
            import shutil
            for name in ("fit.json", "J_table.csv"):
                src = Path(self.applied["run"]) / name
                if src.exists():
                    shutil.copy2(src, out / f"applied_{name}")
                    files.append(f"applied_{name}")
        self.log(f"exported to {out}: {', '.join(files)}")
        return {"directory": str(out), "files": files}
