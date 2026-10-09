"""JSON tool API of the studio for AI agents (and scripts), over local HTTP.

    GET  /api/tools                  tool list with JSON-schema parameters (?format=anthropic|openai)
    GET  /api/state                  the session state
    POST /api/call   {"tool": NAME, "args": {...}}
    POST /api/<NAME> {...args...}    the same, one route per tool

Every call goes to the StudioSession the window shows, so a person and an agent work on one state; the window
redraws when an agent changes it. The server binds to 127.0.0.1 by default. Results of fits are conditional
numerical results (AGENTS.md): report them with their objective, residual and the settings used.
"""
from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Callable, Dict
from urllib.parse import parse_qs, urlparse

import numpy as np

from .session import FIT_DEFAULTS, StudioSession

NUM = {"type": "number"}
STR = {"type": "string"}


def _obj(props: dict, required=()) -> dict:
    return {"type": "object", "properties": props, "required": list(required)}


_JSON_TYPES = {bool: "boolean", int: "integer", float: "number", str: "string", list: "array"}
FIT_OPTIONS = {k: {"type": _JSON_TYPES[type(v)]} for k, v in FIT_DEFAULTS.items()}
FIT_OPTIONS["range"] = {"type": "array", "items": NUM}
FIT_OPTIONS["out"] = STR
FIGURE_OPTIONS = {"title": STR, "wide": {"type": "string", "description": "lo,hi Hz of panel a"},
                  "segments": {"type": "string", "description": "lo,hi;lo,hi detail panels"},
                  "gains": {"type": "string", "description": "display factor per detail panel, e.g. 1,2.5"},
                  "display_window": NUM, "colors": {"type": "string", "description": "JSON {component: colour}"},
                  "insets": {"type": "string", "description": "insets JSON file (needs RDKit)"},
                  "formats": {"type": "string", "description": "png,pdf,svg"}, "dpi": {"type": "integer"},
                  "fid": {"type": "string", "description": "averaged FID .npy when the series has no source_fid"},
                  "out": STR}


# name -> (session method name, description, parameter schema)
TOOLS: Dict[str, tuple] = {
    "get_state": ("state", "Current structure, couplings (Hz), field (nT), decay rate, view, loaded data, fit and "
                           "trace status.", _obj({})),
    "list_motifs": ("list_motifs", "Registered structural motifs with typical one-bond couplings.", _obj({})),
    "set_structure": ("set_structure", "Set the molecule. spec is a structure specification: {'motif': NAME, "
                      "'one_bond': {SITE: Hz}} or {'compound': NAME, 'chain': {'groups': [[SITE, ELEMENT, N_H], ...], "
                      "'bonds': [[A, B], ...]}, 'one_bond': {...}}. Resets coupling overrides unless keep_couplings.",
                      _obj({"spec": {"type": "object"}, "keep_couplings": {"type": "boolean"}}, ["spec"])),
    "set_couplings": ("set_couplings", "Set couplings by key 'J(a,b)' in Hz (site or proton-group labels; a new key "
                      "adds a coupling, e.g. a long-range one that makes another isotopologue visible).",
                      _obj({"values": {"type": "object", "additionalProperties": NUM}}, ["values"])),
    "remove_coupling": ("remove_coupling", "Drop a coupling override (back to the structure's value).",
                        _obj({"key": STR}, ["key"])),
    "set_exchange": ("set_exchange", "Exchangeable N-H/O-H protons: 'fast' (decoupled, dropped) or 'slow' (kept).",
                     _obj({"mode": {"type": "string", "enum": ["fast", "slow"]}}, ["mode"])),
    "set_field": ("set_field", "Static field during evolution in nT: transverse_nt (perpendicular to the detection "
                  "axis) and z_nt (along it). Zero field is the default.", _obj({"transverse_nt": NUM, "z_nt": NUM})),
    "set_linewidth": ("set_linewidth", "Decay rate of every line in 1/s (Lorentzian FWHM = rate / pi Hz).",
                      _obj({"rate_per_s": NUM}, ["rate_per_s"])),
    "set_view": ("set_view", "Frequency window in Hz.", _obj({"lo": NUM, "hi": NUM}, ["lo", "hi"])),
    "set_display": ("set_display", "Shown part of the spectra ('re', 'im', 'abs') and an extra display phase of the "
                    "data (degrees) and delay (ms).", _obj({"part": STR, "phase_deg": NUM, "delay_ms": NUM})),
    "auto_phase": ("auto_phase", "Set the display phase (deg) and delay (ms) of the data automatically: method "
                   "'model' (best match to the current simulation; lines roughly in place) or 'data' (model-free, "
                   "peak phases against frequency). fit_delay=false keeps the delay.",
                   _obj({"method": {"type": "string", "enum": ["model", "data"]}, "fit_delay": {"type": "boolean"},
                         "delay_range_ms": {"type": "array", "items": NUM}})),
    "lock_scale": ("lock_scale", "Freeze the display scale of the simulation (lock=true; value optional, default the "
                   "current automatic scale) or release it (lock=false). Automatic: least squares on the magnitudes.",
                   _obj({"lock": {"type": "boolean"}, "value": NUM})),
    "lines": ("lines", "Line table: every transition of every isotopologue with frequency (Hz), amplitude (times the "
              "natural abundance) and relative amplitude.",
              _obj({"min_relative": NUM, "view": {"type": "array", "items": NUM}})),
    "simulate": ("simulate", "Simulated spectrum in the view (Lorentzian lines), the loaded data there, the matching "
                 "scale and the rms residual. Arrays: f, sim_re, sim_im, data_re, data_im.",
                 _obj({"points": {"type": "integer"}, "view": {"type": "array", "items": NUM}})),
    "load_spectrum": ("load_spectrum", "Load a processed spectrum: series (a series.json from "
                      "scripts/make_series_entry.py) and index, or freq and values (.npy paths).",
                      _obj({"series": STR, "index": {"type": "integer"}, "freq": STR, "values": STR, "label": STR})),
    "fit_command": ("fit_command", "The fit_joint_series command a fit would run (no side effects).", _obj(FIT_OPTIONS)),
    "start_fit": ("start_fit", "Start a complex fit of the loaded series (fit_joint_series) from the current couplings "
                  "and field, in the background. Options as fit_command. Poll fit_status.", _obj(FIT_OPTIONS)),
    "fit_status": ("fit_status", "Running / finished, best objective so far, starts finished, output directory.",
                   _obj({})),
    "stop_fit": ("stop_fit", "Stop the running fit.", _obj({})),
    "apply_fit": ("apply_fit", "Put the couplings, field and decay rate of a finished fit into the session and load "
                  "its trace (default: the last fit).", _obj({"run_dir": STR, "spectrum": STR})),
    "load_trace": ("load_trace", "Load the fit trace (trace.json/npz) of a run directory for the progress slider.",
                   _obj({"run_dir": STR}, ["run_dir"])),
    "trace_frame": ("trace_frame", "Show frame `index` of the loaded trace (couplings and objective of that point of "
                    "the fit); apply=true copies its couplings into the session.",
                    _obj({"index": {"type": "integer"}, "apply": {"type": "boolean"}}, ["index"])),
    "figure_command": ("figure_command", "The scripts/paper_figure.py command for the current state (no side effects). "
                       "If the parameters are an applied fit, its run is drawn; otherwise a parameter file of the "
                       "current state is written and the figure says 'manual parameters (not a fit)'.", _obj(FIGURE_OPTIONS)),
    "make_figure": ("make_figure", "Make the publication figure (paper_figure.py) in the background: PNG (dpi), PDF, "
                    "SVG and caption.txt under the workspace. Poll figure_status.", _obj(FIGURE_OPTIONS)),
    "figure_status": ("figure_status", "Running / finished, directory, files, whether it shows manual parameters.",
                      _obj({})),
    "export_figure": ("export_figure", "Copy the last figure (PNG, PDF, SVG, caption, parameter file) to a directory.",
                      _obj({"directory": STR}, ["directory"])),
    "export": ("export", "Write parameters.json, couplings.csv, lines.csv, spectrum.csv, a session file and the "
               "applied fit's fit.json / J_table.csv to a directory (default under the workspace).",
               _obj({"directory": STR})),
    "set_mode": ("set_mode", "Task mode of the window: simulate (model only), process (scans to spectrum), fit "
                 "(spectrum and a known model) or blind (blind analysis). The model is shared by every mode.",
                 _obj({"mode": {"type": "string", "enum": ["simulate", "process", "fit", "blind"]}}, ["mode"])),
    "set_process_source": ("set_process_source", "Process mode: the averaged FID (.npy) to process; its per-scan "
                           "record (scans.json next to it) fills the scan table.", _obj({"fid": STR}, ["fid"])),
    "set_recipe": ("set_recipe", "Process mode: change the recipe (crop_s, record_s, apodization_per_s, zero_fill, "
                   "grid [lo, hi], sg_window_s (0 = default), phase0_deg and delay_ms (null = calibration and "
                   "switching edge), exclude 'lo,hi;...'); the preview is recomputed.",
                   _obj({"crop_s": NUM, "record_s": NUM, "apodization_per_s": NUM, "zero_fill": {"type": "integer"},
                         "grid": {"type": "array", "items": NUM}, "sg_window_s": NUM, "phase0_deg": NUM,
                         "delay_ms": NUM, "exclude": STR})),
    "process_status": ("process_status", "Process mode: source FID, sampling rate, recipe, switching edge, phase "
                       "and delay, scan counts.", _obj({})),
    "average_selection": ("average_selection", "Process mode: average the chosen scans of the source run again "
                          "(average_scans.py --keep; the run stays read only) and process the new average.",
                          _obj({"keep": {"type": "array", "items": {"type": "integer"}}, "exclude_z": NUM,
                                "out": STR}, ["keep"])),
    "save_processed": ("save_processed", "Process mode: save the spectrum of the recipe as a series with "
                       "recipe.json (make_series_entry.py) and load it for Fit.",
                       _obj({"label": STR, "out": STR, "load": {"type": "boolean"}})),
    "machine_status": ("machine_status", "Cores (performance / efficiency), load average and the analysis processes "
                       "running on this machine with their workers; free and suggested workers. Check it before "
                       "starting a fit or analysis: more workers than cores slow every run down.", _obj({})),
    "job_list": ("job_list", "Every background job of the session (fits, blind analyses, imports, figures), newest "
                 "first: kind, title, running, starts finished / total, stage, elapsed s, best objective, output.",
                 _obj({})),
    "stop_job": ("stop_job", "Stop the job with this index (from job_list).", _obj({"index": {"type": "integer"}},
                                                                                       ["index"])),
    "start_blind": ("start_blind", "Blind analysis of an averaged FID in the background (scripts/analyze_sample.py: "
                    "processing, hypotheses, search, ranked report); fid defaults to the FID of the loaded series. "
                    "With structure (a structure specification) it fits that structure instead; labeling for that "
                    "case. Poll blind_status. The ranking is a list of conditional candidates.",
                    _obj({"fid": STR, "workers": {"type": "integer"}, "structure": {"type": "object"},
                          "labeling": {"type": "string", "enum": ["", "natural", "15N", "2H-exchange", "unknown"]}})),
    "blind_status": ("blind_status", "State of the last blind analysis: running, stage, elapsed, output directory and "
                     "report path (blind.md or structure.md).", _obj({})),
    "import_fid": ("import_fid", "Process an averaged FID (.npy) into a series (make_series_entry: sampling rate "
                   "from scans.json or the .ini next to it, whole-grid fit ranges) and load it when done.",
                   _obj({"fid": STR, "record_s": NUM, "crop_s": NUM, "grid": STR, "exclude": STR, "label": STR,
                         "sampling_rate": NUM}, ["fid"])),
    "import_scans": ("import_scans", "Average an instrument run folder (n.dat, n.ini; read-only) with "
                     "average_scans.py, then import the average.",
                     _obj({"run_folder": STR, "out": STR, "exclude_z": NUM}, ["run_folder"])),
    "open_path": ("open_path", "Open a path by what it holds: session file, series.json, fit run (directory or "
                  "fit.json), directory with series.json, averaged FID (.npy, imported), scan folder (averaged).",
                  _obj({"path": STR}, ["path"])),
    "save_session": ("save_session", "Save structure, couplings, field, line width, view, display and the paths of "
                     "the loaded series and applied fit (.zulfstudio JSON).", _obj({"path": STR}, ["path"])),
    "open_session": ("open_session", "Restore a saved session (.zulfstudio); missing series or fits are reported.",
                     _obj({"path": STR}, ["path"])),
    "read_log": ("read_log", "Last log entries (sources: session, api, fit, terminal).",
                 _obj({"n": {"type": "integer"}, "source": STR})),
}


def _jsonable(x):
    if isinstance(x, dict):
        return {k: _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    return x


class StudioAPI:
    def __init__(self, session: StudioSession):
        self.session = session

    def tools(self, fmt: str = "") -> list:
        rows = [{"name": n, "description": d, "parameters": p} for n, (_, d, p) in TOOLS.items()]
        if fmt == "anthropic":
            return [{"name": r["name"], "description": r["description"], "input_schema": r["parameters"]} for r in rows]
        if fmt == "openai":
            return [{"type": "function", "function": r} for r in rows]
        return rows

    def call(self, tool: str, args: dict = None):
        if tool not in TOOLS:
            raise KeyError(f"unknown tool '{tool}'; see /api/tools")
        args = dict(args or {})
        method = getattr(self.session, TOOLS[tool][0])
        if tool == "set_structure":
            return _jsonable(method(args.pop("spec"), **args))
        if tool == "set_couplings":
            return _jsonable(method(args["values"]))
        result = _jsonable(method(**args))
        if tool not in ("get_state", "list_motifs", "lines", "simulate", "fit_status", "read_log", "fit_command",
                        "figure_status", "figure_command"):
            self.session.log(f"{tool} {json.dumps(args)[:300]}", "api")
        return result


def make_handler(api: StudioAPI):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _send(self, code, body):
            data = json.dumps(body).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            url = urlparse(self.path)
            if url.path == "/api/tools":
                return self._send(200, api.tools(parse_qs(url.query).get("format", [""])[0]))
            if url.path == "/api/state":
                return self._send(200, api.call("get_state"))
            if url.path in ("/", "/api"):
                return self._send(200, {"service": "zulf_studio", "tools": "/api/tools", "call": "POST /api/call",
                                        "state": "/api/state"})
            self._send(404, {"error": "not found"})

        def do_POST(self):
            try:
                n = int(self.headers.get("Content-Length", 0))
                body = json.loads(self.rfile.read(n) or b"{}")
                path = urlparse(self.path).path
                if path == "/api/call":
                    return self._send(200, {"result": api.call(body["tool"], body.get("args", {}))})
                if path.startswith("/api/"):
                    return self._send(200, {"result": api.call(path[len("/api/"):], body)})
                self._send(404, {"error": "not found"})
            except Exception as exc:
                self._send(400, {"error": f"{type(exc).__name__}: {exc}"})
    return Handler


def serve(session: StudioSession, host: str = "127.0.0.1", port: int = 8766, background: bool = True):
    """Start the API server; returns it (server.server_address has the port). background: own daemon thread."""
    server = ThreadingHTTPServer((host, port), make_handler(StudioAPI(session)))
    if background:
        threading.Thread(target=server.serve_forever, daemon=True).start()
    session.log(f"API on http://{host}:{server.server_address[1]}/api/tools", "api")
    return server
