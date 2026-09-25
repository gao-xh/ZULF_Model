"""Tool definitions exposed to AI agents, the CLI and MCP.

Rules the tools enforce and describe:
* experimental inputs are read-only; outputs go to the workspace;
* processing recipes are explicit per dataset and returned with every result;
* refined candidates are conditional numerical results with flags, residuals
  and held-out scores, never a determined assignment.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

import numpy as np

from ..spec import ProblemSpec
from ..spinsystem import Interpretation, SpinSystem
from . import jobs
from .registry import REGISTRY, new_run_dir, workspace

# -- shared schema fragments ----------------------------------------------------------------

SYSTEM = {"type": "object", "description": "Spin system: nucleus symbols, symmetric J matrix in Hz (zero diagonal), "
                                           "optional zero-based magnetic-equivalence groups.",
          "properties": {"isotopes": {"type": "array", "items": {"type": "string"}},
                         "couplings_hz": {"type": "array", "items": {"type": "array", "items": {"type": "number"}}},
                         "groups": {"type": "array", "items": {"type": "array", "items": {"type": "integer"}}}},
          "required": ["isotopes", "couplings_hz"]}
INTERPRETATION = {"type": "object", "description": "One or more components (isotopologues) explaining a spectrum.",
                  "properties": {"components": {"type": "array", "items": {
                      "type": "object", "properties": {"system": SYSTEM, "contribution": {"type": "number"},
                                                       "label": {"type": "string"}},
                      "required": ["system"]}}},
                  "required": ["components"]}
ACQUISITION = {"type": "object", "description": "Processing recipe; every step optional. Omitted fields keep "
                                                "defaults (no crop, no SG, no mean removal).",
               "properties": {k: {"type": t} for k, t in (
                   ("sampling_rate_hz", "number"), ("points", "integer"), ("start_sample", "integer"),
                   ("stop_sample", "integer"), ("sg_window", "integer"), ("sg_order", "integer"),
                   ("remove_mean", "boolean"), ("time_origin_s", "number"), ("apodization_rate_per_s", "number"),
                   ("phase_reference", "string"))}}
FID_SOURCE = {"type": "object", "description": "Averaged FID: {average_npy, ini} or {fid_npy, sampling_rate_hz}.",
              "properties": {"average_npy": {"type": "string"}, "ini": {"type": "string"},
                             "fid_npy": {"type": "string"}, "sampling_rate_hz": {"type": "number"},
                             "label": {"type": "string"}}}
RANGES = {"type": "array", "description": "Sorted disjoint [low_hz, high_hz] bands.",
          "items": {"type": "array", "items": {"type": "number"}, "minItems": 2, "maxItems": 2}}


def _interp(data: dict) -> Interpretation:
    return Interpretation.from_dict({"components": [
        {"system": c["system"], "contribution": c.get("contribution", 1.0), "label": c.get("label", "")}
        for c in data["components"]]})


def _load_fid(src: dict):
    from ..io import load_average
    if src.get("average_npy"):
        return load_average(src["average_npy"], src["ini"])
    from ..io import ExperimentFID
    fid = np.load(src["fid_npy"], allow_pickle=False).astype(float)
    return ExperimentFID(fid, float(src["sampling_rate_hz"]), 0, {"fid_npy": src["fid_npy"]})


def _acquisition(exp, overrides: Optional[dict]):
    from ..render.acquisition import Acquisition
    data = {"sampling_rate_hz": exp.sampling_rate_hz, "points": exp.points}
    data.update({k: v for k, v in (overrides or {}).items() if k not in ("sampling_rate_hz", "points")})
    return Acquisition.from_dict(data)


def _plot(path: Path, curves, xlabel: str, ylabel: str, title: str) -> Optional[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for x, y, label in curves:
        ax.plot(x, y, lw=0.7, label=label)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=0.2)
    if len(curves) > 1:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return str(path)


# -- tools ------------------------------------------------------------------------------------

@REGISTRY.tool("describe_project",
               "Explain the ZULF_Model toolkit: problem spec, conventions, available tools and interpretation rules. "
               "Call first in a new session.",
               {"type": "object", "properties": {"problem": {"type": "string", "description": "Optional ProblemSpec JSON path."}}})
def describe_project(args: dict) -> dict:
    spec = ProblemSpec.load(args["problem"]) if args.get("problem") else ProblemSpec()
    return {"problem_spec": spec.to_dict(), "workspace": str(workspace()), "tools": REGISTRY.names(),
            "long_running_tools": [n for n, t in REGISTRY.tools.items() if t.long_running],
            "rules": ["Experimental inputs are read-only; outputs are written under the workspace.",
                      "Processing recipes differ per dataset: run diagnose_fid and compare candidate recipes on "
                      "held-out acquisitions before trusting parameters.",
                      "Refined candidates are conditional numerical results; report flags, band residuals and "
                      "held-out residuals; never claim a determined molecular assignment.",
                      "A global sign flip of all J is unobservable at zero field; relative signs are observable.",
                      "Use submit_job for long-running tools and poll get_job."],
            "conventions": "docs/CONVENTIONS.md"}


@REGISTRY.tool("simulate_transitions",
               "Zero-field transition list (frequencies in Hz, complex amplitudes) of one spin system under the "
               "sudden-drop protocol, computed by exact diagonalization in magnetic-equivalence sectors.",
               {"type": "object", "properties": {"system": SYSTEM, "max_lines": {"type": "integer", "default": 200}},
                "required": ["system"]})
def simulate_transitions(args: dict) -> dict:
    from ..physics import compute_transitions
    system = SpinSystem.from_dict(args["system"])
    tl = compute_transitions(system)
    order = np.argsort(-np.abs(tl.amplitudes))[: args["max_lines"]]
    out = new_run_dir("simulate")
    np.savez(out / "transitions.npz", frequencies_hz=tl.frequencies_hz, amplitudes=tl.amplitudes)
    return {"n_transitions": len(tl), "groups": [list(g) for g in system.groups],
            "strongest": [{"frequency_hz": float(tl.frequencies_hz[i]), "amplitude": float(abs(tl.amplitudes[i]))}
                          for i in sorted(order, key=lambda i: tl.frequencies_hz[i])],
            "arrays": str(out / "transitions.npz")}


@REGISTRY.tool("render_spectrum",
               "Render the spectrum of an interpretation. mode 'pure' uses a finite record with no processing, "
               "'continuous' gives infinite-record Voigt lines, 'processed' applies the given acquisition recipe.",
               {"type": "object", "properties": {
                   "interpretation": INTERPRETATION, "mode": {"type": "string", "enum": ["pure", "continuous", "processed"],
                                                               "default": "pure"},
                   "acquisition": ACQUISITION, "sampling_rate_hz": {"type": "number", "default": 1000.0},
                   "points": {"type": "integer", "default": 16384}, "rate_per_s": {"type": "number", "default": 1.0},
                   "gaussian_sigma_hz": {"type": "number", "default": 0.0}, "f_min_hz": {"type": "number", "default": 1.0},
                   "f_max_hz": {"type": "number", "default": 400.0}},
                "required": ["interpretation"]})
def render_spectrum(args: dict) -> dict:
    from ..physics import compute_transitions
    from ..render import Acquisition, ContinuousRenderer, Renderer, SpectrumGrid
    interp = _interp(args["interpretation"])
    base = Acquisition.pure(args["sampling_rate_hz"], args["points"])
    acq = base.with_processing(**(args.get("acquisition") or {})) if args["mode"] == "processed" else base
    grid = SpectrumGrid.for_acquisition(acq, args["f_min_hz"], args["f_max_hz"])
    renderer = ContinuousRenderer(acq.n / acq.sampling_rate_hz) if args["mode"] == "continuous" else Renderer(acq)
    total = np.zeros(len(grid), complex)
    for c in interp.components:
        total += renderer.render(compute_transitions(c.system), args["rate_per_s"], grid.frequencies_hz,
                                 c.contribution, 0.0, args["gaussian_sigma_hz"])
    out = new_run_dir("render")
    np.savez(out / "spectrum.npz", frequency_hz=grid.frequencies_hz, spectrum=total)
    figure = _plot(out / "spectrum.png", [(grid.frequencies_hz, np.abs(total), "|S|")], "Frequency (Hz)",
                   "Magnitude (a.u.)", f"Rendered spectrum ({args['mode']})")
    return {"arrays": str(out / "spectrum.npz"), "figure": figure, "acquisition": acq.to_dict(),
            "points": len(grid)}


@REGISTRY.tool("diagnose_fid",
               "Diagnose an averaged experimental FID: first-point anomaly, saturation plateau, ringing end, "
               "baseline exponential fits and candidate processing recipes (not selected automatically).",
               {"type": "object", "properties": {"source": FID_SOURCE}, "required": ["source"]})
def diagnose_fid(args: dict) -> dict:
    from ..diagnostics import diagnose_fid as run
    exp = _load_fid(args["source"])
    report = run(exp.fid, exp.sampling_rate_hz).to_dict()
    out = new_run_dir("diagnose")
    t = np.arange(exp.points) / exp.sampling_rate_hz
    early = t < max(0.3, 3 * report["ringing_end_s"])
    report["figure_early"] = _plot(out / "early_fid.png", [(t[early], exp.fid[early], "average")], "Time (s)",
                                   "ADC", "Early FID")
    report["source"] = exp.source
    (out / "diagnostics.json").write_text(json.dumps(report, indent=2, default=float), encoding="utf-8")
    report["report_path"] = str(out / "diagnostics.json")
    return report


@REGISTRY.tool("process_fid",
               "Process an averaged FID with an explicit recipe and return its complex spectrum in the given bands.",
               {"type": "object", "properties": {"source": FID_SOURCE, "acquisition": ACQUISITION, "ranges": RANGES,
                                                  "zero_fill": {"type": "integer", "default": 1}},
                "required": ["source", "ranges"]})
def process_fid(args: dict) -> dict:
    from ..solver import ObservedSpectrum
    exp = _load_fid(args["source"])
    acq = _acquisition(exp, args.get("acquisition"))
    obs = ObservedSpectrum.from_fid(exp.fid, acq, [tuple(r) for r in args["ranges"]], args["zero_fill"])
    out = new_run_dir("process")
    np.savez(out / "spectrum.npz", frequency_hz=obs.frequencies_hz, spectrum=obs.values, band=obs.band_index)
    figure = _plot(out / "spectrum.png", [(obs.frequencies_hz, obs.values.real, "real"),
                                          (obs.frequencies_hz, np.abs(obs.values), "magnitude")],
                   "Frequency (Hz)", "Spectrum", "Processed spectrum")
    return {"arrays": str(out / "spectrum.npz"), "figure": figure, "acquisition": acq.to_dict(),
            "points": len(obs.values)}


@REGISTRY.tool("generate_samples",
               "Generate ground-truth samples (spin interpretations with provenance) into gzip JSONL shards.",
               {"type": "object", "properties": {"count": {"type": "integer"}, "seed": {"type": "integer", "default": 0},
                                                  "split": {"type": "string", "default": "train"},
                                                  "problem": {"type": "string"}, "generator": {"type": "string"},
                                                  "output_dir": {"type": "string"}},
                "required": ["count"]}, long_running=True)
def generate_samples(args: dict) -> dict:
    from ..generator import build_default_sampler, write_shards
    from ..generator.sampler import load_generator_config
    spec = ProblemSpec.load(args["problem"]) if args.get("problem") else ProblemSpec()
    config = load_generator_config(args["generator"]) if args.get("generator") else {}
    couplings = config.get("couplings_path")
    if couplings and not Path(couplings).is_absolute() and not Path(couplings).exists():
        config["couplings_path"] = str(Path(args["generator"]).resolve().parent / couplings)
    sampler = build_default_sampler(spec, config)
    out = Path(args.get("output_dir") or new_run_dir("samples"))
    manifest = write_shards(sampler.generate(args["count"], args["seed"], args["split"]), out,
                            manifest_extra={"spec_digest": spec.digest(), "split": args["split"], "seed": args["seed"]})
    return {"directory": str(out), "count": manifest["count"], "acceptance_rate": sampler.acceptance_rate}


@REGISTRY.tool("train_model",
               "Train a candidate model from a run configuration (path or inline dict). Long-running: use submit_job.",
               {"type": "object", "properties": {"run_config": {"description": "Path or inline run configuration."},
                                                  "max_steps": {"type": "integer"}},
                "required": ["run_config"]}, long_running=True)
def train_model(args: dict) -> dict:
    from ..training import TrainingSetup
    setup = TrainingSetup(args["run_config"])
    trainer = setup.build_trainer()
    result = trainer.fit(args.get("max_steps"))
    return {"output_dir": str(trainer.output), "steps": result["steps"], "last_eval": result["last_eval"],
            "device": result["device"], "model_path": str(trainer.output / "model.pt")}


@REGISTRY.tool("propose_candidates",
               "Propose top-k spin interpretations for an averaged FID with a trained model checkpoint. Models "
               "trained on phase-corrected spectra need 'phasing' (phase0_rad and delay_s, the first-order phase "
               "as a time; the crop reference is added automatically).",
               {"type": "object", "properties": {"model_path": {"type": "string"}, "source": FID_SOURCE,
                                                  "acquisition": ACQUISITION, "k": {"type": "integer", "default": 5},
                                                  "device": {"type": "string", "default": "cpu"},
                                                  "phasing": {"type": "object", "properties": {
                                                      "phase0_rad": {"type": "number"},
                                                      "delay_s": {"type": "number"}}}},
                "required": ["model_path", "source"]})
def propose_candidates(args: dict) -> dict:
    from ..evaluation import ModelProposer
    from ..models import load_model
    from ..render import SpectrumGrid
    from ..solver import ObservedSpectrum
    model = load_model(args["model_path"], args["device"])
    exp = _load_fid(args["source"])
    acq = _acquisition(exp, args.get("acquisition"))
    grid = SpectrumGrid.from_spec(model.spec.grid, acq.without_processing())
    obs = ObservedSpectrum(grid.frequencies_hz, np.zeros(len(grid), complex), acq, None, "", {}, exp.fid)
    candidates = ModelProposer(model, grid, model.spec.grid.channels, args["device"],
                               phasing=args.get("phasing")).propose(obs, args["k"])
    return {"candidates": [c.to_dict() for c in candidates], "acquisition": acq.to_dict(),
            "note": "Network proposals; refine and validate before interpreting."}


@REGISTRY.tool("refine_candidates",
               "Refine candidate interpretations against an averaged FID (all branches kept). Held-out FIDs give "
               "frozen-prediction ranking. Settings follow RefineSettings (policy for bounds, ties via "
               "'ties', nuisance terms, background order, continuation, and 'search' for a phase-insensitive "
               "global pattern search that supplies the starts). 'protocol' follows physics.Protocol "
               "(for example field_ut for a residual static field). Long-running for large records.",
               {"type": "object", "properties": {
                   "candidates": {"type": "array", "items": INTERPRETATION}, "source": FID_SOURCE,
                   "held_out": {"type": "array", "items": FID_SOURCE}, "acquisition": ACQUISITION, "ranges": RANGES,
                   "settings": {"type": "object"}, "protocol": {"type": "object"}},
                "required": ["candidates", "source", "ranges"]}, long_running=True)
def refine_candidates_tool(args: dict) -> dict:
    from ..solver import ObservedSpectrum, RefineSettings, refine_candidates
    from ..timing import Timer
    exp = _load_fid(args["source"])
    acq = _acquisition(exp, args.get("acquisition"))
    ranges = [tuple(r) for r in args["ranges"]]
    obs = ObservedSpectrum.from_fid(exp.fid, acq, ranges, label="train")
    held = []
    for k, src in enumerate(args.get("held_out") or []):
        h = _load_fid(src)
        held.append(ObservedSpectrum.from_fid(h.fid, acq, ranges, label=src.get("label", f"held_out_{k}")))
    from ..physics.protocol import SUDDEN_DROP, Protocol
    settings = RefineSettings.from_dict(args.get("settings") or {})
    protocol = Protocol.from_dict(args["protocol"]) if args.get("protocol") else SUDDEN_DROP
    timer = Timer()
    results = refine_candidates([_interp(c) for c in args["candidates"]], obs, held, settings, protocol, timer=timer)
    out = Path(args.get("_job_dir") or new_run_dir("refine"))
    summaries = [r.summary() for r in results]
    (out / "refinement.json").write_text(json.dumps(summaries, indent=2, default=float), encoding="utf-8")
    return {"results": summaries, "acquisition": acq.to_dict(), "timing": timer.report(),
            "report_path": str(out / "refinement.json")}


@REGISTRY.tool("identifiability",
               "Local identifiability of an interpretation's couplings on an averaged FID (singular values and the "
               "weakest parameter combination).",
               {"type": "object", "properties": {"interpretation": INTERPRETATION, "source": FID_SOURCE,
                                                  "acquisition": ACQUISITION, "ranges": RANGES},
                "required": ["interpretation", "source", "ranges"]})
def identifiability_tool(args: dict) -> dict:
    from ..evaluation import local_identifiability
    from ..solver import ObservedSpectrum
    exp = _load_fid(args["source"])
    acq = _acquisition(exp, args.get("acquisition"))
    obs = ObservedSpectrum.from_fid(exp.fid, acq, [tuple(r) for r in args["ranges"]])
    return local_identifiability(_interp(args["interpretation"]), obs).to_dict()


@REGISTRY.tool("submit_job", "Run any tool in a background process; returns a job id to poll with get_job.",
               {"type": "object", "properties": {"tool": {"type": "string"}, "arguments": {"type": "object"}},
                "required": ["tool", "arguments"]})
def submit_job(args: dict) -> dict:
    if args["tool"] not in REGISTRY.tools or args["tool"] in ("submit_job", "get_job", "cancel_job", "list_jobs"):
        raise ValueError("Unknown or non-submittable tool.")
    REGISTRY.tools[args["tool"]].validate(args["arguments"])
    return jobs.submit(args["tool"], args["arguments"])


@REGISTRY.tool("get_job", "Status, progress and (when complete) result of a background job.",
               {"type": "object", "properties": {"job_id": {"type": "string"}}, "required": ["job_id"]})
def get_job(args: dict) -> dict:
    return jobs.get(args["job_id"])


@REGISTRY.tool("cancel_job", "Request cancellation of a background job (terminates it after a grace period).",
               {"type": "object", "properties": {"job_id": {"type": "string"}}, "required": ["job_id"]})
def cancel_job(args: dict) -> dict:
    return jobs.cancel(args["job_id"])


@REGISTRY.tool("list_jobs", "Most recent background jobs.",
               {"type": "object", "properties": {"limit": {"type": "integer", "default": 20}}})
def list_jobs(args: dict) -> dict:
    return jobs.list_jobs(args["limit"])
