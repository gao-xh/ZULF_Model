# ZULF_Model

Simulation-trained candidate generation and general J refinement for zero- to
ultralow-field (ZULF) NMR spectra.

A model reads a ZULF spectrum and proposes several candidate interpretations
(isotopologue components, magnetic-equivalence groups, J matrices). Each
candidate is re-simulated with an exact forward model, refined against the data
and ranked by frozen prediction on held-out acquisitions. Training needs no
labelled experimental database; experiments remain the final check.

## Layout

| Path | Contents |
| --- | --- |
| `zulf_core/` | NumPy/SciPy core: spin systems, physics, rendering, phasing, solver, matching, FID I/O, diagnostics |
| `zulf_model/spec.py` | `ProblemSpec`: every problem dimension (spin counts, nuclei, components, J bins, grid) |
| `zulf_core/physics` | Zero-field Hamiltonians, collective-spin and total-M sectors, transition lists, protocol |
| `zulf_core/render` | Optional processing operator, exact analytic and NUFFT renderers, Voigt lines, grids, phasing |
| `zulf_model/render` | Training-data synthesis: perturbations, sample pipeline, features |
| `zulf_model/generator` | Molecule-like graphs, J rules, isotopologues, random-J, splits, shards |
| `zulf_model/codec.py` | Token grammar and set targets |
| `zulf_model/models` | CNN set baseline and CNN+Transformer with constrained beam search |
| `zulf_model/training` | Datasets, trainer (CUDA / Apple MPS / CPU), metrics, curriculum |
| `zulf_core/solver` | General refinement: ties, nuisance terms, backgrounds, guarded continuation, global search, sign variants, held-out ranking |
| `zulf_core/evaluation` | Matching of interpretations, identifiability, solver basin |
| `zulf_model/evaluation` | Proposers and benchmarks |
| `zulf_model/finetune` | Failure mining and active-learning rounds |
| `zulf_core/diagnostics.py` | Per-dataset FID diagnostics and candidate processing recipes |
| `zulf_model/agent` | Tool registry for AI agents: JSON CLI, MCP server, Anthropic/OpenAI tool export, background jobs |
| `configs/` | Problem, generator, couplings, processing, perturbation, protocol, model and run configurations |
| `skills/zulf-model` | Agent skill guide (Claude skill format and OpenAI agent card) |
| `skills/zulf-blind-analysis` and category skills | Experimental workflow skills: FID processing, spectrum interpretation, phasing, hypothesis refinement; case notes and derivation paths |
| `zulf_studio/` | ZULF Studio: real-time simulator (structure, couplings, static field, line width), fitting with a progress slider, log, terminal, Python console (PySide6) and a JSON API for AI agents |
| `docs/` | Architecture, conventions, decisions, plan ledger, references |

## Install

On the owner's machine: the conda environment `zulf` (docs/ENVIRONMENT.md). Elsewhere:

```bash
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[dev]"                            # numpy, scipy, torch, h5py, matplotlib
pip install "mcp>=1.12,<2"                         # optional, for the MCP server
```

PyTorch picks CUDA, Apple MPS or CPU automatically (`train.device` = `auto`).
Simulation always runs in float64 on the CPU.

## Use

```bash
python -m unittest discover -s tests              # full test suite
python scripts/analyze_sample.py FID.npy --id ID  # new sample: overview, blind search, report (runs/blind/ID)
python scripts/analyze_sample.py FID.npy --id ID --structure '{"motif": "ethyl", "one_bond": {"C1": 131, "C2": 125}}'
ZULF_DATA_DIR=... python scripts/regression_confirmed.py --mode both   # confirmed-sample regression
python scripts/make_series_entry.py --fid FID.npy --id NAME --out runs/series/NAME   # FID -> processed spectrum
python scripts/fit_joint_series.py --series runs/series/NAME/series.json ...        # complex known-structure fit (docs/WORKFLOW.md)
python scripts/band_diagnosis.py <fit options> --fit runs/processed/NAME/fit.json    # which band is bad, which parameters act on it
python scripts/smoke_pipeline.py                  # generate -> train -> propose -> refine
python scripts/run_studio.py --series runs/series/NAME/series.json --fit runs/processed/NAME   # ZULF Studio (docs/STUDIO.md): sliders, fit, figures, AI assistant; JSON API on :8766
zulf-model tools list                             # agent tools
zulf-model diagnose AVERAGE.npy 0.ini             # per-dataset processing diagnostics
zulf-model train configs/run_cnn_set_v1.json      # training run (live rendering)
zulf-model prerender configs/run_cnn_set_phased_v1.json runs/shards/phased --count 200000 --workers 8
                                                   # pre-render shards, then set "data": {"prerendered": ...}
zulf-model tools export --format anthropic        # tool schemas for the Claude API
python -m zulf_model.agent.mcp_server             # MCP server for Claude Code / Codex
```

Outputs go to `$ZULF_MODEL_WORKSPACE` (default `runs/`, ignored by git).
Experimental inputs are read-only and never committed.

## Data

Experimental data is not in git. On the local machine of the owner it is stored under
`~/research/zulf/data/` (layout, rules and inventory with sha256 in its `README.md`). Every file sits in a
measurement folder (`<YYYY-MM-DD>-<compound>`, or `<compound>-<sample id>` when the date is unknown), created
first under each area:

| Folder | Contents |
| --- | --- |
| `original/<measurement>/` | Instrument output of one run, every scan (NMRduino `<n>.dat`, `<n>.ini`), as given, read-only |
| `raw/<measurement>/` | Other received data, original names: averaged FID `<id>-<name>.npy`, `selected-scans/` (legacy ScanSelector), `spectra/` (other groups or software) |
| `processed/<measurement>/` | Data products of our code (scan-selection records and their averages) |

`ZULF_DATA_DIR` points at `raw/` and is passed per command, e.g.
`ZULF_DATA_DIR=~/research/zulf/data/raw python scripts/regression_confirmed.py ...`; the scripts find each
FID by file name in any subfolder (`zulf_core.io.find_data_file`, exactly one match required). Cloud sessions
have no copy of the data (docs/HANDOFF.md, "Data").

## Status

Start with `docs/HANDOFF.md` (state as of 2026-10-05: data, current results, how to fit a known structure, open
work), then `docs/WORKFLOW.md` (map of every workflow) and `docs/PLAN.md`. Results from the solver are conditional numerical
candidates, reported with their flags, residuals and held-out scores.
