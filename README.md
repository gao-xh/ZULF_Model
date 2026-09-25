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
| `zulf_model/spec.py` | `ProblemSpec`: every problem dimension (spin counts, nuclei, components, J bins, grid) |
| `zulf_model/physics` | Zero-field Hamiltonians, collective-spin sectors, transition lists, protocol |
| `zulf_model/render` | Optional processing operator, exact analytic and NUFFT renderers, Voigt lines, perturbations |
| `zulf_model/generator` | Molecule-like graphs, J rules, isotopologues, random-J, splits, shards |
| `zulf_model/codec.py` | Token grammar and set targets |
| `zulf_model/models` | CNN set baseline and CNN+Transformer with constrained beam search |
| `zulf_model/training` | Datasets, trainer (CUDA / Apple MPS / CPU), metrics, curriculum |
| `zulf_model/solver` | General refinement: ties, nuisance terms, backgrounds, matched continuation, held-out ranking |
| `zulf_model/evaluation` | Proposers, identifiability, solver basin, benchmarks |
| `zulf_model/finetune` | Failure mining and active-learning rounds |
| `zulf_model/diagnostics.py` | Per-dataset FID diagnostics and candidate processing recipes |
| `zulf_model/agent` | Tool registry for AI agents: JSON CLI, MCP server, Anthropic/OpenAI tool export, background jobs |
| `configs/` | Problem, generator, couplings, processing, perturbation, protocol, model and run configurations |
| `skills/zulf-model` | Agent skill guide (Claude skill format and OpenAI agent card) |
| `docs/` | Architecture, conventions, decisions, plan ledger, references |

## Install

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
python scripts/smoke_pipeline.py                  # generate -> train -> propose -> refine
zulf-model tools list                             # agent tools
zulf-model diagnose AVERAGE.npy 0.ini             # per-dataset processing diagnostics
zulf-model train configs/run_cnn_set_v1.json      # training run
zulf-model tools export --format anthropic        # tool schemas for the Claude API
python -m zulf_model.agent.mcp_server             # MCP server for Claude Code / Codex
```

Outputs go to `$ZULF_MODEL_WORKSPACE` (default `runs/`, ignored by git).
Experimental inputs are read-only and never committed.

## Status

See `docs/PLAN.md`. Results from the solver are conditional numerical
candidates, reported with their flags, residuals and held-out scores.
