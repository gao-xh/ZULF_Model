# Architecture

ZULF_Model learns to propose candidate spin interpretations of a zero-field
NMR spectrum and refines each candidate with a physical forward model. The
network never has the final word: every candidate is re-simulated, refined and
checked against held-out acquisitions.

```
                         offline training line
  +-----------+   +---------+   +----------+   +---------+   +---------+
  | generator |-->| physics |-->| renderer |-->|  codec  |-->| trainer |
  | (graphs,  |   | (trans- |   | (process |   | (targets|   | (models)|
  |  J rules) |   |  itions)|   |  + noise)|   |  tokens)|   |         |
  +-----------+   +---------+   +----------+   +---------+   +---------+
        ^                                                         |
        | focused resampling                              weights |
  +-----------+                                                   v
  | finetune  |<-- failures --+  +-----------+   +--------+  +----------+
  +-----------+               +--| evaluation|<--| solver |<-| proposer |
                                 +-----------+   +--------+  +----------+
                         online inference line (experimental spectrum)
```

## Two packages (D29)

- `zulf_core` (NumPy, SciPy): nuclei, spin systems, zero-field physics, the
  acquisition operator and exact rendering, grids and phasing, the refinement
  solver (parameterization, forward model, global search, sign variants),
  matching and identifiability, FID I/O and diagnostics. Everything needed to
  simulate and to fit experimental spectra.
- `zulf_model` (adds PyTorch): problem spec, generator, codec, training-data
  synthesis (perturbations, sample pipeline, features), models, training,
  fine-tuning, proposers and benchmarks, the agent tool layer and CLI.

`zulf_model` imports `zulf_core`; `zulf_core` never imports `zulf_model` or
torch (`tests/test_spinsystem.py` checks this), so training data and
experimental fits always share one forward model. Short module names in these
docs (`physics.transitions`, `render.phasing`, `solver.search`) refer to
`zulf_core`; `render.pipeline`, `render.perturb` and `render.features` are in
`zulf_model.render`, which also re-exports `zulf_core.render`.

## Package map

| Module | Responsibility | Depends on |
| --- | --- | --- |
| `zulf_core.nuclei` | Extensible nucleus registry: symbol, spin quantum number, gamma. | - |
| `zulf_core.spinsystem` | `SpinSystem`, `Component`, `Interpretation`; validation, equivalence groups, permutations, matching, JSON. | nuclei |
| `zulf_core.physics` | Spin operators, collective-spin sectors, zero-field Hamiltonian, transition lists (frequency, complex amplitude), observation protocol. | spinsystem |
| `zulf_core.render` | Acquisition and the single preprocessing operator; exact analytic and NUFFT rendering of transition lists through that operator; grids; phasing. | physics |
| `zulf_model.render` | Training-data synthesis: perturbations, sample pipeline, model input features (re-exports `zulf_core.render`). | zulf_core |
| `zulf_model.generator` | Molecule-like heavy-atom graphs, rule-based J assignment, isotopologue enumeration, random-J mode, sample specification, family-based splits, shard storage. | spinsystem, physics, render |
| `zulf_model.codec` | Canonical ordering, J binning, token vocabulary and grammar, fixed-size set targets with masks. | spinsystem |
| `zulf_model.models` | `CandidateModel` interface, CNN spectrum encoder with absolute-frequency features, CNN set-prediction baseline over equivalence groups, CNN+Transformer encoder-decoder with constrained beam search. | codec (torch) |
| `zulf_model.training` | Torch datasets over generator or shards, pre-rendered feature shards, permutation-aware losses, metrics, `Trainer` with checkpoints, curriculum and logging. | models, generator |
| `zulf_core.solver` | Observed-spectrum container, general parameterization (free, fixed, tied J; per-component or per-family rates), variable-projection forward model, phase-insensitive global pattern search for starts, bounded multistart refinement, batch refinement, frozen held-out prediction. | physics, render |
| `zulf_core.evaluation` | Matching of interpretations (per-block sign freedom), local identifiability, solver basin measurement. | solver |
| `zulf_model.evaluation` | Candidate proposers (model, random multistart, graph search) and the benchmark runner. | zulf_core, generator |
| `zulf_model.finetune` | Failure classification, focused resampling that respects frozen test families, active-learning loop around the trainer. | evaluation, training |
| `zulf_core.io` | Loading averaged FIDs (npy, legacy NMRduino DAT decoding, INI parsing, group averaging). | - |
| `zulf_core.diagnostics` | Per-dataset FID diagnostics (first point, saturation plateau, ringing end, baseline fits) and candidate processing recipes. | - |
| `zulf_model.agent` | Tool registry shared by the JSON CLI, the MCP server and exported Anthropic/OpenAI tool schemas; background jobs. | all |
| `zulf_model.device`, `zulf_core.timing` | CUDA/MPS/CPU policy; built-in timers. | - |
| `zulf_model.cli` | Command line entry points (thin layer over the agent registry). | agent |

Only `models`, `training`, `finetune` and parts of `evaluation` import torch.
Physics, rendering, generation and refinement run on NumPy/SciPy alone.

## Core data objects

- `SpinSystem(isotopes, couplings_hz, groups=None)`: isotopes is a tuple of
  registry symbols; couplings is a real symmetric matrix with zero diagonal in
  Hz; groups are optional zero-based magnetic-equivalence groups. The spin
  count is `len(isotopes)`; nothing assumes eight.
- `Component(system, contribution, label, metadata)`: one isotopologue or spin
  fragment with a nonnegative relative contribution (natural abundance times
  count of equivalent positions, before experimental response).
- `Interpretation(components, score, metadata)`: an explanation of a whole
  spectrum. Candidate lists are lists of interpretations.
- `TransitionList(frequencies_hz, amplitudes, dc, metadata)`: positive
  frequencies with complex amplitudes `a` such that the detected real signal is
  `sum Re(a exp(2 pi i f t))` plus `dc`.
- `Acquisition(sampling_rate_hz, points, start_sample, stop_sample, sg_window,
  sg_order, remove_mean, time_origin_s)`: the complete experimental processing
  recipe. `render.preprocess.process_record` is the only implementation.
- `SpectrumGrid(frequencies_hz)`: frequencies at which spectra are evaluated.
  Native bins or zero-filled bins are both exact evaluations of the processed
  finite record, never interpolation.
- `ObservedSpectrum(grid, values, band_index, acquisition, label)`: complex
  spectrum values used by the solver.

## Interfaces meant for extension

| Interface | Location | Extend by |
| --- | --- | --- |
| Nuclei | `nuclei.NucleusRegistry.register` | Adding 19F, 31P or quadrupolar nuclei; the physics layer handles any spin quantum number. |
| Observation protocol | `physics.protocol.Protocol` | New preparation or detection weights, ideal pulses. Finite pulses and fields are planned extensions. |
| Renderer backend | `render.renderer.Renderer` | A torch or GPU backend with the same `render(transitions, rates, grid)` contract. |
| Perturbations | `render.perturb.Perturbation` | New nuisance processes; each takes an RNG and returns time-domain additions or amplitude transforms. |
| Coupling rules | `generator.couplings.CouplingRules` | New coupling families; ranges live in configuration. |
| Sample sources | `generator.sampler.SampleSource` | Random J, labeled, natural-abundance or focused sources, mixed by weights. |
| Candidate model | `models.base.CandidateModel` | New architectures implementing `loss` and `propose`. |
| Candidate proposer | `evaluation.proposers.CandidateProposer` | Any strategy that turns an observed spectrum into interpretations. |
| Solver objective | `solver.forward.MixtureForward` | Magnitude or windowed objectives, background terms. |
| Global search | `solver.search.PatternObjective`, `global_search` | Other phase-insensitive objectives or search strategies that return distinct starting vectors for `refine(initial_points=...)`. |
| Evolution field | `physics.protocol.Protocol.field_ut` | Residual static field during evolution; exact in the sector decomposition. |
| Input phasing | `spec.GridSpec.phasing`, `render.phasing` | "none" (random phase, real+imag) or "corrected" (manual-style 0/1-order phasing, real part). |
| Training data source | `training.prerender.PrerenderedDataset` | Live rendering (default) or pre-rendered shards; any iterable of items in the `make_item` format. |

## Data flow for one training sample

1. A `SampleSource` draws a molecule graph (or a random spin system), assigns
   J by rules, enumerates isotopologues and returns an `Interpretation` plus a
   family identifier used for splitting.
2. `physics.transitions` computes each component's transition list under the
   frozen protocol, exploiting magnetic-equivalence sectors.
3. `render` evaluates every component through the acquisition operator with
   random linewidths, gains and phases, then adds preprocessed noise, drift and
   interference generated in the time domain.
4. `render.features` turns the complex spectrum into model channels with a
   recorded normalization scale.
5. `codec` converts the interpretation into set targets and a token sequence.

## Inference flow

1. Load an averaged FID, process it with the frozen acquisition recipe and
   evaluate it on the model grid.
2. A proposer returns the top-k interpretations.
3. The solver refines each candidate independently, keeping every branch.
4. Frozen parameters are evaluated on held-out acquisition groups. Ranking uses
   held-out prediction, not training residual.

## AI agent interface

`zulf_model.agent.tools` registers typed tools (JSON schema in, JSON out):
describe_project, simulate_transitions, render_spectrum, diagnose_fid,
process_fid, generate_samples, train_model, propose_candidates,
refine_candidates, identifiability, submit_job, get_job, cancel_job, list_jobs.
The same registry is served over MCP (`python -m zulf_model.agent.mcp_server`),
exported for the Anthropic and OpenAI APIs (`zulf-model tools export`), and
called by the CLI. Long-running tools run as background jobs whose status
survives client reconnects. `skills/zulf-model/SKILL.md` tells an agent the
workflow and what it may and may not conclude.

## Reuse of earlier projects

The physics and rendering layers are rewritten from ideas in
`ZULF_Analysis_Tools` (`simulation.py` sectors, `jfit.ProcessedSpectrum`
templates, `decay.fit_modes` budgets and variable projection) and
`ZULF_NMR_Suite` (`TwoD_simulation.py` sequences). They are not imported, so
this package has no dependency on either repository. Known legacy errors (the
15N gamma in `TwoD_simulation.py`) are not carried over.
