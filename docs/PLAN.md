# Implementation plan and ledger

Status legend: `[x]` done and tested, `[~]` implemented with open items,
`[ ]` not started. Each completed item names the test module that covers it.
Update this file in the same commit as the work.

Source plan: "ZULF spectrum-to-J self-trained model" (2026-09-24). Phase
numbers below follow that plan.

## Phase A: foundation (this repository's first milestone)

The goal is a complete, tested skeleton of every component, with interfaces
fixed, so later phases are experiments rather than restructuring.

| Step | Deliverable | Status | Tests |
| --- | --- | --- | --- |
| A1 | Docs: architecture, conventions, decisions, agent rules | [x] | `scripts/check_ascii.py` |
| A2 | `ProblemSpec` single source of problem dimensions | [x] | `tests/test_spec.py` |
| A3 | Nucleus registry, `SpinSystem`, `Component`, `Interpretation`, equivalence, permutation matching, JSON | [x] | `tests/test_spinsystem.py` |
| A4 | Physics: operators, collective sectors, transitions, protocol | [x] | `tests/test_physics.py` |
| A5 | Rendering: optional processing operator, analytic + NUFFT time-domain renderers, continuous Voigt route, perturbations, features, sample pipeline, timers | [x] | `tests/test_render.py` |
| A6 | Generator: graphs, coupling rules, isotopologues, random J, sources, splits, storage | [x] | `tests/test_generator.py` |
| A7 | Codec: canonical order, J bins, tokens, grammar, set targets | [x] | `tests/test_codec.py` |
| A8 | Models: interface, CNN encoder, set baseline, CNN+Transformer, beam search | [x] | `tests/test_models.py` |
| A9 | Training: datasets, losses, metrics, trainer, curriculum, checkpoints | [x] | `tests/test_training.py` |
| A10 | Solver: observed spectra, parameterization, forward, refine, batch, held-out, matched linewidth continuation | [x] | `tests/test_solver.py` |
| A11 | Evaluation: proposers, identifiability, basin, benchmark | [x] | `tests/test_evaluation.py` |
| A12 | Fine-tuning: failure classification, focused sampling, loop | [x] | `tests/test_finetune.py` |
| A13 | CLI, configs, end-to-end smoke pipeline, throughput script | [x] | `tests/test_agent.py`, `scripts/smoke_pipeline.py` |
| A14 | AI agent layer: tool registry, MCP server, Anthropic/OpenAI export, jobs, skill guide | [x] | `tests/test_agent.py` |
| A15 | Per-dataset diagnostics and solver nuisance terms (exponential, damped sinusoid, template) | [x] | `tests/test_solver.py` |

## Phase 0: basis (plan week 1)

- [ ] Freeze protocol after Q1 (`configs/protocol_v1.json`).
- [ ] Freeze acquisition and grid after Q2 (`configs/acquisition_v1.json`).
- [ ] Throughput table for 8-spin CHN components, CPU and GPU, with and
      without equivalence sectors; 9-12 spin interface check
      (`scripts/throughput.py`, results recorded here).
- [x] Cross-check transitions against `ZULF_Analysis_Tools.jfit.transitions`:
      both isopropylamine isotopologues agree (16 and 108 transitions,
      frequencies and normalized weights to 1e-12; `tests/test_crosscheck.py`,
      skipped without that repository).

## Phase 1: data generator (plan weeks 2-3)

- [ ] Verify coupling ranges against literature (Q3); record sources in
      `configs/couplings_v1.json`.
- [ ] Generate and freeze `gen-v1` train/validation/test shards with family
      splits and an unseen-topology test set; record digests here.

## Phase 2: M0 feasibility (plan weeks 4-5)

- [ ] Measure solver basin of attraction per coupling type
      (`evaluation.basin`), set J bin widths from it.
- [ ] Identifiability pairs (`evaluation.identifiability`).
- [ ] Freeze pass criteria before training (coverage@k, J tolerance,
      held-out threshold, budget).
- [ ] Train CNN set baseline; benchmark against random multistart and graph
      search at equal simulation budget.

## Phase 3: general solver (plan weeks 4-6, parallel)

- [~] Reproduce the isopropylamine staged-fit result with `solver.refine`
      and explicit ties. First run on the uploaded 10000-scan average
      (recipe crop 100 ms + SG801 + band quadratic background, bands 110-150
      and 230-275 Hz, starting J from the earlier tool): relative residual
      0.063, methine rate at its 20/s bound (flag `search_boundary`), vicinal
      H-H 8.2 Hz. Not accepted: needs recipe comparison (explicit nuisance
      baseline vs SG), held-out group averages and comparison with the
      earlier tool's fitted values.
- [ ] Recovery tests on random generated systems from perturbed starts.

## Phase 4: architecture comparison (plan weeks 6-9)

- [ ] Train CNN+Transformer; compare (a) graph search, (b) CNN set, (c)
      CNN+Transformer under identical budgets.

## Phase 5: closed loop (continuous)

- [ ] Dev-set failure mining and focused resampling rounds; each round records
      generator version, data digest and metrics.

## Phase 6: blind test (plan weeks 10-11)

- [ ] Freeze model, generator and protocol; run on molecules from Q9.

## Phase 7 (optional): J to structure

- [ ] Not started. The generator already records graphs for paired data.

## Experimental data inventory (user Google Drive, "Metabolites/Low Gamma")

"Good Scan Data" and "Raw Data" folders hold zipped NMRduino runs (0.2-2 GB each):
isopropylamine, tert-butylamine, diethylamine (neat and 50% in d-benzene),
dipropylamine (several dilutions), triethylamine, ethylenediamine,
N-ethylmethylamine, N,N-dimethylethylenediamine, pyridine series including 1%
15N-enriched, L-alanine, L-lactic acid, WC in CDCl3, Nsime3, ethanol,
methyl dimethylphosphonate; a "py_dil" folder adds 10 and 33 mol% pyridine.
The Drive connector can list these files but returns content inline, which does
not scale to zips of this size, and direct download hosts are blocked by the
environment network policy. Averaged FIDs (npy plus ini, a few MB each) can be
uploaded instead.
Isopropylamine is the development molecule (never a blind test); blind-test
candidates are chosen from the rest before models are frozen (Q9).

Isopropylamine average (10000 scans, fs 4000 Hz, 65516 decoded points,
16.4 s): the ADC plateau (about 28848) is present in every acquisition for the
first about 3 ms, followed by about 500 Hz ringing decaying by about 25 ms and a
slow baseline of about -1.5e4 ADC decaying over seconds. Three exponentials
from 30 ms leave 73 ADC RMS versus 116 for SG801 (single test, not yet a
chosen recipe). Processing parameters are chosen per dataset.

## Measurements

CPU container (4 cores, no GPU), commit after "NUFFT renderer":

| Quantity | Value |
| --- | --- |
| Transition list, 8-spin CHN component, median over generated samples | about 140 transitions |
| NUFFT synthesis, 2000 transitions x 16384 samples | 8.5 ms, relative error 6e-12 |
| Rendered training sample (mixture, noise, randomized processing, 16384-point record, 6537-point grid) | about 60 ms, of which about 39 ms transition lists (cached per system) |
| Continuous (infinite-record) route, same grid | about 0.5-0.9 s per sample (exact reference route, not for bulk training) |
| Sampler acceptance (8-spin natural isotopologues by rejection) | about 0.27 |
