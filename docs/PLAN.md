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
| A16 | Global pattern search for refinement starts; optional residual field in the protocol | [x] | `tests/test_search.py` |
| A17 | Phase-corrected real-part input (`grid.phasing`), free component ratios, rendered-weight targets | [x] | `tests/test_phasing.py` |
| A18 | Pre-rendered training shards (`zulf-model prerender`, `data.prerendered`) | [x] | `tests/test_prerender.py` |

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
- [~] CPU verification of the training stack (`configs/run_verify_cpu.json`:
      4-6 spins, up to 2 components, phased real input, 48k pre-rendered
      spectra, CNN set v1 with the group-level head (D23), batch 32, 3000
      steps at about 1 s per step). Validation (128 samples) at steps
      500 / 1000 / 1500 / 2000 / 2500 / 3000: loss 1.06 / 0.86 / 0.74 / 0.67 /
      0.67 / 0.59; structure coverage@10 0.06 / 0.08 / 0.06 / 0.06 / 0.09 /
      0.07; J coverage within 1 Hz 0. All loss terms were still falling. At
      step 2500 the dominant component is right in composition 27 percent and
      group structure 22 percent of top-1 proposals; the main confusions are
      13C versus 15N and the proton count; the model almost always proposes
      one component (component-count accuracy 0.61 equals the one-component
      fraction). Median J errors after structure match at step 3000: 13C-1H
      strong 6.5 Hz, 1H-1H 4.0 Hz, 13C-1H weak 37 Hz. Conclusion: the stack
      learns, the scale is far too small; real runs need a GPU and 1e5+ steps.
      The earlier spin-level head reached structure coverage 0.008 (D23).
- [x] Same data and steps for the CNN+Transformer (`model_cnn_transformer_v1`,
      1.66 M parameters, about 1.2 s per step on CPU). Validation at steps
      500 / 1000 / 1500 / 2000 / 2500 / 3000: structure coverage@1 0.04 /
      0.07 / 0.13 / 0.28 / 0.35 / 0.41; structure coverage@10 0.15 / 0.33 /
      0.31 / 0.42 / 0.40 / 0.45; J coverage (all J within 1 Hz)@1 0 / 0.02 /
      0.05 / 0.09 / 0.13 / 0.17; @10 0.03 / 0.10 / 0.16 / 0.26 / 0.27 / 0.30;
      component-count accuracy 0.79 at step 3000. Median J errors after
      structure match at step 3000: 13C-1H strong 1.2 Hz, 13C-1H weak 1.9 Hz,
      1H-1H 0.2 Hz, 15N-1H weak 0.7 Hz, but 15N-1H strong 115 Hz (relative
      sign of 1J(15N,1H), which is negative, versus 1J(13C,1H), positive).
      At the step-1000 checkpoint the soft-target entropy floor of the token
      loss is 1.07 nats (KL 1.74 above it), structure-token accuracy 0.73,
      J-token accuracy within one bin 0.12. Decision: the sequence model is
      the main model; the set model stays as a baseline.
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
- [~] Global pattern search (D18) places the strongest methine line at the
      observed 133.4 Hz, where every local fit from the earlier start failed;
      the joint search-then-refine run on the real average is pending. With
      JHH fixed at 6.4 Hz the methine multiplet intensities still do not match
      (pattern cost about 0.2 on 127-141 Hz); a residual field (D19) improves
      it only marginally at its bound. Paused at the user's request in favour
      of training work; next questions are the pulse sequence and shield field.
      Automatic phasing of the real average (crop 100 ms, SG801) gives
      phase0 -64 deg and delay -6.1 ms relative to the first recorded sample
      (misfit 0.065 over 14 peaks): spin evolution starts about 6 ms after the
      record starts. Across 110-150 Hz this is about 1.5 rad of first-order
      phase that the refinement (shared phase, no delay) could not follow; the
      next refinement should fix the delay at this value or fit it
      (`policy.fit_phase_delay`).
- [ ] Recovery tests on random generated systems from perturbed starts.

## Phase 4: architecture comparison (plan weeks 6-9)

- [ ] Train CNN+Transformer; compare (a) graph search, (b) CNN set, (c)
      CNN+Transformer under identical budgets.

## Model optimization roadmap (from the second survey, docs/REFERENCES.md)

Ordered by expected gain per effort; each item is measured on the same
verification data before adoption.

- [ ] Rerank the top-k beam by short solver fits (residual), adding mutation
      neighbours (13C/15N swap, group size +/-1) for the observed confusions.
- [x] Sign variants in refinement (D25, `RefineSettings.sign_variants`); to be
      enabled for every refinement against experimental spectra.
- [~] Noise-free training control: same systems and steps as the noisy
      CNN+Transformer run, to separate missing information from model limits.
- [~] Faster rendering: total-M blocks inside each collective sector (done,
      exact, transition lists 42 ms -> 15 ms per 8-spin sample); next a
      batched GPU renderer with blocks of at most 32.
- [ ] J decoding: local expectation over bin probabilities, a Wasserstein
      term over bin centres, then coarse-to-fine J tokens.
- [ ] Weak components: explain-away second pass (fit, subtract, re-run),
      noise-component augmentation and end-of-sequence down-weighting.
- [ ] Encoder: peak-token stream with high-resolution frequency features and
      pairwise frequency-difference attention; diverse beam search by
      structure prefix.

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

CPU container, after the operator and NUFFT speedups (commit "Phased input,
free ratios, pre-rendered shards"):

| Quantity | Value |
| --- | --- |
| Transition lists per 8-spin CHN sample (uncached) | 110 ms before, 55 ms after (pair operators built as Kronecker products) |
| NUFFT, 150 lines x 65516 samples | 43 ms before, 13 ms after (fast FFT length; relative error 1.6e-11) |
| Rendered sample, experimental route (instrument nuisances, SG, 65516-point record, phased) | 338 ms before, 160 ms after |
| Rendered sample, verification problem (4-6 spins, pure, 8192 points, 3268-point grid) | about 20 ms |
| Pre-rendered shard size (float16, one channel, 6.5k-point grid) | about 18 KB per item |
| Training step, CNN set v1 (1.47 M parameters), batch 32, 3268-point grid, CPU 4 threads | about 0.55 s |
