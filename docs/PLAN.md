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
- [x] Signal-focused weighting with model cores from the transition-list
      envelope (D32; d5543fc, 196e14f, 9e46dbd, 6595bad) and fixed amplitude
      ratios (D33; 2a16dfe). Measured on e3d282da / b683220d, see
      docs/ANALYSIS_LOG.md.

## Phase 3b: programmatic hypotheses (D34-D37, package `zulf_hypothesis`)

- [x] Stage 1: fragments with symmetry, isotopologue builder with automatic
      ties/ratios/omissions, checks registry, extension-move interface with
      `AddCoupledProton`, reference couplings (L-alanine, lactic acid,
      measured). Replay on b683220d reproduces the hand-built fits.
- [x] Stage 2 (D35): band inventory, group candidates, fragment enumeration,
      `propose_hypotheses`; neural models as hint providers with an insight
      report. Both confirmed samples rank the correct CH3-CH first; the
      set-model checkpoint agrees with the data in one band per sample.
- [x] Minor isotopologues (D36): `Labeling` (natural default, enriched
      uniform or site-specific), exact label-set weights, minor sets tied to a
      parent via `amplitude_map`, SNR gate with omission reasons.
- [ ] Better hints: train or fine-tune a group-level head on the
      hypothesis vocabulary (group type, 1J) so model output maps directly
      onto group hints; calibrate hint confidence on the synthetic benchmark.
- [x] Stage 3 (D38): `search_hypotheses` (parallel refinement in free and
      fixed variants, common yardstick and BIC, checks, triggered extension
      moves with acceptance by BIC, knowledge matching, ranked table and
      log).
- [x] Motif library scan (D39) and recognition benchmark
      (`python -m zulf_hypothesis.benchmark`): motif rank 1 in 11 / 14 cases,
      including slow-exchange amines (13th: ethylenediamine, fde3fbb2;
      14th: triethylamine, 4322bdfc, rank 2).
- [ ] More motifs as structures are confirmed (CH2-CH, rings with
      substituents, carbonyl-bridged groups); a stage-3 benchmark column.
- [x] Lessons from 4322bdfc (triethylamine), D42: motif screen (short
      refinement of the first 8 scan proposals before the cutoff), move
      change_proton_count (CH <-> CH2 <-> CH3), check collapsed_component
      and demotion of such fits in the ranking. Rerun on the five confirmed
      real samples: see ANALYSIS_LOG.
- [ ] Automatic fitting without hand-proposed variants (e66a4b08,
      N-ethylmethylamine; user: "the fit algorithm should do these tests
      itself"). Steps:
  - [x] A1 variant "ratios" (fixed natural ratios, free rates) next to free
        and fixed, in the search and every known-structure fit.
  - [x] A2 exchange variants: for a fragment with N-H / O-H groups, fit
        fast exchange (groups dropped) and slow exchange (kept)
        automatically.
  - [x] A3 round-0 starts: optional multi-start plus global pattern search
        for the round-0 fits (on by default for a known structure).
  - [x] A4 `fit_structure(fragment, observed)`: a thin entry point over
        `search_hypotheses` (the solver is unchanged) running A1-A3,
        ranking, checks and a phased result figure.
  - [x] B change_proton_count rescales 1J so the main X-Hn line stays in
        place (XH at J, XH2 at 1.5 J, XH3 at J); e66a4b08 turned a CH at
        204 Hz into a CH2 at 204 Hz instead of 136 Hz.
  - [x] C automatic phasing as a library function, maintained: per-line
        complex line fits on resolved lines, then a robust fit of line phase
        against frequency (modulo pi) for phase0 and delay; each spectrum
        phased on its own. Validation: synthetic test, noise-free rendered
        e66a4b08 model, and the six confirmed samples against the phase of
        their known-structure complex fits. (The minimum-entropy estimator
        failed the synthetic test and is not in the library.)
  - [x] D rerun e66a4b08 with A4 (timed, both routes: complex fit and
        phased-first), tests, benchmark, decision entry, skills.
  - [x] E line shape: Gaussian width (Voigt) as a triggered variant; width
        range from the data (no instrument prior known).
  - [x] F remote couplings: a triggered move that frees the couplings built
        as 0 (4J, 5J) between labelled sites and proton groups.
  - [x] G intermediate exchange (D45, 2026-09-29): Liouville-space model
        with a fitted exchange rate per exchangeable group
        (zulf_core.physics.exchange, analytic derivatives, optional torch/GPU
        backend); fit_staged.py --exchange. Earlier: both limits were fitted
        on e66a4b08 and neither explained the misfit. Applied first to
        pyridinium (Blake 6 M HCl).
  - [x] H automatic report per fit: ranked table, phased fit figure, J
        matrix file, findings.
  - [x] R regression script over the confirmed samples (blind search and
        known-structure fit; skeleton found, chi2, time, 1J, phase), run
        after each change.
  - Result (ANALYSIS_LOG, final regression): blind 5/6 at rank 1, N-ethylmethylamine rank 2;
    C done as fit-phased route (model-free phasing unreliable on dense multiplets, D43).
    Open: known-structure budget for 11-spin forms (10-22 min), blind two extension rounds,
    physics double diagonalisation (compute_transitions + transition_derivatives), delay prior values.
  - Budget agreed with the user: blind sample <= 30 min, known structure
    <= 10 min (4 cores). Variants beyond the first round are triggered by
    findings, not run for every hypothesis. Autonomous 6 h session started
    2026-09-27 (UTC, see ANALYSIS_LOG); new problems found on the way are
    fixed and logged.
- [ ] More moves: add 15N site, change equivalence (CH3 <-> CH(CH3)2), add a
      remote proton group, split a group; more templates as cases arrive.
- [ ] Validation: replay the four blind datasets end to end; synthetic
      benchmark from the generator (fraction of correct fragments within a
      budget), compared with the neural proposers.

## Phase 3c: experimental processing (D44, package `zulf_processing`)

- [x] Package split: raw, diagnostics, plan, phase, dataset; tests
      (tests/test_processing.py); validation script
      (scripts/validate_processing.py) against the known-structure fits.
- [x] Delay per dataset from the switching edge of the raw FID.
- [x] Diagnostics speed: 546 s -> 0.1 s per dataset (exponential fits
      skipped in diagnose_raw).
- [x] Default phase: instrument calibration + per-dataset edge
      (`calibrated_phase`, config "phase_calibration"; within 14 deg
      leave-one-out on 5 samples).
- [ ] Zero-order phase per dataset from the data alone: better criteria
      (validation: entropy / lines err 2-86 deg, worse than the
      calibration); ideas: a phase reference recorded with every acquisition
      (test signal), the ringing after the edge, joint use of several
      processing windows, criteria restricted to +-20 deg of the calibration.
- [x] Window length per dataset (user: a too-short window loses spectral
      information): `signal_extent` + plan window_mode "signal_extent"
      (commit 160c33d; default still "fixed"). Known fits (ANALYSIS_LOG):
      lactic acid and triethylamine better (delay at the edge, phased-route
      chi2 halved), N-ethylmethylamine mixed; cost 4-6x.
- [ ] Renderer / SG cost for long records (needed before the long window
      becomes the default).
- [ ] Soft crop window (flat middle, Gaussian roll-off at both ends) and a
      weighted SG kernel of the same shape (user proposal; design in
      docs/HANDOFF.md), inside zulf_core's single operator, default off.
- [ ] Processing parameters per dataset from the data: SG window from the
      lowest band / baseline rates, apodization vs stop taper; each choice
      validated on the confirmed samples (J stability, uncertainty, c_hat).
- [ ] Instrument-line detection per dataset instead of a fixed list.
- [x] Wire `process_dataset` into the regression:
      `regression_confirmed.py --processing dataset [--window-mode ...]`.
- [~] Scan averaging with per-scan screening (`scripts/average_scans.py`, tests/test_diagnostics.py
      AverageScansTests; commit of the acetonitrile analysis): robust-z deviation metric, even / odd half
      averages, a record (scans.json). Open: the full rewrite of the legacy ScanSelector (signal-band metric,
      leave-one-out reference, held-out check of the selection).
- [ ] Phase calibration of the 2 kHz sequence (standard_zf_2000Hz_no_dead.seq); the config holds only the
      4 kHz one (acetonitrile 2026-09-28: calibrated spectrum far from absorptive).

## Phase 3d: processed spectra from other groups (Blake pyridine series)

- [x] Known-structure fits of processed real spectra (fit_processed_spectrum, fit_staged), coupling
      uncertainties and start agreement (zulf_hypothesis.uncertainty), Gaussian coupling priors in refine.
- [x] Joint monotone fits of a concentration series (fit_joint_series: direction and shape free, shared
      nodes, raw FIDs as complex members, leave-one-out prediction, parallel and prior-drawn starts).
- [x] Robustness (16 starts, leave-one-out, two weightings): 1J(C2,H2), 1J(C3,H3) determined; small
      couplings not unique from these spectra.
- [x] Joint fit with 7ad4aafb (volume reading x 0.18 consistent, mole reading not), 48-start fit,
      synthetic recovery tests: the small couplings of these spectra are at an information limit (a
      plausible truth and fits several Hz away score equally); only 1J(C2,H2), 1J(C3,H3) determined.
- [x] Reliability on the data-only score (scripts/reliability_series.py; prior term removed, all fits
      pooled, 3 % set): 5 solutions of the 48-start run; reliable 1J(C2,H2), 1J(C3,H3), J(H2,H5),
      J(H2,H4); literature-start fits +9 % / +26 % on the data
      (docs/analysis/2026-09-30_blake-pyridine_reliability.md).
- [x] Raw FIDs of the series (7 points incl. x 0.10): one processing for all (crop 0.1-8.1 s, 0.3 1/s,
      phase 176.9 deg at the switching edge); the x 0.66 anomaly was Blake's baseline correction. Complex
      joint fit with the model through the same processing (`--real-only false`, no AsLS): residuals
      0.15-0.24, 1J(C2,H2), 1J(C3,H3), 1J(C4,H4) as on Blake's spectra; the AsLS route biases 1J(C4,H4)
      by -3.5 Hz. Small-coupling jump at x 0.50-0.75 remains, with a jump of the free model delay
      (docs/analysis/2026-10-02_blake-pyridine_fid-joint-fit.md).
- [x] One model delay for the series (`--shared phase_delay`): total +0.3 %, delay -1.86 ms; the
      small-coupling jump at x 0.50-0.75 stays (not caused by the per-spectrum delays).
- [x] Without the monotone constraint (`--shape free`, shared delay): total -20 %, jump stays. Two-branch test:
      every spectrum alone prefers branch A (smooth in x, far from the literature small couplings); the
      series-average priors split the series between branches A and B: the jump is a prior artefact.
- [x] Monotone fit started in branch A: total 0.2321 (-22 %, data -25 %), smooth in x, no jump; 1J fall
      smoothly. Several H-H couplings far from literature (3J(H3,H4) 10-14 Hz): best data description, not an
      assignment of the small couplings.
- [ ] Fit with the H-H couplings held at the literature values (solvent-independent within a few 0.1 Hz).
- [x] Uncertainty budget around the branch-A reference (noise, near-equivalent set, six processing variants,
      four synthetic recoveries with the real residual): 1J +-0.1-0.6 Hz, small couplings +-0.3-2.9 Hz, recovery
      dominates; model error (H-H far from literature) not covered.
- [ ] Joint fit without series-average priors (or per-spectrum priors).
- [ ] Pyridinium (single spectrum, same limits expected).
- [ ] Real-data model mismatch: the staged joint fit recovers a synthetic truth (median error 0.7 Hz) but
      not the best known real-data minimum; test line-shape models (Voigt / apodized record) and the C2
      isotopomer (14N relaxation) on the real spectra; ask for the raw FIDs and processing of the series.
- [ ] Correct long-range C-H literature values (Shiner and Wyllie 1973) as prior centres.

## Phase 3e: fine structure (local refinement, then global; user 2026-10-04)

Motivation: on isopropylamine (docs/analysis/2026-10-03_isopropylamine_complex-fit.md) the global fit with one
rate per line cluster reproduces both bands, but not the 121.69 / 122.28 Hz pair with its near-baseline dip
(full record): the fitted transitions 121.85 / 122.00 / 122.68 Hz form a doublet 0.3 Hz too high with the
intensities reversed. Such fine structure carries the small couplings; a global residual is dominated by the
tall lines and barely sees it.

- [x] Residual-driven region finder (JointSeries.find_residual_peaks; peak_sources with shift / split / selectivity; commits 8430327, 7b6417b, 737265a): windows where the residual is structured (sharp peaks or dips the model
      misses, e.g. a data valley the model fills), with the transitions that fall there and their df/dJ
      (zulf_core.physics.derivatives) to name the couplings that can move them.
- [x] Local refinement (JointSeries.local_fit, window forward with line_band_hz and jacobian_only: 26 x faster Jacobian; commits 109cf89, 1e7894b): fit a narrow window (1-3 Hz around the feature) on the highest-resolution data (full
      record, no window, more zero fill), free only the sensitive couplings plus local widths / phase, everything
      else held; multi-start inside the local window (local minima are cheap to enumerate there).
- [x] Local fit with gains, phase and background held at the global values (or one free scale): the window
      forward re-solves them and overstates the local gain (isopropylamine NH2 test: 0.0198 local, 0.122 in
      the full-model window).
      Done: `JointSeries.local_fit(hold_gains=True)`, `--component-search-hold-gains` (commit 04b9df6;
      HeldGainLocalFitTests: the local cost equals the full forward's window residual).
- [~] Propagate (residual-peak stage in fit_joint_series; to fix: sticky windows judged on the starting transitions, cap on the plain-objective loss; local-candidate -> global refit as an option): each local candidate becomes a global start; accept it only if the global objective (all
      bands) improves, else report the conflict (the feature wants J that the rest of the spectrum rejects:
      a model error such as a missing coupling or species).
- [x] Line table (zulf_core.physics.lines, scripts/line_table.py; commits 2dc4c16, 5de437e): exact split of every
      line into coupling contributions (Euler), df/dJ, second-order and cross terms, near-degenerate flag.
- [x] Fit trace and slider viewer (fit_joint_series --trace, scripts/trace_view.py; commit b5ef7fa).
- [x] Interactive J tuner (scripts/j_tuner.py, tests/test_j_tuner.py; user request 2026-10-05): coarse / fine
      sliders with the live objective, view cost, residual, lines and residual peaks; peak sources; suggested
      steps; stoppable refinement of chosen parameters; tuned.json readable by --from-joint.
      To do: show several spectra of a series at once; a rate / delay panel.
- [x] Publication figure (scripts/paper_figure.py, figure_tools.py, zulf_processing.anchor_spline_baseline,
      skills/zulf-figures; user request 2026-10-05). To do: an inset placer that finds empty regions.
- [x] Fit monitor (scripts/fit_monitor.py, fit_joint_series --monitor on by default; user request 2026-10-05):
      objective of every evaluation of every start, stage, best couplings, console; page and terminal view;
      results identical with and without (tests/test_j_tuner.py). To do: the same hook in fit_structure.
- [x] Filled-dip rows and sharp-extrema rows (fit_joint_series --dip-penalty, --peak-max-width; user request
      2026-10-05: one broad line must not cover a resolved splitting, a broad data line may be refined).
- [x] Band diagnosis (scripts/band_diagnosis.py, tests/test_j_tuner.py BandDiagnosisTests; user request 2026-10-05:
      "one look shows the CH2 group is bad; know which part to adjust"): misfit share per isotopologue band,
      levers per bad band (wanted step, band gain, cost elsewhere, selectivity; verdict free knob / conflict /
      weak), and whether the band's own parameters are exhausted (then the model lacks something there).
- [x] Component search in fit_joint_series (--component-search start|end|both, default both; user request
      2026-10-05: "at the start and at the end; the start pass guides the following fit"). Every isotopologue on
      the window it dominates, its own couplings (shared H-H held) and rate families refit on the window from
      random starts. The start pass result becomes the first centre of the multi-start; the end pass searches from
      the best solution and refits globally, keeping the result only if it wins. Test: tests/test_j_tuner.py
      ComponentSearchTests (+6 Hz wrong basin found back). Motivated by N-ethylmethylamine 13CH2 (whole
      objective -29 % from a basin the multi-start never reached).
- [x] Joint fitter: model-line weighting pass (the single-spectrum solver's signal_model_passes), so a model
      line where the data show none is fully penalised (eda_slow83 / eda_slow65: unobserved 15N lines at
      weight 0.2). Done: `--model-line-passes N`, MixtureForward.add_signal_cores (tested against a forward
      built with signal_extra_hz), commit 0ecce84. First use: ethylenediamine two-species (model-only
      177.5 Hz line), eda_2sp_mlp.
- [~] (in progress: claude/eager-franklin-tp94gf, after the fitter features) 14N (spin 1) in spin systems, with quadrupolar relaxation: the worst bands of the amines are the
      N-bonded carbons (N-ethylmethylamine CH2, ethylenediamine CH2); test partially resolved 1J(13C,14N) and
      2J(14N,H). Update 2026-10-05 overnight: after a second species (ethylenediamine) and narrow rate families
      (triethylamine, N-ethylmethylamine), the open misfits are methyl 13C lines near 121-122 Hz with no model
      line (isopropylamine 121.67, N-ethylmethylamine 121.33).
- [ ] Amines: test a common minor species (CO2 carbamate / ammonium) as a second component for all four samples;
      ethylenediamine needs one (eda_2sp_fam198, ratio 0.26), isopropylamine rejects a same-molecule copy.
- [x] Several molecules with free ratios (`--structure` JSON list, commit aa02966); asymmetric motif
      X-CH2-CH2-Y (ABCD, commit 518d2c8); `--from-joint` keeps structure values for couplings a previous fit lacks
      (commit 0ecce84).
- [x] N-H protons kept for symmetric motifs: `fit_staged.group_index` located no group in an AA'BB' component
      (GroupIndexTests, commit a532271). Ethylenediamine with static N-H: one species 0.049, two species 0.0090
      against 0.0051 decoupled; the N-H couplings fit to about 0 (fast exchange favoured).
- [ ] Ethylenediamine species 2 with one tied decay rate (refit): do its -18 Hz geminal and line intensities
      survive? (Fixed-rate check without refit: 0.0103.)
- [ ] Intermediate N-H exchange for systems the exact model cannot hold (ethylenediamine with N-H: 9 spins, no
      reducible group, 35 GiB): an approximation such as averaging over the N-H spin states, tested against the
      exact model on a smaller system.
- [x] Automatic rate-family edges from the line clusters and sharp data lines (now chosen by hand, D46).
      Done: `--family-edges auto` (edges_from_lines, auto_family_edges; commit 8f5a09a). To check on real data
      (local line): auto against the hand edges of tea_fam121, nema_fam121_final, eda_2sp_fam198.
- [x] Fitter: one tied decay rate per chosen component (code for the ethylenediamine species-2 test;
      the real-data refit is the local line's); rates remapped automatically when --from-joint comes from a fit with
      other family edges (fit.json records its edges); scripts/snapshot_fit.py (best monitor vector -> --from-joint
      file). Replaces the scratchpad helpers remap_families and snapshot_fit. Done: rate remap and recorded
      edges (cbf7655), `--tie-rates REGEX` (8caf49b), scripts/snapshot_fit.py (64bcc2e); tests in
      tests/test_j_tuner.py (RateRemapTests, TiedRatesTests, SnapshotTests).
      Figure helpers ported: scripts/run_problem.py, plot_runs.py, plot_components.py (commit d6adbff, RunToolsTests).
- [ ] 2J(C,H) sign: multi-start refit on the positive side with wider bounds (the flip test used one start and
      the +3.7 Hz template bound).
- [ ] Objective terms for fine structure: missing-valley rows (the counterpart of the missing-peak rows) and/or
      a derivative-spectrum term in chosen windows; check that they do not trade the tall lines.
- [ ] Validation: synthetic spectra with known close splittings (does the local-then-global loop recover them
      where the global fit alone does not?), and frozen prediction of independent scan groups (as in the 9.22
      analysis).
- [x] Fixed-by-default couplings: list every coupling the template leaves unspecified (fixed at 0) with its
      |df/dJ| (line_table reports zero couplings); free those that move lines by more than a line width
      (isopropylamine methyl-methyl 4J(H,H) freed: +0.19 Hz, total -5 %).

## Phase 3f: measurements in a field (user 2026-10-06)

- [x] Fitted static field, default zero field (D47; `ParameterPolicy.fit_field`, `fit_joint_series --fit-field`;
      tests/test_solver.py FieldFitTests, tests/test_j_tuner.py FieldOptionTests; commits 913bce5, 4fbac0a).
- [x] Acetonitrile fringe-field run 2026-09-28 (docs/analysis/2026-10-06_acetonitrile_field-fit.md): 13CH3 with
      a fitted field and rate families, objective 0.58 -> 0.044, |B| about 0.06 uT, 1J 136.28 Hz (halves agree
      within 0.007 Hz and 3 nT).
- [ ] Field distribution (gradient over the sample) instead of rate families for field-sensitive lines.
- [ ] Analytic field derivatives (Zeeman directions in physics.derivatives) instead of central differences.
- [ ] Field in the global pattern search (solver.search) and in scripts/line_table.py.
- [ ] Acetonitrile open lines: 83-85, 154.1, about 627 Hz; nitrile 13C isotopologue below 22 Hz; 14N.
- [x] Coupling precision option (D55): fits stop at --precision (default 0.01 Hz; 0.001 finer; 0 = old); refine,
      fit_joint_series, J tuner, Studio. Acetonitrile: 17 instead of 30 evaluations, same result.
- [x] Labelling hypotheses (D54): ProtonGroup.isotope, deuterate_exchangeable, labeling_variants,
      fit_structure_labelings, analyze_sample.py --labeling; synthetic methylamine natural / 15N / 2H all recovered.
- [x] Fitted gyromagnetic ratio in a known field (D53): Protocol.gamma_overrides, fit_gamma / field_fixed,
      fit_joint_series --field / --fit-gamma, gamma_identification; tests GammaFitTests, GammaOptionTests.
- [x] J tuner: live sliders for spectrum parameters (field in nT, decay rates on a log scale, delay in ms;
      tests/test_j_tuner.py TunerSpectrumParameterTests; commit 0035e2b).
- [x] ZULF Studio (D48, zulf_studio/, scripts/run_studio.py; tests/test_studio.py): real-time simulator with
      field and line-width sliders, overlay of a processed spectrum, line table and export, fit_joint_series from
      the sliders with a trace (progress) slider, log, terminal, Python console, JSON API for AI agents.
      First use: a 4-start acetonitrile fit from a 20 / 60 nT field start ended in a worse basin (0.158 against
      0.044): the field start matters (try several, as in the analysis).
- [x] Studio display scale continuous (least squares on the magnitudes, lock option) and figures from Studio
      (paper_figure.py: sampling rate and phase from the series entry, PNG/PDF/SVG/caption, manual-parameter
      label; series entries record source_fid), auto phase of the shown data (model / data). D49;
      tests/test_studio.py DisplayScaleTests, FigureTests, AutoPhaseTests.
- [x] Studio AI assistant (D50): Anthropic (claude-opus-5-5) and OpenAI (Responses API) drive the session through
      its tools; docs/STUDIO.md (user guide, API, credentials) and docs/ENVIRONMENT.md. tests/test_studio.py
      AssistantTests. Not yet run against the live APIs (no keys on this machine).
- [ ] Studio: field-start grid for fits; rate families in the quick-look simulation; MCP wrapper of the API.

## Phase 4: architecture comparison (plan weeks 6-9)

- [ ] Train CNN+Transformer; compare (a) graph search, (b) CNN set, (c)
      CNN+Transformer under identical budgets.

> Note (D27): every coverage number above was measured with the
> hash-seed-dependent split and is inflated by train/validation leakage;
> clean re-runs after the fix are recorded below.

- [x] Clean re-run after D27 (same verification problem, 48k freshly
      pre-rendered spectra, CNN+Transformer, 3000 CPU steps, 3826 s). Validation
      on unseen topology families, steps 1000 / 2000 / 3000: structure@10
      0.031 / 0.039 / 0.039; J@10 (1 Hz) 0.008 at every step; component-count
      accuracy 0.70. Re-evaluation of the saved checkpoint in a fresh process
      reproduces the logged numbers exactly (val loss 2.872). J errors after
      permutation matching, 13C-1H strong: median 2.2 Hz, p90 18.7 Hz (the
      18.7 Hz quoted here earlier was the p90, not the median). The
      earlier structure@10 of 0.45 therefore reflected memorized families.
      At this scale the network does not yet generalize to unseen spin-graph
      topologies; the decisive experiment is a GPU run with far more families
      and steps, with the fresh-process re-evaluation as a standing check.

## Model optimization roadmap (from the second survey, docs/REFERENCES.md)

Ordered by expected gain per effort; each item is measured on the same
verification data before adoption.

- [~] Rerank the top-k beam by short solver fits (residual). Measured with
      `scripts/rerank_eval.py` on 40 simulated validation observations (noisy
      FID, random phase and delay, automatic phasing, CNN+Transformer after
      3000 CPU steps, k = 10, 10 s solver budget per candidate):

      | Decoding | structure@1 | structure@10 | J@1 (1 Hz) | J@10 (1 Hz) |
      | --- | --- | --- | --- | --- |
      | beam, offset head | 0.225 | 0.30 | 0.05 | 0.20 |
      | beam, local expectation | 0.225 | 0.30 | 0.05 | 0.20 |
      | solver rerank | 0.30 | 0.30 | 0.225 | 0.225 |
      | solver rerank, J within 0.1 Hz | 0.30 | 0.30 | 0.20 | 0.20 |

      Reranking moves every structure-correct candidate to rank 1 and raises
      J@1 4.5-fold (two-component samples 0 -> 0.26); expectation decoding
      changes nothing. Coverage here (structure@10 0.30) is below the
      training-time validation (0.45), which uses exact phasing and spectral
      noise; automatic phasing and FID-domain noise cost part of it. About
      100 s per sample on a contended CPU. Next: mutation neighbours (13C/15N
      swap, group size +/-1) and sign variants inside the rerank.
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
- [x] Informal blind tests on NMRduino averages (docs/ANALYSIS_LOG.md):
      7dc9a043 no signal; 7ad4aafb pyridine-consistent (not revealed);
      e3d282da final pair lactic acid / alanine, confirmed L-alanine;
      b683220d final guess lactic acid, confirmed. Models contributed no
      usable structure at these SNRs; hypotheses plus the solver did.

## Phase 7: J to structure (docs/J_TO_STRUCTURE.md; user 2026-10-05: "both routes", keep the direct readout)

- [x] Keep the hypothesis pipeline (spectrum -> fitted fragments) unchanged as the first route.
- [x] Route A, rule-based direct readout (`zulf_hypothesis.j_structure.rank_structures`): enumerate heavy-atom
      graphs over the carbon copies and up to two unseen atoms, with bond orders, and rank them by
      p(J | bond count) and p(1J | hybridization). CLI `scripts/j_structure.py`; real networks in
      `configs/j_networks/` (all four amines rank first). Tests: `tests/test_j_structure.py`.
- [x] Paired data from the generator (`zulf_model.structure.observations.observation_from_graph`, with bond counts
      and hybridization labels).
- [x] Route B, learned likelihood (`zulf_model.structure.edge_model`): bond-count and hybridization MLPs trained on
      generator pairs and plugged into the same ranking. Benchmark `scripts/j_structure_benchmark.py`: synthetic
      top-1 0.91 / top-3 1.00 for both routes; residual failures are equal-bond-count ties.
- [ ] Route B as a GNN over the whole network; compare with A.
- [ ] Per-coupling sigma from the fit uncertainty; 15N and 13C-13C observations.
- [ ] Readout on a structure-free J network from the hypothesis pipeline (closes spectrum -> J -> structure).
- [x] Reduced coupling constants K (D51): K next to J in fit.json and J_table.csv; conversion functions for
      15N -> 14N and 1H -> 2H starts (zulf_core.nuclei). Tests: ReducedCouplingTests, RunToolsTests.
- [x] Route B on J vs K (D52; docs/analysis/2026-10-07_j-to-structure_j-vs-k.md): equal on 13C / 1H data; K alone
      transfers to deuterated couplings (0.843 vs 0.298). scripts/j_k_learning_comparison.py.
- [x] Route B in both modes every run, J and K side by side (D56): one model file per mode, `isotopes` in the
      observation, scripts/j_structure.py prints A, B-J, B-K. Tests: BothModesTests.
- [ ] 15N units and N couplings in the J observation (needed for K to matter on 15N-labelled samples).

## Phase 8: ZULF Studio in four modes (owner 2026-10-08)

Studio is organised by task; the model (structure or spin system, couplings, field, line width) is shared by every
mode, each mode shows only what its task needs (D58).

| Mode | Input | Left | Centre | Right | Output |
|---|---|---|---|---|---|
| Simulate | model only | model, couplings, field and line width | simulated spectrum (each component, weighted sum, overlays) | simulation settings (quick look / exact rendering through an acquisition, window, zero fill, Gaussian broadening), export | spectrum, line table, figure |
| Process | scans or an FID, no model | processing recipe (crop, record, window, zero fill, phase, delay) | FID and spectrum, live with the recipe | scan selection (keep checkboxes, quality metric), save | averaged FID, processed spectrum (series.json), recipe |
| Fit | spectrum and a known model | data, model, couplings, field | data, model, residual | fit, figure | couplings, field, figures |
| Blind analysis | spectrum or FID, unknown model | data | spectrum, ranked candidates | analysis settings, report | candidates; one goes to Fit |

Defaults chosen with the owner (2026-10-08, changeable): in a custom spin system the component weights are fitted
amplitudes (can be fixed); several couplings may share one variable; acquisition settings and Gaussian broadening
live in Simulate, not on the Fit screen; ZULF_NMR_Suite stays as it is, its molecule folders can be imported.

Steps:
- [x] 8a. Mode bar; Simulate, Fit, Blind analysis on the present widgets; the session records its mode (done, commit below; Process shows the import entry points until 8b).
- [x] 8b. Process mode (done): scan table from average_scans metrics, recipe with live preview (debounced), averaging
  of the kept scans as a job, save the spectrum and the recipe (never into the data folder).
- [x] 8c. Custom spin system (done; Studio editor, API, tests): isotopes, J matrix (upper-triangle table or text) with variables, component list
  with weights; sliders per variable; save / load models; import ZULF_NMR_Suite molecule folders.
- [x] 8d. Fits of a custom spin system (done, D60) (fit_joint_series model input {"spin_system": ...}; variables as fit
  parameters; weights as amplitudes).
- [ ] 8e. Exact rendering in Simulate (the processed forward model of zulf_core.render with an acquisition);
  overlays of several parameter sets.

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
