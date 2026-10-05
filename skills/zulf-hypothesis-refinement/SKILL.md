---
name: zulf-hypothesis-refinement
description: Turn model proposals and chemical hypotheses into refined, physically checked ZULF candidates - running and refining model candidates, building multi-isotopologue hypotheses with tied couplings, complex fitting settings, sign variants, blind global search, physical consistency checks, fair model comparison and practical compute rules. Use whenever candidates are fitted and ranked against experimental data.
---

# Refining and ranking hypotheses

Checklist of red flags and comparison rules: `references/checks.md`.

## Model proposals

- Run every checkpoint (`ModelProposer`) with a solver-fitted phase and with
  the automatic one; record each checkpoint's training range. A model cannot
  propose spin counts, nuclei or component counts outside it, and it may
  saturate at the maximum spin count.
- Deduplicate by groups and couplings >= 1 Hz, refine every distinct
  candidate, rank by residual, then apply the checks. At CPU scale the models
  gave skeletons but no small couplings, and nothing on weak data. Trust the
  models only as far as their clean (leak-free) validation says.

## Automatic fitting (D43)

- A known or proposed structure: `fit_structure(fragment, observed,
  report=prefix, route="both")`. It fits every exchange regime of N/O/S-H
  protons (dropped = fast, kept = slow), the variants free / fixed / ratios
  (ratios: natural amplitudes, one rate per isotopologue), thorough starts
  (perturbed + global pattern search), then the triggered moves (coupled
  proton, proton count with 1J rescaled, free remote couplings, Gaussian
  line width) on one yardstick, and writes the report (ranked table, J
  matrices with held couplings marked, fit-phased figure). Do not propose
  these variants by hand; read the table and the findings.
- After any change to zulf_hypothesis or the solver run
  `ZULF_DATA_DIR=... python scripts/regression_confirmed.py --mode known`
  (and `--mode blind` for search changes) and compare with the previous
  summary in ANALYSIS_LOG; add every newly confirmed sample to
  `configs/confirmed_samples.json`.
- Solver rendering uses `Acquisition.local_record()` (exact); a fit of a
  9-spin, 3-isotopologue model takes about 10 s per start on one core.

## Building hypotheses

- Start with `propose_hypotheses(observed, instrument_hz, providers=...)`:
  it lists bands, X-Hn group candidates per band, and ranked fragment
  proposals built into isotopologue models (weak bands may stay open and are
  listed). Neural checkpoints join as hint providers
  (`zulf_hypothesis.adapters.model_hints`); read `ProposalSet.insight` for where
  they agree or disagree with the data, and never rank by hint confidence.
- Labelling: natural abundance is the default (`Labeling.natural`); for an
  enriched sample pass `Labeling.enriched({"13C": 0.99})` or site-specific
  levels, and ask the user for them. Doubly labelled isotopologues are
  weighted by abundance, tied to a parent (no free amplitude or rate) and
  omitted when their lines cannot reach 2 sigma; read `model.omitted`.
- `propose_hypotheses` also runs the motif scan (whole motifs such as
  ethyl, CH2-CH2, benzene and pyridine rings, with generic couplings and 1J
  from band positions); read `motif_table()` next to `table()`. When a new
  structure type is confirmed, register it as a motif and as a benchmark
  case, and run `python -m zulf_hypothesis.benchmark`.
- Then `search_hypotheses(observed, proposals, SearchSettings(...))`: it
  refines the top proposals in free and fixed-abundance variants (parallel
  workers), scores all on one yardstick (data cores, BIC), runs the checks,
  extends the best clean ones by triggered moves (accepted only when the BIC
  improves by 6), matches known compounds, and returns a ranked table and a
  log. Report the table with its findings; `best` prefers hypotheses without
  warnings. Fits that switch an isotopologue off (`collapsed_component`)
  rank last (`demoted` in the table; 4322bdfc: triethylamine lost to a
  CH-CH3 fit with a 0.06 CH component before D42). Motif-scan proposals get
  a short screening refinement before the `top_motifs` cutoff (read the
  `screen` log entries), and `change_proton_count` turns CH <-> CH2 <-> CH3.
- Build isotopologue sets with `zulf_hypothesis` instead of by hand:
  describe the fragment by labels (`Fragment`, or a template from
  `fragments.TEMPLATES`), then `build_model(fragment, ranges=...)` gives the
  interpretation, ties, natural-abundance ratios and the omitted
  isotopologues with reasons; `model.settings(base, fixed_ratios=...,
  shared_rate=...)` gives the RefineSettings; `model.named_couplings(params)`
  reports couplings by label ('J(Ca,Ha)'); `run_checks(model, summary,
  peers)` applies the red flags; `propose_all(fragment, findings)` gives the
  triggered one-step extensions. New motifs, checks, moves and reference
  compounds are registry entries (see docs/ARCHITECTURE.md).

- One component per isotopologue; tie every coupling that the isotopologues
  share (all proton-proton couplings) and every symmetry-equivalent coupling
  (`RefineSettings.ties`). Map group indices per component explicitly: the
  heteronucleus changes the group order.
- For a pure sample, fix the isotopologue amplitude ratios to natural
  abundance (`RefineSettings.amplitude_ratios`, for example 1 : 1 : 0.34 for
  two 13C sites and one 15N site) and tie the decay rates across
  isotopologues (`ties` on `cN.log_rate0`). With free amplitudes and rates
  one component can turn into broad background (e3d282da: C-alpha at
  7.7 1/s and 7.5x amplitude cut a 139 Hz hole that cost little outside the
  cores). Free ratios remain a check, not the default.
- Start from literature-like values, say where they come from, and give
  uncertain small couplings wide bounds (coupling_margin_hz 6-12).
- Include every isotopologue that can put lines in the fitted range; state
  which ones are omitted because their lines fall outside it.

## Fitting

- Default: complex data, `gain_model="shared_phase"`, delay fitted,
  background_order 1, 3 starts, analytic (Kaufman) Jacobian.
- Prefer `band_weighting="signal"`: noise-scaled chi-square; the weight is 1
  at peak cores and falls off smoothly (Gaussian, `signal_taper_hz` 2 Hz) to
  `signal_outside_weight` (0.2) far from any core - no hard edges. Cores come
  from the data first (narrow excess above `signal_threshold` sigma). The
  model's cores come from its transition lists, not from peak picking the
  rendered spectrum: the incoherent line envelope
  E(f) = sum_k abs(g_c a_k) P_R(f - f_k) (P_R the rendered unit-line
  magnitude profile for the line's rate; `MixtureForward.model_envelope`)
  joins the cores where it exceeds `signal_model_threshold` (2) sigma. It
  covers isolated lines of either sign, many broad overlapping lines, and
  holes where lines cancel. The fit is repeated until the cores are stable
  (`signal_model_passes=3`; flag `signal_mask_model_pass`). Per-line peak
  heights alone missed a broad component (e3d282da: 50 lines of about
  0.2 sigma each summing to 9 sigma and cutting a hole at 139 Hz).
  The result reports `signal_region_residual` (data and model cores),
  `data_region_residual` (data cores only: compare models on this one, it
  is the same region for every model) and `signal_model_lines_hz`. Check the mask: every region must be either an
  assigned line or an explicitly unresolved feature. A data-only mask cannot
  tell molecular lines from sharp background (e3d282da: 72-76 Hz, still
  unresolved between instrument background and a 15NH3+ line near
  abs(1J(N,H)) = 73 Hz); keep such regions in the fit and report them as
  unresolved.
  Plot the fit over the whole range, not only the mask: a data-only mask with
  a low outside weight lets a model put spurious lines outside the mask
  (e3d282da: 263-265 Hz); the cores must cover the model's own lines.
- Weak data: fit the truncated window from `zulf-fid-processing`; starts at
  narrow lines (initial rate 0.5-1/s) with continuation (1, 0).
- Sign variants: only when signs are the question; they multiply the cost up
  to 16x.
- Blind global search from a model skeleton measures the pipeline; it found
  an envelope fit with unphysical couplings, not the structure.
- Save full results (parameters, gains, prediction) for every fit.

## Spectra processed elsewhere and concentration series (Blake pyridine, ANALYSIS_LOG)

- Processed real spectra (frequency and value arrays): `scripts/fit_processed_spectrum.py` (fit_structure) or
  `scripts/fit_staged.py`; record=None gives ideal Lorentzian lines with analytic derivatives. With
  `--record fs,points` every derivative needs an FFT of the record length (88457 = 53 x 1669 points: 8x
  slower); the Gaussian-width move needs the record route, so leave it out for long records.
- Without a FID there is no global pattern search: the starting couplings decide the minimum. Use accepted
  values as starts, hold the H-H couplings first and free them afterwards (staged): a direct all-free fit from
  the same values ended in a worse minimum (pyridine x 1.00: J(H3,H4) 3.7 Hz vs 7.2 Hz staged, lower chi2).
- The user decides what is held: couplings are the result, so do not fix them; keep them free with Gaussian
  priors (RefineSettings.priors, weight ~ reduced chi2) when a direction is weakly determined.
- A concentration series: `scripts/fit_joint_series.py`, every coupling monotonic in the concentration with
  fitted direction and free shape (softmax steps). Never a linear dependence (user). Spectra at one
  concentration (two instruments) share a node; raw FIDs can join as complex data.
- Judge couplings by reproducibility over near-best starts, jackknife over leave-one-out fits, robustness to
  the weighting and plausibility, not by the residual and not by linearised errors (they understated the
  spread 5-20x). On the Blake series only 1J(C2,H2) and 1J(C3,H3) passed every test.
- Compare fits with different prior centres on the data alone: the joint-fit score includes the prior term,
  so a fit started from other centres can look worse only because of its prior. `scripts/reliability_series.py`
  removes the prior term from every stored solution (checked against a direct recomputation), pools all fits
  and takes the solutions within 3 % of the best data score (Blake series: the literature-start fits were
  +9 % and +26 % on the data, the staged fit +8 %; the set was 5 solutions of the 48-start run).
- Weighting of weak peaks: the 4 sigma core threshold left weak lines at the outside weight 0.2; the user chose
  2.5 sigma and a 4 Hz fall-off (`--signal-threshold`, `--signal-taper`). Plot the weights (figure) before
  a long run.
- Before believing a coupling, run a synthetic recovery test (`scripts/synthetic_recovery_series.py`: a
  plausible truth, the measured grids and noise, the same fit). Compare the truth's score with the best found:
  a best found above the truth's score is a search failure. On the Blake-like series the plain fit ended in
  compensating minima (errors up to 6 Hz) because the 1J starts were constant while 1J changes with the
  concentration; `--hold-small-first` (fit 1J, rates and delays with the small couplings held, then free all)
  recovered the truth (median error 0.7 Hz, 1J within 0.6 Hz). Smoothing continuation, coordinate scans and
  differential evolution did not help. Fix the large couplings first, then the small ones.
- Protonated forms keep the molecular symmetry (fixed in `protonated`): without it the duplicated
  isotopomers over-fitted and a neutral pyridine spectrum chose the protonated form.

## Checks before ranking

Reject or flag fits with parameters at bounds (especially decay rates at the
lower bound), abundance ratios that contradict natural abundance, unequal
couplings that must be equal, implausible sizes, collapsed components, or
gains that come from background regions or band edges. Prefer the model that
explains every clear line with consistent physics; test nested models with an
F-type comparison (zero-filled points are correlated). A lower residual alone
is not evidence.

## Compute practice

- Parallel fits: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`, at most one
  process per core; order queues so hypotheses run before long model lists.
- Stop runs by PID (process and children), saved when the run starts. Never
  select processes by command-line text (`pkill -f`, `ps | grep | kill`): the
  calling shell's own command line contains the same text and gets killed.
- Queue follow-up runs from a small launcher script file, one queue at a
  time; two waiters on the same trigger start duplicate runs.
- Module-level settings (for example fit ranges) that several imported helper
  scripts assign must be set after the last import, and the run must print
  its effective settings (ranges, signal mask) so mistakes show at once.
- 8-spin isotopologues with 3 starts take about 10 minutes each on one core.

## Chemical exchange (N-H, O-H) with a fitted rate

- Limits first: fast exchange (protons dropped) and slow exchange (static
  spins) are fitted by `fit_structure`. For an intermediate regime, fit the
  exchange rate: `scripts/fit_staged.py ... --protonate --exchange HA1
  --kex-start 10` (label = the proton group; one rate tied across the
  isotopologues; extra starts at 1, 100, 3000 1/s are added).
- The exchange model is a Liouville-space one (D45): about 3 s per
  evaluation for 7 spins on CPU, 65 % in one eigen-decomposition. On a
  machine with a GPU use `--device cuda` (or ZULF_LINALG_DEVICE=cuda); torch
  must be installed. 8 spins are about 10x slower; more are impractical.
- The incoming solvent spin carries the preparation polarization by default
  (the solvent was prepolarized too); use solvent_weight=0 in
  `exchange_transitions` for an unpolarized pool.
- Linearised uncertainties are not computed for exchange fits yet.
- Fast vs slow with the series fitter (`fit_joint_series.py --shape free --exchange fast|slow`, one spectrum):
  start the slow variant also nested, from the fast solution with every coupling to the exchangeable protons at 0
  (`--start-couplings`). On isopropylamine 64 generic slow starts (N-H couplings at their motif defaults, e.g.
  3J(H,N-H) 5.5 Hz) ended 26 % above the fast fit, while the nested start ended 13 % below it. Then check what
  the extra freedom is used for: N-H couplings below 1 Hz and a broad 15N line without a clear data counterpart are
  not evidence of a resolved NH group (2026-10-03_isopropylamine_complex-fit.md).


## Logging an analysis

- Start `docs/analysis/YYYY-MM-DD_<sample>_<topic>.md` from
  `docs/analysis/TEMPLATE.md` when the analysis starts; fill in commands,
  numbers, figure paths and the conclusion as results arrive; add its row to
  `docs/analysis/README.md` and a one-line entry to `docs/ANALYSIS_LOG.md`.
- The fitting scripts write `RUN_LOG.md` (command, commit, times, results)
  into their output directory; cite it in the analysis log.

## Fine structure: from a misfit to the couplings (docs/WORKFLOW.md, section 8)

Use occasionally or as the last step, not inside every fit.

- Start with `scripts/band_diagnosis.py <fit options> --fit RUN/fit.json --figure OUT.png`. It shows which
  isotopologue band carries the misfit and, per bad band, the parameters that act on it:
  - **free knob**: the band gains and the rest barely changes. Refit it.
  - **conflict**: the band and the rest of the spectrum want different values of a shared parameter (an H-H
    coupling, a 1J that also places other lines).
  - **weak**: no real gain.

  "Already at the band optimum" can be a local optimum. Isolate the component's window and search its own
  couplings globally from random starts (`JointSeries.local_fit` with the shared H-H couplings held), then refit
  globally from the best basin. N-ethylmethylamine 13CH2: the whole objective dropped 29 % (3J(C2,H3) 12.1 -> 5.8,
  1J 138.4 -> 131.6 Hz). Only when that search also fails, and every lever is weak or a conflict, will no
  coupling change fix the band. The model lacks something there: a fixed or tied
  coupling, a line, another nucleus. Do not tune further; change the model.
- **Slow N-H exchange check.** Slow exchange implies a natural-abundance 15N-H line: at |1J(N,H)| for NH, and
  near 1.5 |1J| for NH2. Check that the data show it before taking a slow-exchange gain. The joint fitter weights
  a model line where the data show none at only 0.2. Amines, 2026-10-05: N-ethylmethylamine has no line at
  65 or 83 Hz, and ethylenediamine slow fits at 65 and 83 Hz score alike (1J not determined).

- `scripts/line_table.py --structure ... --fit RUN/fit.json --band lo,hi [--second-order]`: every line with df/dJ
  and the exact split f = sum J df/dJ (zero and tied couplings included). Couplings fixed at 0 by the template
  show up here with their df/dJ: free them when they move lines by more than a line width.
- Second-order terms are for the report; near-degenerate lines (another line within 0.1 Hz) need an exact
  recomputation, not a Taylor step (isopropylamine 122.13 / 122.15 Hz: 1 Hz step off by 0.15 Hz at 2nd order).
- `JointSeries.find_residual_peaks` (both signs, prominence over a robust local noise, narrow, assignable only
  with a model line within 0.6 Hz), `peak_sources` (rank couplings by `split` for a dip or a doublet, by `shift`
  for a displaced line; a one-bond coupling only shifts), `local_fit` (window-only forward, few free parameters,
  fast), then a global refit that must improve the global objective.
- `local_fit` re-solves gains, phase and background on the window: compare candidates by the full-model window
  cost, not the local cost. Couplings that start at 0 between an uncoupled group and the rest (protons added
  when switching fast -> slow exchange) sit at a stationary point: give them nonzero starts.
- Before blaming couplings for a doublet, narrow the rate family of the feature: one peak where the data have
  two means a position problem. A species added for a test (e.g. 15N at natural abundance) must not create
  lines the data do not show; check the cost per frequency bin against the baseline.
- Long fits: `python scripts/fit_monitor.py runs/processed` (or `--text`) shows every start live (objective,
  stage, best couplings, console); stop a run early by PID when all starts sit in the same basin.
- One broad line over a doublet (rate families make broadening cheap): add `--dip-penalty` and
  `--peak-max-width 1` with a low `--peak-prominence` (0.02) so every sharp line and valley gets a row.
- Hand tests: `scripts/j_tuner.py` (the fit's command line + `--fit RUN/fit.json`) shows the fitter's objective,
  the view cost and the residual peaks live while a coupling slider moves; refine the ticked couplings from
  there and save tuned.json as a --from-joint start. Judge a move by the objective, not the view alone.
- The residual-peak stage (`--residual-peaks`) still lets assignability flip between rounds and may accept a
  large plain-objective loss; check `fit.json["residual_peaks"]` (rounds, scores) before using its choice.
- `--trace N` and `scripts/trace_view.py RUN` to watch a fit step by step (slider).


## From a fitted J network to a structure (Phase 7, docs/J_TO_STRUCTURE.md)

- Write the fitted couplings as an observation, following configs/j_networks/*.json. It needs units (13C sites
  with proton count and copies), proton groups, and J(unit, group) / J(group, group).
- Rank with route A: `python scripts/j_structure.py OBS.json --explain 1`. Add route B with
  `--model runs/models/j_edges_v1.json`; train it first with `scripts/j_structure_benchmark.py --model ...`.
- Read the margin.
  - A margin equal to a prior difference (1.0, 2.0) with equal log L means the J values cannot separate the two
    graphs. This typically happens when an unseen atom X is placed differently, or a ring passes through X
    instead of two separate X atoms. Report both graphs, not one.
  - The readout says where the carbons are and that something sits at a free valence, not what it is (N, O or
    a quaternary C).
- A network from a fit of the true structure is biased toward it. For a real test, fit the J network without a
  motif and then rank it.
- The hypothesis pipeline (spectrum -> fragments) stays the primary route. The readout proposes candidates for
  it and checks it.
