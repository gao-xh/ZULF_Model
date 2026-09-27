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
