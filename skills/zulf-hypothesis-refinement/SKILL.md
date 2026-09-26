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

## Building hypotheses

- One component per isotopologue; tie every coupling that the isotopologues
  share (all proton-proton couplings) and every symmetry-equivalent coupling
  (`RefineSettings.ties`). Map group indices per component explicitly: the
  heteronucleus changes the group order.
- Start from literature-like values, say where they come from, and give
  uncertain small couplings wide bounds (coupling_margin_hz 6-12).
- Include every isotopologue that can put lines in the fitted range; state
  which ones are omitted because their lines fall outside it.

## Fitting

- Default: complex data, `gain_model="shared_phase"`, delay fitted,
  background_order 1, 3 starts, analytic (Kaufman) Jacobian.
- Prefer `band_weighting="signal"`: noise-scaled chi-square, points outside a
  data-driven signal mask down-weighted (default 0.2), and the result reports
  `signal_region_residual`. Check the mask: a sharp instrument or background
  feature (e.g. 72-76 Hz on the NMRduino) is detected as signal; exclude it
  with the fit ranges until reference-based exclusion exists.
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
- Stop runs by PID (process and children); `pkill -f` can match the shell.
- 8-spin isotopologues with 3 starts take about 10 minutes each on one core.
