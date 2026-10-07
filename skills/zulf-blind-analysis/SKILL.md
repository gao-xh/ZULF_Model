---
name: zulf-blind-analysis
description: End-to-end workflow for identifying an unknown sample from an averaged ZULF NMR FID - orders the processing, interpretation, phasing and refinement skills, sets the ground rules and the report format, and keeps the case notes and full derivation paths of earlier blind tests (including the confirmed L-alanine case). Use when a user hands over an experimental FID and asks what it is.
---

# Blind analysis of an experimental ZULF FID

This skill orders four category skills; follow each one's own SKILL.md:

| Step | Skill | Output |
|---|---|---|
| 1. Diagnose and process | `zulf-fid-processing` | processed spectra at several windows, instrument lines and background marked, residual scale |
| 2. Read the spectrum | `zulf-spectrum-interpretation` | line list, group assignments, isotopologue expectations, hypotheses with predictions |
| 3. Propose and refine | `zulf-hypothesis-refinement` | refined model candidates and hypotheses, physical checks, ranking |
| 4. Phase (display and cross-check) | `zulf-phasing` | unbiased phased spectrum, delay uncertainty, phased cross-check |
| 5. Report | this skill | ranked candidates, refined couplings, figures, what data would decide |

Earlier cases: `references/case_notes.md` (outcomes) and
`references/derivation_path.md` (the step-by-step reasoning, including the
steps that were wrong). Read the derivation path of a similar case before
starting; follow the same observation -> inference -> test -> outcome format.

## Ground rules

- Every number comes from a tool; every result is a conditional candidate.
- Separate data from prior knowledge; label values "from memory".
- Report residuals with their scale (no-lines residual and noise floor of the
  same window, range and mode).
- Test a hypothesis only with data processed without it (no phase taken from
  its own fit).
- Screen every feature against instrument lines and a no-signal dataset
  before interpreting it.
- Prefer physically consistent candidates over marginally lower residuals.
- Report what the data cannot decide, and which measurement would.
- Ask for the sample state when it can be known (solvent, pH, pure or
  mixture, labelled or natural abundance); it fixes exchange models and
  isotopologue abundance ratios.

## Working practice

- Start every new sample with `python scripts/analyze_sample.py FID --id ID`
  (standard processing, overview figure, blind search, report with J
  matrices and a fit-phased figure); after the reveal run it again with
  `--structure` for the confirmed structure, add the sample to
  `configs/confirmed_samples.json` and run the regression.

- Start with an overlay on confirmed samples and a no-signal dataset,
  processed identically; then fit the simplest fragment, then extend it one
  step at a time on the same data, window and weighting.
- Run independent starts or variants as separate processes (one per core);
  a three-start 8-spin fit takes about 10 minutes that way.
- Check a claimed mechanism (a weight, a mask, a penalty) on the actual fit
  before reporting it; one claim here ("spurious lines are fully
  penalised") was wrong until the model cores came from the transition
  lists.
- Log each step in docs/ANALYSIS_LOG.md as it happens, and update the skill
  when it teaches something reusable.

## Isotope labelling (D53, D54)

- Ask the sample provider whether the sample is labelled (15N, 2H, 13C) or in D2O: it does not reveal the
  structure, and labelling is deliberate. Natural abundance means 13C (1.1 %) and weaker 15N (0.37 %); 2H is not
  visible.
- Unknown labelling: `analyze_sample.py --structure ... --labeling unknown` fits natural abundance, 15N enrichment
  (if the structure has N) and deuterated exchangeable protons (if it has N-H / O-H) and ranks them on one scale
  (labelings.json). In fast exchange the 1H and 2H hypotheses fit identically (decoupled either way): a tie, not a
  result.
- A spin-1/2 heteronucleus can be named from the data in a known field: fit its gamma (`fit_joint_series --field
  ... --fit-gamma 13C`, `gamma_identification`); a change of spin (1H / 2H, 15N / 14N) needs the hypothesis
  comparison above.

## Report format

1. Processing and assumptions (sampling rate, crop, SG, window).
2. Line list with SNR; instrument lines and background excluded.
3. Ranked table: residual on the stated scale, parameters at bounds,
   isotopologue amplitude ratios, one-line verdict.
4. Refined couplings of the leading candidates, with their physical checks.
5. Figures: data and fits (magnitude or complex; phased real when the phase
   is sound), and any check that rejected an alternative.
6. Conclusion as a short candidate list and the decisive next measurement.

## Lessons from the confirmed case (pure L-alanine in water, sample e3d282da)

- The workflow kept alanine in the final pair (with lactic acid); the
  remaining ambiguity was a data limit, not a workflow error.
- Decisive steps were the truncated window and the natural-abundance
  amplitude test (rejected isopropyl).
- The 72-76 Hz feature stays unresolved (instrument background or a 15NH3+
  line under intermediate exchange). Lesson: a background check must compare
  the shape at the same position, not a band-averaged level; keep features
  "unresolved" until a blank or a chemistry change (pH, D2O) decides.
- The CPU-scale models contributed nothing on this sample; chemistry-guided
  hypotheses plus the solver did the work.

## Lessons from the confirmed case (lactic acid, sample b683220d)

- Final guess correct. Decisive steps: two bands with free amplitudes at
  1 : 1 (one CH, one CH3; isopropyl rejected again by abundance), the
  CH-CH3 couplings (3J(HH) 7.06 Hz), and a same-instrument comparison with
  a confirmed sample: 1J(CH) 147.1 vs 145.4 Hz in alanine (O vs N on the CH)
  and no 72-76 Hz feature (no 15N-H).
- Keep a library of confirmed samples processed identically and overlay new
  data on them first: line shifts of 1-2 Hz between related compounds are
  obvious in an overlay and are good substituent probes.
- A decay rate on one isotopologue much larger than on the other (H7 on
  lactic acid: 4.7 vs 1.9 1/s on the methyl carbon) marks missing couplings
  on that isotopologue; adding one weakly coupled spin (H8) removed it and
  halved the misfit. The extra spin's couplings were not determined (two
  fits gave different sets), so report such a spin as "required, identity
  open" unless its couplings are stable.
- Solver defaults used here: band_weighting="signal" with envelope model
  cores, free amplitudes first (they are the abundance test), then fixed
  natural ratios with a shared rate as the physical check.

## Lessons from the confirmed case (ethylenediamine, sample fde3fbb2)

- The automatic program found the correct core (symmetric CH2-CH2) with no
  hand input; the hand guess added a wrong substituent (carboxyl).
- 1J(CH) identifies the carbon type (CH2, sp3), not the neighbour: CH2 next
  to NH2 and next to COOH differ by a few Hz at most. Name a substituent
  only with direct evidence (15N lines, N-H couplings, an overlay with a
  confirmed reference); otherwise report "CH2-CH2 core, substituents open".
- Amines with fast N-H exchange look like bare carbon skeletons; ask for pH
  and solvent, and propose a low-pH or D2O measurement to expose N-H.
- A large c_hat and unassigned weak bands mean the model is incomplete:
  list them as open features rather than explaining them away.

## Lessons from the confirmed case (triethylamine, sample 4322bdfc)

- Read every band's X-Hn position before trusting the ranked table: a weak
  band at 1.5 J marks a CH2 carbon even when the top hypothesis has a CH.
- A free-variant winner with an isotopologue near zero (abundance finding)
  is not a structure; compare the fixed-abundance variants among
  themselves and report them next to the free winner.
- Quick-scan motif ranks rest on 1J picked from band positions and generic
  couplings; a motif just below the refinement cutoff can be the answer.
- Extra couplings through a heteroatom (3J(C,N,C,H)) distinguish a
  tertiary amine's N-CH2 from ethanol or ester ethyls.

