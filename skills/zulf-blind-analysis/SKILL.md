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

## Report format

1. Processing and assumptions (sampling rate, crop, SG, window).
2. Line list with SNR; instrument lines and background excluded.
3. Ranked table: residual on the stated scale, parameters at bounds,
   isotopologue amplitude ratios, one-line verdict.
4. Refined couplings of the leading candidates, with their physical checks.
5. Figures: data and fits (magnitude or complex; phased real when the phase
   is sound), and any check that rejected an alternative.
6. Conclusion as a short candidate list and the decisive next measurement.

## Lessons from the confirmed case (L-alanine, sample e3d282da)

- The workflow kept alanine in the final pair (with lactic acid); the
  remaining ambiguity was a data limit, not a workflow error.
- Decisive steps were the truncated window, the natural-abundance amplitude
  test (rejected isopropyl) and the no-signal background check (withdrew a
  spurious 15N improvement; L-alanine in neutral water has fast NH3+
  exchange, so no 15N-H coupling is expected).
- The CPU-scale models contributed nothing on this sample; chemistry-guided
  hypotheses plus the solver did the work.
