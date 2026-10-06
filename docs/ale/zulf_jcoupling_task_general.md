# ALE-V2 task (general version, paste-ready): J-coupling network from a raw zero-field NMR FID

One workflow and one grading logic, applied to several molecules as variants (ALE: "one main.py is one
workflow plus one general grading logic; variants reuse that logic over different concrete inputs").
Each section below is meant to be pasted into the matching field of the submission form
(https://agents-last-exam.org/submit). Status: draft, not submitted; open points at the end.
The single-molecule draft is docs/ale/2026-10-06_zulf_triethylamine_task.md.

---

## Title

Determine a molecule's J-coupling network from a raw zero-field NMR free-induction decay

## Summary

Zero- to ultralow-field (ZULF) NMR measures J couplings directly, without chemical shifts, but there is no
standard software to analyse it. Given the raw signal-averaged free-induction decay (FID) of a small organic
liquid at natural 13C abundance and its structure, the agent must process the FID, build a zero-field spin
simulation from first principles, fit the visible 13C isotopologue spectra and report the one-bond and
long-range J couplings with their relative signs. The same task applies to any molecule; the variants use four
experimentally measured samples.

## Domain

Chemistry / physical chemistry: NMR spectroscopy (zero-field J-spectroscopy), spin physics, spectral fitting.

## Software

Python 3.10 with NumPy, SciPy and Matplotlib (preinstalled). No NMR or spin-simulation packages, no network.

## Task description (prompt shown to the agent)

You are working on a Linux VM.

A zero-field NMR spectrometer recorded the free-induction decay of a neat liquid sample at natural isotopic
abundance. Determine the J-coupling network of the molecule from this FID.

Visible input:
- `input/fid.npy`: averaged FID, float64, sampled at the rate given in `input/acquisition.md`.
- `input/acquisition.md`: how the data were taken. Prepolarization in a permanent magnet, transfer into a
  magnetic shield, sudden switch-off of the guiding field, detection of the z magnetization by an atomic
  magnetometer. The initial polarization and the detected signal of each nucleus are proportional to its
  gyromagnetic ratio. Sample 0 is the first sample after the field-switch command; the detector saturates for
  a few milliseconds and rings for tens of milliseconds after the switch. Known interference: mains harmonics
  (multiples of 60.06 Hz) and listed instrument lines.
- `input/structure.md`: the molecule (name and SMILES), the labels of its carbon sites and proton groups, the
  list of couplings to report, and the output schema.

Produce:
- `output/results.json` with every requested coupling in Hz, with sign (one-bond 13C-1H couplings are positive
  by convention; report the other couplings with their signs relative to them), the isotopologues you modelled,
  and short notes on your processing and fitting choices.
- `output/fit_figure.png`: the processed experimental spectrum and your simulated spectrum over the signal
  bands.

Rules: do not modify `input/`; write only under `output/`; do not install packages; the coupling values must
come from your analysis of this FID, not from literature or databases.

## What the agent must do

1. Diagnose the raw FID (field-switch edge, saturated and ringing samples, baseline drift, noise) and choose
   the crop, record length, baseline removal and apodization.
2. Compute and phase the spectrum (zero-order phase and the delay between the field switch and the start of
   free evolution).
3. Decide which natural-abundance isotopologues are visible (one 13C per molecule; 13C-13C, 15N and 2H species
   negligible), which proton groups are magnetically equivalent, and whether exchangeable N-H / O-H protons are
   coupled or decoupled by fast exchange.
4. Implement the zero-field Hamiltonian H = 2 pi sum J_ij I_i . I_j and the sudden-drop observable to obtain
   line frequencies and amplitudes for each isotopologue.
5. Fit couplings, line widths, amplitudes and phase of all visible isotopologues jointly to the processed
   spectrum, rendering the model through the same crop and window as the data.
6. Determine the signs of the long-range couplings relative to the one-bond couplings.
7. Report the couplings and the figure.

## Inputs (per variant)

| File | Content |
|---|---|
| `input/fid.npy` | Averaged raw FID: 65516 samples at 4000 Hz (16.4 s), ADC units |
| `input/acquisition.md` | Instrument, sequence (sudden field drop), sampling rate, known artifacts, instrument lines |
| `input/structure.md` | Name, SMILES, site and proton-group labels, couplings to report, output schema |

## Hidden reference (per variant)

| File | Content |
|---|---|
| `reference/couplings.json` | Reference value and tolerance of each graded coupling; which signs are graded |
| `reference/peaks.json` | Positions of the strongest resolved peaks of the processed spectrum (fixed recipe, above 10 noise sigma, instrument lines excluded) |
| `reference/zulf_sim.py` | Standalone zero-field line simulator used by the grader (exact diagonalization) |

## Evaluation (deterministic, same logic for every variant)

Score in [0, 1]. A missing or malformed `results.json` scores 0.

- Tier 1, one-bond couplings (0.4): each one-bond 13C-1H coupling of the variant within its tolerance; the 0.4
  is shared equally among them. Sites that are indistinguishable from the data alone are compared as an
  unordered set.
- Tier 2, signs and long-range couplings (0.3): for each long-range coupling listed in
  `reference/couplings.json`, the correct sign relative to 1J and the value within its tolerance; shared
  equally. Variants whose reference does not establish a sign grade the value only.
- Tier 3, spectral consistency (0.3): the grader simulates the line lists of the variant's 13C isotopologues
  from the reported couplings with its own simulator and checks the fraction of hidden experimental peaks that
  have a simulated line (relative amplitude >= 0.05) within 0.15 Hz; 0.3 x min(1, fraction / 0.8).

Fixtures per variant: the reference fit (scores 1.0) and a textbook guess (typical 1J values, positive 2J;
scores below 0.3).

## Variants

Reference values are best complex fits of each FID (ZULF_Model, github.com/gao-xh/ZULF_Model, private). Graded
items and tolerances are set to at least twice the spread over refits and processing variants (to be measured
before submission, see open points).

| Variant | Molecule | Isotopologues | Graded one-bond couplings (Hz) | Graded long-range |
|---|---|---|---|---|
| alanine | L-alanine | 13CH3, 13CH | 1J(CH3) 129.87, 1J(CH) 145.39 | to be measured |
| lactic_acid | L-lactic acid | 13CH3, 13CH | 1J(CH3) 129.06-129.15, 1J(CH) 146.79-147.04 | to be measured |
| pyridine | pyridine | 13C2, 13C3, 13C4 | 1J(C2) 178.68; 1J(C3), 1J(C4) 163.08 / 163.24 (unordered) | to be measured |
| triethylamine | triethylamine | 13CH2, 13CH3 | 1J(CH2) 130.69-131.02, 1J(CH3) 124.98-125.02 | 2J(C_CH2,H_CH3) -4.86 and 2J(C_CH3,H_CH2) -3.07 (sign + 1 Hz), 3J(H,H) 7.17 +- 0.5 |

Two further measured samples (ethylenediamine, N-ethylmethylamine) are left out: their best models are still
open (an unidentified second species; a coupling that moves between fit basins), so their reference values are
not yet reliable.

## Why this is hard

- No library performs zero-field NMR analysis: the agent must derive the Hamiltonian, the sudden-drop
  observable and the isotopologue set, and implement them correctly.
- Real data with real artifacts: detector saturation and ringing after the switch, baseline drift, a
  first-order phase from the switching delay, mains harmonics and instrument lines. Processing changes the
  lineshapes, so the model must pass through the same processing as the data.
- Overlapping multiplets of several isotopologues near 1J and 1.5 x 1J; small couplings and their signs only
  come out of a joint fit with many local minima.
- An expert with a dedicated fitting code needs days per molecule; the reference analyses took several
  sessions each.

## Notes for the ALE team

- Data: unpublished NMRduino measurements (4 kHz, sudden-drop sequence) of <lab, institution>. Release
  under CC-BY-4.0 needs the PI's approval (open point 1).
- The grader's simulator and peak lists will be provided as standalone files; no proprietary software.

---

## Open points before submission (not for the form)

1. Data release approval for the four FIDs (CC-BY-4.0, public).
2. Literature leakage: check whether zero-field couplings of these molecules are published (e.g. Andrews et
   al., arXiv 2604.26071, 13 molecules at natural abundance; alanine and lactic acid are classic ZULF examples).
   Drop or down-weight variants whose values can be looked up; Tier 3 stays data-bound in any case.
3. Robustness of each reference: refit from independent starts and with several processing recipes; record
   the spread of every graded coupling; set tolerances to at least twice it; grade a sign only where a
   sign-flip refit is clearly worse (done for triethylamine only: 14x).
4. Build `peaks.json` per variant from the FID with one fixed recipe and check that the reference passes Tier 3
   and the textbook fixture fails.
5. Write `reference/zulf_sim.py` without importing ZULF_Model and test it against
   `zulf_core.physics.compute_transitions` on every isotopologue.
6. Fill in the lab and institution in "Notes for the ALE team".
