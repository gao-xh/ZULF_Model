# ALE-V2 task draft: J-coupling network of triethylamine from a raw zero-field NMR FID

Draft for the Agents' Last Exam (ALE-V2) chemistry subdomain, written in the fields of an ALE task card
(`task_card.json`: title, summary, software, taskPrompt, agentMustDo, inputFiles, referenceFiles, evaluation,
vm). Submission: https://agents-last-exam.org/submit (the portal turns a description plus raw data into a
runnable `main.py`). Status: draft, not submitted. Open points for the author are at the end.

## Title

Extract the J-coupling network of triethylamine from a raw zero-field NMR free-induction decay

## Summary

Given a raw, signal-averaged zero- to ultralow-field (ZULF) NMR free-induction decay of neat triethylamine at
natural 13C abundance and the molecular structure, the agent must build a zero-field spin simulation from first
principles, process the FID, fit the 13C isotopologue spectra, and report the one-bond and long-range
heteronuclear and homonuclear J couplings with their relative signs.

## Category

physical_sciences / chemistry (NMR spectroscopy, spin physics, spectral fitting)

## Software

- Python 3.10, NumPy, SciPy, Matplotlib (preinstalled on the VM)
- No NMR or spin-simulation packages; no network access needed

## Task prompt (visible to the agent)

You are working on a Linux VM.

### Your task
A zero-field NMR spectrometer recorded the free-induction decay (FID) of neat triethylamine, N(CH2CH3)3, at
natural isotopic abundance. Determine the J-coupling network of the molecule from this FID.

### Visible input
- `input/fid.npy`: averaged FID, float64, 65516 samples, sampling rate 4000 Hz (16.4 s). Sample 0 is the first
  sample recorded after the field switch command.
- `input/acquisition.md`: how the data were taken: prepolarization in a permanent magnet, shuttling into a
  magnetic shield, sudden switch-off of the guiding field, detection by an atomic magnetometer sensitive to the
  z magnetization; the initial magnetization and the detected signal of each nucleus are proportional to its
  gyromagnetic ratio; the detection electronics saturate for the first few milliseconds after the switch and
  ring at about 500 Hz for about 50 ms; mains interference at multiples of 60.06 Hz and an instrument line near
  180 Hz are present.
- `input/structure.md`: the structure (SMILES `CCN(CC)CC`) and the convention for naming couplings (below).

### What you must produce
Write `output/results.json`:

```json
{
  "J_hz": {
    "1J(C_CH2,H_CH2)": ..., "1J(C_CH3,H_CH3)": ...,
    "2J(C_CH2,H_CH3)": ..., "2J(C_CH3,H_CH2)": ...,
    "3J(H_CH2,H_CH3)": ...
  },
  "other_couplings_hz": { "...": ... },
  "isotopologues_used": ["..."],
  "notes": "..."
}
```

- All couplings in Hz with sign. Use the convention that one-bond 13C-1H couplings are positive; report every
  other coupling with its sign relative to them.
- `C_CH2` is the methylene carbon of one ethyl group, `C_CH3` the methyl carbon of the same group; `H_CH2` and
  `H_CH3` are the protons on them. Couplings to the protons of the other two ethyl groups go into
  `other_couplings_hz` (optional).
- Write `output/fit_figure.png` showing the processed experimental spectrum and your simulated spectrum over
  the signal bands.

### Constraints
- Do not modify files under `input/`. Write only under `output/`.
- Do not install packages and do not use external NMR software or web resources for coupling values; the
  values must come from your analysis of this FID.

## Agent must do

1. Inspect the raw FID: find the switching edge, the saturated and ringing samples, the baseline drift and the
   noise level; choose a crop start, record length, baseline removal and apodization.
2. Compute the spectrum and phase it (zero-order phase and the time delay between the field switch and the
   start of evolution).
3. Work out which natural-abundance isotopologues are visible (each 13C site with the protons of the molecule;
   N-H is absent; 14N, 15N and 13C-13C isotopologues negligible or not observable) and which proton groups are
   magnetically equivalent.
4. Implement a zero-field spin-Hamiltonian simulation (H = sum J_ij I_i . I_j) with the sudden-drop
   observable (gamma-weighted initial state and detection), giving line frequencies and amplitudes; exploit
   equivalence (collective spins) or brute force.
5. Fit the couplings, decay rates, amplitudes and phase of all visible isotopologues jointly to the processed
   spectrum, with the model passed through the same processing (crop and window) as the data.
6. Determine the signs of the two-bond couplings relative to the one-bond couplings.
7. Write `results.json` and the figure.

## Input files

| Name | Format | Path | Description |
|---|---|---|---|
| FID | NumPy array (.npy) | `input/fid.npy` | Averaged raw FID, 65516 samples at 4 kHz, ADC units |
| Acquisition notes | Markdown | `input/acquisition.md` | Instrument, sequence (sudden field drop), sampling, known artifacts |
| Structure | Markdown | `input/structure.md` | SMILES, coupling naming convention, output schema |

## Reference files (hidden)

| Name | Format | Path | Description |
|---|---|---|---|
| Reference couplings | JSON | `reference/couplings.json` | Reference values and tolerances of the graded couplings |
| Experimental peak list | JSON | `reference/peaks.json` | Positions (Hz) of the K strongest resolved peaks of the processed spectrum in 110-260 Hz, measured from the data with a fixed recipe |
| Grader simulator | Python | `reference/zulf_sim.py` | Independent zero-field line simulator used by the grader (exact diagonalization) |

Reference values (best complex fit of this FID, ZULF_Model, docs/analysis/2026-10-05_amines_overnight.md;
spread over the alternative fitted models in parentheses):

| Coupling | Reference (Hz) | Spread over models (Hz) | Tolerance (Hz) |
|---|---|---|---|
| 1J(C_CH2,H_CH2) | +130.69 | 130.69 - 131.02 | 0.6 |
| 1J(C_CH3,H_CH3) | +125.00 | 124.98 - 125.02 | 0.4 |
| 2J(C_CH2,H_CH3) | -4.86 | to be measured | sign + 1.0 |
| 2J(C_CH3,H_CH2) | -3.07 | to be measured | sign + 1.0 |
| 3J(H_CH2,H_CH3) | +7.17 | to be measured | 0.5 |

The signs of the two-bond couplings are supported by a refit with both signs flipped, which was 14 times worse
(objective 0.527 against 0.0383); the literature value for ethanol (12CH3-13CH2-OH) at zero field is
2J(C,H) = -4.6 Hz.

## Evaluation (deterministic)

Score in [0, 1], sum of three tiers; `results.json` missing or malformed scores 0.

- Tier 1 (0.4): the two one-bond couplings within tolerance (0.2 each).
- Tier 2 (0.3): the two two-bond couplings have the correct sign and lie within 1.0 Hz (0.1 each); the vicinal
  3J(H,H) within 0.5 Hz (0.1).
- Tier 3 (0.3): spectral consistency. The grader simulates the line list of the 13C-CH2 and 13C-CH3
  isotopologues from the agent's couplings with its own simulator (`reference/zulf_sim.py`) and checks that
  at least 80 % of the hidden experimental peaks (`reference/peaks.json`) have a simulated line with relative
  amplitude >= 0.05 within 0.15 Hz. This rewards a coupling set that reproduces the measured spectrum even
  where the individual small couplings are correlated.

Replay fixtures (to be built): `output_test_pos` (the reference fit, scores 1.0) and `output_test_neg` (a
textbook guess: 1J 135 / 125, 2J +4, 3J 7; scores below 0.3).

## VM

`cpu-free-ubuntu`, 4 vCPU, 16 GB, timeout 4 h (the reference analysis of one known-structure fit takes about
10-40 min of compute; most of the work is building and validating the simulation and the processing).

## Why this is hard (for the reviewers)

- No library does zero-field NMR: the agent must derive the spin Hamiltonian, the observable under a sudden
  field drop and the natural-abundance isotopologue set, and implement them correctly.
- Real data: saturation and ringing after the switch, baseline drift, a first-order phase from the switching
  delay, mains harmonics and an instrument line; processing choices change the lineshapes, so the model must go
  through the same processing.
- The spectrum is a superposition of two isotopologues with dense multiplets near 1J and 1.5 x 1J; small
  couplings and their signs are only reachable by a joint fit, and the fit has many local minima.
- Expert time: days (the reference analysis took several sessions of a dedicated fitting code).

## Open points before submission

1. Data release: ALE task data are published under CC-BY-4.0. The FID (4322bdfc-average_fid.npy) is
   unpublished lab data; the PI must agree to its release.
2. Literature leakage: check whether zero-field J values of triethylamine are published (e.g. the
   natural-abundance 13C ZULF set of Andrews et al., arXiv 2604.26071). If they are, the prompt's "no external
   values" rule is not enforceable; prefer a sample without published zero-field values, or grade mainly on
   Tier 3.
3. Reference robustness: the reference couplings are conditional results of one model (fast N-H exchange,
   rate families). Before submission, measure the spread of the graded couplings over refits from independent
   starts and processing variants, and set each tolerance to at least twice that spread.
4. Build `reference/peaks.json` from the FID with a fixed recipe (crop 0.1 s, 8 s record, 0.3 1/s, zero fill 3,
   calibrated phase; peak tops above 10 noise sigma in 110-260 Hz excluding 179.6-180.5 Hz) and check that the
   reference fit passes Tier 3 and the negative fixture fails it.
5. Write `reference/zulf_sim.py` as a standalone file (no ZULF_Model import) and test it against
   `zulf_core.physics.compute_transitions` on the two isotopologues.
