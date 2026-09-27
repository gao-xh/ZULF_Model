---
name: zulf-phasing
description: Phase experimental ZULF spectra correctly and without bias - project conventions, when phasing is unnecessary (complex fits), model-free phasing from sign-agnostic peak coherence, delay/phase ambiguity in narrow bands, and the ripple that first-order correction creates from record-start baseline. Use before showing or fitting a phased (absorption) spectrum.
---

# Phasing ZULF spectra

Details and pitfalls: `references/phasing.md`.

## Per dataset (D44)

The phase comes from the instrument and the processing, never from the
molecule. Default (`process_dataset(..., phase_criterion="calibration")`,
used by `scripts/analyze_sample.py`): delay = this dataset's switching edge
(raw FID, half height of the field switch-off) + a calibrated offset, zero-
order phase from the instrument calibration (`calibrated_phase`, config
`processing.phase_calibration`). Leave-one-out on the confirmed samples it is
within 14 deg of each sample's own complex fit at 125-250 Hz.
Refresh the calibration (median of fitted delay - edge, circular mean of
fitted phase0 carried to that delay) whenever the pulse sequence, hardware or
processing defaults change; leave out fits whose delay is > 1 ms from the
others (harmonic-band ambiguity).

Model-free alternative (`phase_criterion="entropy"` or `"lines"`,
`zulf_processing.phase_dataset`): global search over phase0 (0-180 deg) and
delay (+-0.5 ms around the edge), then Nelder-Mead fine-tune;
`register_phase_criterion` adds criteria. On J-spectra these err 2-86 deg
(lines of both signs overlap) and pull the delay to the window bound; never
use one without running `scripts/validate_processing.py` (about 7 s) and
reporting the method next to the phased figure. Without the edge, bands
near harmonics (J, 2J) make delays about 4 ms apart equivalent.

## Rules

1. Prefer complex fitting (real and imaginary) with zero-order phase and
   delay fitted; no phasing is needed and nothing is assumed.
2. A phase taken from a hypothesis fit must not be used to test that
   hypothesis: rank hypotheses on complex fits. For an unbiased phased
   spectrum a model-free phase is needed, but it only works when lines are
   resolved: `python scripts/phase_coherence.py FID.npy --start 200 --stop
   4200` or `zulf_processing.phase_dataset(..., criterion="lines")`. On dense J multiplets
   (lines of either sign closer than one width) every model-free estimator
   tried failed (D43; e66a4b08: minimum-entropy, doubled-angle and per-line
   fits, delay off by 0.4-4.7 ms, 10-40 deg line errors even noise-free).
   Then show the data phased with each hypothesis's own complex fit
   (`fit_structure(..., route="both")`, `report.fit_phasing`) and say that
   the phase came from that fit; the phased real-only refit is a display and
   cross-check, not independent evidence.
3. Check the coherence curve: with all strong peaks in one narrow band the
   delay is weakly determined (plateau of several ms); say so, and show what
   changes outside the band (for example the shape of a weak line at 2J).
4. ZULF lines can be negative. Phase criteria must be sign-agnostic (doubled
   angle), and the overall sign is a convention (strongest peak positive).
5. A ripple of period 1/(crop + delay) after first-order correction is the
   record-start baseline, not signal. Confirm with a no-signal dataset or by
   changing the crop. In solver fits of phased spectra keep background_order
   >= 1; the background is modelled in the record frame.
6. Real-only fits use half the information and depend on the fixed delay; use
   them as a cross-check on the signal bands only.

## Conventions

`phase_correct(values, f, phase0, delay, acq)` multiplies by
exp(-i (phase0 + 2 pi f (delay + crop reference))). The solver's shared gain
phase and `phase_delay` use the same convention, so a complex fit's phase can
be applied directly, and `ObservedSpectrum.from_spectrum(..., phasing=...)`
applies it to the model.
