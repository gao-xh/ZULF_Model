---
name: zulf-fid-processing
description: Diagnose and process an averaged experimental ZULF FID before any interpretation - acquisition assumptions, first-point/plateau/ringing diagnostics, crop and SG baseline choice, instrument-line and background screening against reference datasets, choice of time window, and the residual scale (no-lines residual and noise floor) that every later fit must be reported against.
---

# Processing an experimental ZULF FID

All processing goes through `zulf_core.render.acquisition` (`Acquisition`,
`process_record`, `evaluate_spectrum`); never re-implement it.

## Acquisition assumptions

- Without an ini, match the array length to earlier datasets of the same
  instrument (NMRduino: 65536 acquired, 65516 after averaging, 4000 Hz) and
  state the assumption in the report.

## Diagnostics

- `zulf_core.diagnostics.diagnose_fid(fid, fs)`: first-point anomaly,
  saturation plateau end, ringing end, baseline exponentials, late noise RMS.
- Crop after the ringing end (0.05 s worked; 0.1 s is safe for strong data).
- SG window from the baseline rates: a fast component (about 20/s) needs a
  short window. SG 801 at 4 kHz left a broad background and an 11 Hz ripple;
  SG 201 (order 2) worked on every sample so far. Mean removal on.
- Frequency grids that are not on the zero-filled FFT grid fall back to a
  slow direct DTFT; build grids as integer multiples of fs / n / zero_fill.

## Instrument lines and background

- Process one or more reference datasets of the same instrument identically;
  a line present in every dataset is instrumental. Known here: 60 Hz and
  harmonics, 294.0, 322.9, 442.9, 922.9 Hz. Exclude them from fit ranges with
  small gaps (for example ranges 62-119, 121-179, 181-293, 295-320 Hz).
- Below about 100 Hz there is baseline leakage and structured background, and
  a slow 1.6-3 Hz decaying oscillation after the baseline exponentials in
  every dataset. Compare any low-frequency feature with a no-signal dataset
  before interpreting it.
- Show the magnitude spectrum on a log scale (0-1500 Hz) and a linear
  zoom with instrument lines blanked, plus windowed versions.

## Time window

- Signals that decay within about a second are best fitted on a truncated
  window (for example samples 200-4200, 0.05-1.05 s) with zero_fill=4: the
  rest of the record is noise. The solver models the finite record exactly.
- Show full, about 1 s and about 0.3 s windows; the 0.3 s window has the
  best SNR for broad features, the full record the best resolution.

## Truncation and apodization

- Shortening the record lowers the noise per spectral point and costs
  resolution; with an exact model and a noise-scaled chi-square the late,
  noise-only samples barely change the estimates, so the window is chiefly a
  robustness and diagnostic knob.
- Choose an exponential apodization rate from the data: scan the peak SNR of
  the strongest lines versus rate at a long window. Rates taken from a fit can
  be inflated by unresolved fine structure (e3d282da: fitted 3-5 1/s, SNR best
  at 0-0.3 1/s); over-apodization merges close lines.
- Scan the window length for the leading model: parameters that drift with the
  window indicate a model problem at late times; fitted decay rates that fall
  with longer windows indicate unresolved structure.
- Use short windows (about 0.3 s) for detection, masks and model-free phasing;
  long windows (1-4 s) for small couplings.

## Residual scale

- `scripts/window_references.py FID --start 200 --stop 4200 --ranges ...`
  prints the no-lines residual (per-band linear background only) and the noise
  floor (residual expected from noise, from signal-free points) for a window
  and range, for complex or phased real-only data.
- Report every fit as (no-lines - r) / (no-lines - floor). The floor is rough
  where the range contains structured background.
