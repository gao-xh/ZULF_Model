---
name: zulf-fid-processing
description: Diagnose and process an averaged experimental ZULF FID before any interpretation - acquisition assumptions, first-point/plateau/ringing diagnostics, crop and SG baseline choice, instrument-line and background screening against reference datasets, choice of time window, and the residual scale (no-lines residual and noise floor) that every later fit must be reported against.
---

# Processing an experimental ZULF FID

All processing goes through `zulf_core.render.acquisition` (`Acquisition`,
`process_record`, `evaluate_spectrum`); never re-implement it.

## Package zulf_processing (D44): use it for every dataset

- `process_dataset(fid, fs, sample_id)` gives the per-dataset diagnostics,
  a `ProcessingPlan` (each choice with its reason), the complex observation
  for fits, the phase and the phased observation; `.figure(path)` and
  `.save_record(path)` document it. Processing parameters may differ between
  datasets (crop after this dataset's ringing, later SG and window from the
  data), so never copy another dataset's phase: the crop reference and the
  line shapes change with the plan.
- The raw FID starts with the field switch-off: plateau (2.75 ms), a steep
  edge (3.41-3.51 ms on the NMRduino no-dead-time sequence), ringing to about
  50 ms. `switching_edge` measures the edge per dataset; minus its time is
  the delay of the processed spectrum.
- After a change to the package run
  `ZULF_DATA_DIR=... python scripts/validate_processing.py` and compare with
  the last summary in ANALYSIS_LOG (phase errors against the known-structure
  fits, next to the leave-one-out median calibration). It takes about 7 s;
  if it is slow, something in diagnostics or phase regressed (the unused
  exponential fits of zulf_core.diagnose_fid once took 546 s per dataset).
- Phase default: instrument calibration + this dataset's edge
  (`phase_criterion="calibration"`, config `processing.phase_calibration`);
  see the zulf-phasing skill.

## Tune the parameters for every new FID set (required)

Defaults are starting points, not answers: re-derive crop, window, phase and
baseline for each new set of FIDs and record the choice and its evidence in
the analysis log (lesson of the Blake pyridine series, 2026-10-02,
docs/analysis/2026-10-01_blake-pyridine_data-consistency.md).

1. Sampling rate from the data (CSV time column, ini), not 4000 Hz by habit;
   check the frequency axis of any supplied processed spectrum against your
   own (the Blake spectra sat 0.1-0.2 Hz high).
2. Crop start: after the ringing (diagnostics) is only the lower bound. Scan
   the start (0.05, 0.1, 0.2, 0.3, 0.5 s) and plot the main-peak FWHM and the
   depth of a valley between close lines; take the earliest start where both
   stop changing. Pyridine: 52.5 ms left a fast-decaying broad component
   (valleys 0.2-0.3 of the maximum, broad negative regions); 0.1-0.2 s gave a
   plateau (FWHM 1.11 / 1.01 / 0.98 Hz at x 0.33 / 0.50 / 0.75, valleys ~0)
   and matched the supplied processed spectra. Later starts keep narrowing the
   lines when the decay is not single-exponential: say so instead of choosing
   the narrowest.
3. Do not start right after the switching step: the ringing (497 Hz, 3-50 ms,
   1e4 counts) swamps the spectrum (magnitude 50x the lines).
4. SG drift removal runs on the full record before the crop
   (process_record). Check that the broad features do not depend on the SG
   window/order or on SG-before-crop vs after; on pyridine they did not.
5. Window: compare apodization 0, 0.3, 0.6, 1 1/s (noise/max and line
   width); 0.3 1/s cut the noise 8x for +0.1 Hz width.
6. Phase: see the zulf-phasing skill (edge delay, data-only zero-order phase
   by peak symmetry with a model-spectrum bias check, one phase per set unless
   the per-spectrum values differ beyond their errors).
7. Baseline: do not remove what the processing itself creates. A crop that
   starts at t_c after the switching edge, corrected back to the edge, gives
   every positive line oscillating wings (period 1/t_c, first lobe about -0.12
   of the height at 52.5 ms and -0.21 at 0.1 s; one-line test, SG irrelevant).
   On pyridine these wings looked like broad negative regions or negative
   lines; with sudden drop and gamma-weighted preparation and detection all
   line amplitudes are positive. Positive-peak AsLS (and the supplied
   processed spectra) flattened the wings. Fit such spectra with the model
   rendered through the same processing: give `record` (the Acquisition dict)
   and `phasing` ({phase0_rad, delay_s}) per series entry
   (fit_joint_series / reliability_series pass them to from_spectrum); keep
   the record short (signal extent, e.g. 4 s: jacobian 1.8 s instead of 12 s
   for 16 s). One rule for all spectra of a series: per-spectrum hand
   baselines created the x 0.66 outlier band.
8. Peak-top noise floor for missing-peak rows (`--peak-min-sigma`): raise it
   for unwindowed or weak spectra (x 0.02 had 49 noise tops at 2 sigma).

## Averaging scans of a run

- `scripts/average_scans.py RUN OUT [--exclude-z Z]`: decodes every `<n>.dat`, a per-scan deviation from the
  mean (0.1-4.1 s, own mean and trend removed) with a robust z, even / odd half averages and scans.json.
  Look at OUT/scans.png before choosing Z. Acetonitrile 2026-09-28: 11 % of the scans had z > 5 (scattered
  disturbances); leaving them out changed no fitted number, so check this before trusting a selection.
- The half averages are the held-out check: fit both with the same settings and compare couplings and
  nuisance parameters (field, delay).
- Use the run's own sampling rate (`--sampling-rate` in make_series_entry.py): the 2 kHz sequence
  (standard_zf_2000Hz_no_dead.seq) exists next to the 4 kHz one, and the config's phase calibration is for
  4 kHz only; a complex fit absorbs the phase, a phased figure does not.

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
- Fit the whole spectrum, not only the bands where the data have lines
  (`make_series_entry.py` default since 2026-10-08): a model line in an
  unfitted gap is not penalised, and with rate families its rate can run to
  the lower bound. Ethanol 2026-09-13: ranges 198-221 and 244-260 Hz let the
  model put tall sharp lines at 232-236 Hz where the data are empty. Check the
  whole-spectrum figure (scripts/paper_figure.py) for model lines outside the
  fit ranges. The 83 Hz dispersive feature of the NMRduino runs (acetonitrile,
  ethanol, isopropanol) is excluded with `--exclude 81.5,86`.
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
- Start of the record (2026-10-10, acetonitrile LF_1): right after the field
  switch there can be a transient that the SG drift filter does not remove. It
  shows as broad troughs or humps under every band, scales with the signal and
  grows with the SG window. The model has no such term; on high-SNR data a fit
  bends other parameters to imitate it (LF_1: a 117 nT field instead of 64 nT).
  Check: reprocess with crop 0.3 s; if the broad structure goes away, fit with
  `make_series_entry.py --crop 0.3` (the solver models the crop exactly).

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
- Default for refinement: a window of 1 s or more with weak matched
  apodization chosen from the SNR curve (e3d282da: 0.6 1/s gave identical
  parameters at 1, 2 and 4 s and a lower residual than hard truncation).
- An isotopologue amplitude ratio that drifts with the window (e3d282da:
  Cb/Ca 1.1 at 0.3 s, 0.45 at >= 1 s) is a model-deficiency signal: the
  short window is right about abundance, the long window exposes missing
  couplings.

## Residual scale

- `scripts/window_references.py FID --start 200 --stop 4200 --ranges ...`
  prints the no-lines residual (per-band linear background only) and the noise
  floor (residual expected from noise, from signal-free points) for a window
  and range, for complex or phased real-only data.
- Report every fit as (no-lines - r) / (no-lines - floor). The floor is rough
  where the range contains structured background.
