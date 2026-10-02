# Blake pyridine series: model-free consistency of the spectra (2026-10-01)

- Data: Blake `6_3_26_pyridine_{02,33,50,66,75,100}` processed real spectra
  (frequency and amplitude arrays only; no acquisition or processing metadata
  in the upload), x = pyridine mole fraction 0.02 ... 1.00 (not committed).
- Question: every joint-fit solution so far (with and without the
  missing-peak penalty) moves many small couplings by several hertz in one
  step between x 0.66 and 0.75. Is there something in the spectra themselves
  at that point (acquisition, processing), independent of any model?
- Code: scratchpad data_check.py, interp_check.py, wide_check.py (data only,
  no fitting).

## Checks and results

Frequency grid: identical for all six (0 ... 4000 Hz, 88458 points, step
0.04522 Hz). Noise / maximum (5-point high-pass outside 135-205 Hz): 0.0092,
0.0031, 0.0012, 0.0024, 0.0027, 0.0028.

Main peak (data maximum in 140-200 Hz):

| x | 0.02 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|
| position (Hz) | 168.71 | 168.35 | 167.95 | 167.58 | 167.36 | 166.82 |
| slope from previous (Hz per unit x) | | 1.2 | 2.4 | 2.3 | 2.4 | 2.2 |
| FWHM (Hz) | 1.15 | 1.10 | 1.00 | 0.96 | 0.96 | 0.88 |

Peak tops (prominence >= 12 % of the maximum), offset from the main peak (Hz):

| peak | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|
| ~172 Hz | +3.85 | +3.93 | +3.89 | +3.84 | +3.84 |
| ~165 Hz | -3.21 | -3.21 | -3.30 | -3.21 | -3.31 |
| ~156 Hz | -11.89 | -11.90 | -11.98 | -12.21 | -12.30 |
| ~176 Hz | +7.74 | +7.68 | +7.69 | +7.82 | +8.14 |
| ~186 Hz | +17.77 | +17.72 | +17.78 | +17.95 | +18.35 |

Splitting of the ~159 Hz doublet 0.59-0.64 Hz throughout (shallower dip from
x 0.75). Small peak at +1.5 Hz: height 0.35, 0.31, 0.20, 0.15, absent at
x 0.33 ... 1.00 (smooth decrease). Peak positions show no step at 0.66 / 0.75.

Spectra aligned on the main peak (offset grid, normalised to the maximum):
change between neighbours (relative norm of the difference) per unit x:
1.09, 1.02, 1.02, **1.80** (0.66 -> 0.75), 0.80. Mean level per band
(data / max):

| band (Hz from main) | 0.02 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|
| +10 ... +14 | 0.204 | 0.196 | 0.203 | **0.120** | 0.192 | 0.188 |
| other bands | smooth in x | | | | | |

46 % of the squared 0.66 -> 0.75 difference lies in the +10 ... +14 Hz band
(177.6-181.6 Hz at x 0.66). There the x 0.66 spectrum lacks most of the
broad rising feature that every other spectrum has (figure
runs/processed/joint_peakpen_smooth_ms48/data_066_band.png); the sharp tops at
181.5 and 185.4 Hz are present. A feature near 294 Hz (outside the fit range)
decreases with x (0.29, 0.35, 0.12, 0.10, 0.04, 0.06): probably a solvent
line, not specific to x 0.66.

Figures: runs/processed/joint_peakpen_smooth_ms48/data_aligned.png,
data_066_band.png.

## Conclusion

Peak positions and spacings change smoothly with x; nothing in them steps at
0.66 / 0.75. The x 0.66 spectrum is an outlier in one region: the broad
feature at +10 ... +14 Hz from the main peak is about 40 % lower than in all
other spectra, including both neighbours. With monotone couplings the joint
fit cannot follow a feature that drops at 0.66 and returns at 0.75; this is a
likely cause of the one-step jumps of the small couplings between 0.66 and
0.75. Without metadata the origin (acquisition, sample, processing) cannot be
told from the arrays.

## Open points

- Joint fit with x 0.66 left out (`--leave-out 3`): do the jumps disappear or
  move to 0.50 / 0.75?
- Ask for the acquisition / processing records of the x 0.66 spectrum.

## Follow-up (2026-10-02): raw FIDs

The user supplied the raw FIDs as CSV (row 1 time, row 2 signal; 65515 points,
dt 2.5000382e-4 s, 16.4 s) for x 0.02, 0.10 (new), 0.33, 0.50, 0.66, 0.75,
1.00 (uploads, not committed). All seven: first 10 points pinned near 28840
(saturation), switching edge 3.18-3.25 ms, ringing ends 44-52 ms
(zulf_processing.diagnose_raw); nothing sets x 0.66 apart in the acquisition.

Same processing for all (zulf_processing.process_dataset; crop 200-16200,
SG 201/2, apodization 0.6 1/s, zero fill 2) and the instrument phase
calibration (configs/confirmed_samples.json: phase0 176.3 deg, delay = edge -
0.033 ms; NMRduino 4000 Hz sequence). Scratchpad fid_load.py, fid_compare.py,
fid_calphase.py. Figure:
runs/processed/joint_peakpen_smooth_ms48/fid_calphase_vs_blake.png.

Band +10 ... +14 Hz from the main peak, mean level (data / max):

| x | 0.02 | 0.10 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|---|
| from FID, magnitude, 2 s window | 0.404 | 0.373 | 0.302 | 0.267 | 0.266 | 0.279 | 0.309 |
| from FID, calibrated phase, real | -0.065 | 0.213 | 0.121 | 0.136 | 0.113 | 0.132 | 0.105 |
| Blake processed (real) | 0.204 | - | 0.196 | 0.203 | 0.120 | 0.192 | 0.188 |

From the FIDs, x 0.66 is not an outlier in this band. Blake processed minus
our calibrated real spectrum (both normalised to the maximum, 140-200 Hz): rms
0.17-0.29, of which 0.17-0.28 is smooth on a 2 Hz scale and only 0.020-0.036
is not. Sharp lines agree; Blake's spectra carry no broad negative regions,
ours do (down to -0.4 at 160-166 and 172-178 Hz). So the processed spectra
are, to within 0.02-0.04, our phased spectra plus a smooth baseline, and that
baseline correction differs between spectra: in this band it raised the
level by about +0.06 to +0.08 in the other spectra (0.12-0.14 -> 0.19-0.20)
and by about +0.01 in x 0.66 (0.113 -> 0.120).

Conclusion: the x 0.66 band deficit comes from the baseline correction of the
processed spectrum, not from the measurement. Every joint fit so far used
these baseline-corrected spectra; the small-coupling jumps between 0.66 and
0.75 are probably fitted to that processing difference.

Open points: fit the series from the FIDs with one processing for all and the
solver's own background terms (complex or calibrated real), x 0.10 added;
check whether the broad negative regions are signal (fast-relaxing component)
or a processing effect.

## AsLS baseline, the same for all seven (2026-10-02)

zulf_processing.asls_baseline (new, commit 4da54bb; lam = (smooth_hz / step)^4)
on the calibrated real spectra above, p 0.01, smooth 1.5 and 2 Hz (scratchpad
fid_asls.py; figure runs/processed/joint_peakpen_smooth_ms48/fid_asls_1.5_vs_blake.png).

| x | 0.02 | 0.10 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|---|
| band +10 ... +14 Hz, smooth 1.5 Hz | 0.287 | 0.368 | 0.306 | 0.291 | 0.278 | 0.289 | 0.282 |
| band +10 ... +14 Hz, smooth 2 Hz | 0.291 | 0.382 | 0.298 | 0.288 | 0.274 | 0.289 | 0.289 |
| Blake - ours rms (1.5 Hz) | 0.133 | - | 0.080 | 0.075 | 0.088 | 0.091 | 0.107 |

With one baseline rule for all, the band is smooth in x (0.27-0.31 from
x 0.33 on). Remaining differences to Blake's spectra: our valleys between
close lines are shallower (e.g. between the main peak and the +3.9 Hz peak
about 0.3 against 0.05) and our lines slightly wider; x 0.02 keeps a broad
hump at 143-150 Hz (lowest signal).

## Phase from these data only (2026-10-02)

Requested: no earlier calibration values; FID cropped to its first part, no
apodization. Crop 200 ... 200 + 1 s (and 2 s), apodization 0, zero fill 4.

- zulf_processing `lines` criterion (complex Lorentzian fits of resolved lines):
  not usable here; 5-11 lines per spectrum, mostly overlapping, delays ran to
  the search bounds (+-0.5 ms near the edge, +-6 ms free, -1.5 ms joint).
- Peak tops (scratchpad apex_phase.py): first-order phase from each FID's own
  switching edge (3.18-3.25 ms), then the phase of the complex spectrum at every
  magnitude peak top (prominence >= 15 %, 145-195 Hz, 180 Hz excluded), pooled
  over the seven spectra (weights height^2). Zero-order phase 175.7 deg (1 s,
  43 tops) and 169.3 deg (2 s, 55 tops), coherence 0.676 / 0.644. An extra
  delay on top of the edge does not help (best at the +2 ms search bound,
  coherence 0.673 / 0.645). This agrees with the instrument calibration
  (176.3 deg, edge - 0.033 ms) without using it.
- Systematic residual phases at the same peaks in every spectrum: main peak
  (166.5-168.5 Hz) +23 to +50 deg, peak near 171-172.7 Hz -38 to -69 deg. One
  zero- and first-order phase cannot make all lines absorptive.

Figure runs/processed/joint_peakpen_smooth_ms48/front1s_phased_vs_blake.png
(phase0 175.7 deg, delay = own edge; real part without baseline and with AsLS
1.5 Hz).

Full record, only the switching step removed (user request; no apodization,
zero fill 1, plan built without the ringing rule): crop from the first sample
after the edge (14-15) to 65515. The 3-50 ms switching ringing then dominates:
magnitude 1-2 across 130-210 Hz against about 0.03 for the pyridine lines
(S/N 25-35, only ringing tops at 145.0 and 180.8 Hz). From 210 (52.5 ms, after
the ringing in all seven) to 65515: lines clear, S/N 8 (x 0.02) to 59
(x 0.50); peak-top phase 177.6 deg with the delay at each FID's own edge
(coherence 0.664; an extra delay does not help); the same residual pattern
(main peak +37 to +57 deg, ~171-172 Hz peak -26 to -77 deg). Figures
runs/processed/joint_peakpen_smooth_ms48/crop_ringing_effect.png,
full_phased_vs_blake.png.

Commit: see git log (this file).
