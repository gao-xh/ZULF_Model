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

Blake's processed x 0.10 spectrum (uploaded 2026-10-02, same grid): main peak
168.62 Hz, band +10 ... +14 Hz 0.213, in line with x 0.02 (0.204) and 0.33
(0.196); x 0.66 (0.120) stays the only outlier among the seven processed
spectra. Added to full_phased_vs_blake.png.

## Phase from peak symmetry, checked on a simulated spectrum (2026-10-02)

User suggestion: use the symmetry of the peaks. Absorption is symmetric about
the line centre, dispersion antisymmetric. Centre from the magnitude top
(parabola, phase independent); asymmetry of the real part over +-0.4 Hz:
sum (R - L)^2 / sum (R + L)^2; per-peak phase = least asymmetric (absorption
sign); global phase0 = minimum of the height^2-weighted sum (scratchpad
sym_phase.py; spectra 52.5 ms to the end, no window, first-order phase from
each FID's own edge).

Data: global phase0 174.0 deg (extra delay does not help: +0.40 ms gives the
same mean asymmetry 0.371). Same check on a noise-free model spectrum
(sym_model.py: couplings of the penalty-run best at that x, all rates 2.3 1/s,
natural abundance ratios, no phase error, same acquisition and first-order
step): global 353-357 deg, i.e. the criterion is biased by about -5 deg. So
the data phase is about 179-180 deg with the delay at the switching edge.

Per-peak phase minus global (deg):

| x | source | ~151 | ~159.5 | main | ~172 | others |
|---|---|---|---|---|---|---|
| 0.33 | model | +5 | +21 | +9 | -11 | 178.6: +92 |
| 0.33 | data | +40 | -12 | +62 | -16 | 164.8: +109, 176.7: -178 |
| 0.50 | model | +15 | +17 | +10 | -11 | 178.2: +95 |
| 0.50 | data | +29 | -17 | +42 | -47 | 164.4: +95, 176.9: +147 |
| 0.75 | model | +8 | -50 | +26 | +5 | |
| 0.75 | data | +14 | -58 | +54 | -7 | 176.4: +156 |
| 1.00 | model | -37 | +12 | +8 | -1 | 177.4: +107 |
| 1.00 | data | | -47 | +44 | -78 | |

Overlap alone shifts per-peak symmetric phases by up to about 25 deg on the
main lines (and 90-110 deg at a model feature near 178 Hz). The data exceed
this: the main peak +42 to +62 deg (model +8 to +26), the ~172 Hz peak -47 and
-78 deg at x 0.50 and 1.00 (model -11, -1); the data features near 164.5 Hz
(+95 to +109) and 176.5 Hz (about 180 deg, a negative line or a dispersive
pair) have no counterpart among the model's peak tops.

Complex fits with couplings held (complex_fit.py; one phase vs one phase per
isotopologue; delay, rates, gains and a linear background free): unstable.
Delays wander from -3.1 to +3.9 ms, the 13C2 component takes rate 9-13 1/s
and up to 7 times the 13C3 amplitude (it absorbs a broad component), relative
residuals 0.27-0.46; per-isotopologue phases improve the residual by 2-10 %
only. Not usable as a phase reference until the model covers the broad part.

Conclusion: one global phase (about 180 deg at the switching edge) is
supported by the data and by the simulation check; the remaining per-peak
departures are partly overlap and partly real model-data differences (line
phases or line structure the sudden-drop model with these couplings does not
have).

Centre choice (user question: why not the real-part apex?). scratchpad
sym_centre.py: per trial phase the apex of the rotated real part (parabola,
within +-0.3 Hz of the magnitude top) as centre, against the magnitude top.

| x | model error, magnitude top | model error, apex | data, magnitude top (model-corrected) | data, apex (model-corrected) |
|---|---|---|---|---|
| 0.33 | -5 | +4 | 185 | 157 |
| 0.50 | -7 | 0 | 177 | 177 |
| 0.75 | -6 | -7 | 181 | 203 |
| 1.00 | -3 | -2 | 171 | 188 |

On the noise-free model both centres find the phase within 7 deg (mean
per-peak scatter 20-29 deg vs 20-30 deg). On the data the magnitude-top
centre gives a phase that agrees between spectra (171-185 deg) and the apex
centre does not (157-203 deg). Per-peak values depend on the centre choice
(main peak at x 0.50: +46 vs +14 deg), so only the global phase is robust.
Kept: magnitude-top centre.

Light apodization (2026-10-02; scratchpad apod_plot.py): exponential 0, 0.3,
0.6, 1.0 1/s (+0, 0.10, 0.19, 0.32 Hz FWHM) on the 52.5 ms - end crop, phase
179 deg, AsLS. Noise / maximum at x 0.50: 0.008 without, 0.001 with 0.3 1/s.
The valley between the main peak and the ~172 Hz peak stays at 0.18-0.29 for
every apodization (Blake about 0). Figures
runs/processed/joint_peakpen_smooth_ms48/apod06_phased_vs_blake.png,
apod_series_x050.png.

Frequency axis: cross-correlation (0.6 1/s spectra against Blake's) needs
Blake's spectra shifted down by 0.105-0.160 Hz over 145-195 Hz (145-170 Hz:
0.095-0.155; 170-195 Hz: 0.12-0.235). Our axis uses the sampling rate from the
CSV time column (3999.939 Hz); a 4000 Hz axis would shift 168 Hz by only
+0.003 Hz. Blake's lines are therefore about 0.1-0.2 Hz higher than in the raw
data, for reasons not visible in the arrays; absolute line positions (and so
the large couplings) of fits to the processed spectra carry this offset.

Phase refinement on the 0.3 1/s spectra (2026-10-02; scratchpad phase_opt.py):
per spectrum the symmetry phase (magnitude-top centre, noise from point
differences, prominence >= 15 % and >= 8 sigma: 5-12 peaks), a jackknife
error over peaks, and the criterion bias from a model spectrum at the same x
(same acquisition and apodization, all rates 3 1/s, no phase error).

| x | 0.02 | 0.10 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|---|
| symmetry phase (deg) | 177.5 | 172.5 | 171.5 | 169.5 | 172.5 | 170.0 | 167.0 |
| jackknife SE (deg) | 17.6 | 39.2 | 19.6 | 20.1 | 25.3 | 29.2 | 34.6 |
| model bias (deg) | +10.0 | +19.5 | -7.0 | -8.0 | -14.5 | -7.5 | -6.5 |
| corrected (deg) | 167.5 | 153.0 | 178.5 | 177.5 | 187.0 | 177.5 | 173.5 |

Weighted mean 174.8 deg, chi2 0.8 for 6 dof: the spread between spectra is
within the errors, no sign of spectrum-specific phases. Error of the mean
about 9 deg; the per-peak spread (overlap) limits this criterion. Joint scan
over all spectra: delay at the edge gives phase0 173.5 deg (cost 3.885); an
extra +0.25 ms gives 158.0 deg at cost 3.882, i.e. phase0 and delay trade off
on a flat surface. The joint fit refines a common phase per spectrum anyway
(shared_phase gain model), so a residual of this size is absorbed. Next fit
input: 0.3 1/s spectra phased at 174.8 deg, AsLS (series_fid03.json in the
scratchpad).

Commit: see git log (this file).
