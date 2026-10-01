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

Commit: see git log (this file).
