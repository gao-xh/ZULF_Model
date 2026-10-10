# acetonitrile fringe field LF_1 vs 250 uL: decay times per line (2026-10-10)

- Data (drive only): `original/2026-09-23-acetonitrile-fringe-field/ACN_fringe_field_ms1_LF_1/`, 2000 scans,
  2026-09-23 to 09-25, 2000 Hz; sample volume 2 mL (owner, 2026-10-10; not in the name), 8 x the 250 uL sample. z5 average (--exclude-z 5): 1744 scans,
  ~/research/zulf/data/processed/2026-09-23-acetonitrile-fringe-field/z5/.
- Question (Xuehan): fit it as the 250 uL data and compare the decay rates per line.
- Series runs/series/acn_lf1_full (same defaults as acn250_*_full); fit runs/processed/acn_lf1_full_fit with the
  command of docs/analysis/2026-10-10_acetonitrile-250ul_runs.md (uniform static field, one rate per model line,
  rates 0.1-40 1/s, 8 starts).
- Code: dc4baf2. Figure: runs/figures/acn_lf1_vs_250ul.png (J band 126-147 Hz and 2J band 262-285 Hz; top LF_1,
  bottom 250 uL both acquisitions; scripts/overlay_figure.py, display baseline "model").

## Results

| | LF_1 | 250 uL run 1 | 250 uL run 2 | 250 uL both |
|---|---|---|---|---|
| 1J(C1,H1) (Hz) | 136.248 | 136.282 | 136.289 | 136.286 |
| B transverse / z (nT) | 108.4 / 44.1 | 38.4 / 43.1 | 39.8 / 45.0 | 39.4 / 43.9 |
| \|B\| (nT) | 117.0 | 57.7 | 60.1 | 59.0 |
| score | 0.0152 | 0.233 | 0.214 | 0.165 |

Decay time T = 1 / rate (s) of the model lines (relative intensity > 0.05; FWHM = rate / pi):

| line | LF_1 | 250 uL run 1 | run 2 | both |
|---|---|---|---|---|
| J line 136.3 Hz (field-insensitive) | 2.63 (rel 0.05) | 2.33 | 2.22 | 2.26 |
| J sidebands 134.7 / 137.9 Hz | - (model puts them at 133.2 / 139.4, T 0.025 = bound) | 1.03 / 0.27 | 1.17 / 0.38 | 1.05 / 0.32 |
| 2J triplet 271.5-273.7 Hz | not in the model (see below) | 0.39-0.48 | 0.36-0.47 | 0.37-0.48 |
| 2J outer lines | 268.5 / 270.3 / 274.7 / 276.6 Hz: 0.11 | 270.6 / 274.6: 0.21 / 0.03 | 0.025 / 0.085 | 0.13 / 0.065 |

## Observations

- The field-insensitive J line is as narrow in LF_1 as in 250 uL (T 2.6 vs 2.2-2.3 s): the intrinsic decay is
  the same.
- The 2J band of LF_1 holds the same sharp triplet as 250 uL (271.5 / 272.5 / 273.6 Hz) plus two broad wings
  centred near 267.5 and 277.5 Hz (about 3 Hz FWHM) that 250 uL does not have. A single uniform field cannot
  make both: the fit chose a large transverse field (108 nT) for the wings and smears the triplet into one broad
  bump; the J sidebands at 134.7 / 137.9 Hz are not reproduced either (model lines at 133.2 / 139.4 Hz at the rate
  bound).
- So the LF_1 per-line decay times of the field-sensitive lines and its field are not meaningful as fitted: the
  sample sees (at least) two field regions, one like the 250 uL sample (~59 nT) and one about twice as strong with
  a wide spread. Plausibly a larger sample volume reaching further into the fringe-field gradient (the name gives
  no volume; to confirm with the collaborator).

## Diagnosis and refit with crop 0.3 s (owner: "the wings come from the baseline correction; higher SNR makes
its effect larger")

- The real part of the processed spectra (no display baseline; runs/figures/acn_lf1_vs_250ul_experiment.png) has no
  positive wings: the sharp lines sit on broad negative troughs, in LF_1 and (smaller against the noise) in 250 uL.
  The wings of the first figure were made by the "model" display baseline from the wrong fit.
- Reprocessing LF_1 (scratch test, zulf_processing.series_spectrum): the broad structure grows with the drift-filter
  window (0.05 / 0.1 s equal, 0.4 / 2 s much larger) and is gone when the record starts at 0.3 s instead of 0.1 s.
  It is an early transient after the field switch that the drift filter (SG, 0.1 s, mirror edges) does not
  remove; the model has no such term, and with the high SNR of LF_1 (J line 15 x the 250 uL height) the fit
  bent the field (108 nT transverse) to imitate it.
- Refit of both with `--crop 0.3` (make_series_entry; series runs/series/acn_lf1_c03, acn250_combined_c03; fits
  runs/processed/<series>_fit, same fit command; script runs/processed/acn_c03_fits.sh). Figure:
  runs/figures/acn_lf1_vs_250ul_crop03.png (fitted spectra and the fits' own models, real part, no display
  baseline; residual offset): both fit to the noise.

| | LF_1, crop 0.3 | 250 uL both, crop 0.3 | 250 uL both, crop 0.1 |
|---|---|---|---|
| 1J(C1,H1) (Hz) | 136.283 | 136.280 | 136.286 |
| B transverse / z (nT) | 30.0 / 56.1 | 41.7 / 44.9 | 39.4 / 43.9 |
| \|B\| (nT) | 63.6 | 61.3 | 59.0 |
| score | 0.0040 | 0.109 | 0.165 |

Decay time T (s) per model line (FWHM = 1 / (pi T)):

| line (Hz) | LF_1, crop 0.3 | 250 uL, crop 0.3 | 250 uL, crop 0.1 |
|---|---|---|---|
| 136.3 (J line, field-insensitive) | 2.57 | 2.85 | 2.26 |
| 134.6-134.7 (J sideband) | 0.48 | 1.29 | 1.05 |
| 137.9-138.0 (J sideband) | 0.31 | 0.52 | 0.32 |
| 270.4-270.5 (weak) | 0.33 | 0.59 | 0.13 |
| 271.4 / 271.6 | 0.43 | 0.54 | 0.44 |
| 272.6 | 0.43 | 0.45 | 0.48 |
| 273.6 / 273.7 | 0.39 | 0.38 | 0.37 |
| 274.7-274.8 (weak) | 0.16 | 0.37 | 0.07 |

- The crop changes the 250 uL decay times by up to 25 % on the strong lines (J line 2.26 -> 2.85 s) and much more
  on the weak outer lines: the crop-0.1 values are biased by the transient too.
- With crop 0.3: the J line decays alike (2.6 vs 2.9 s); the 2J triplet decays alike (0.39-0.43 vs 0.38-0.54 s);
  the J sidebands and the weak outer 2J lines decay faster in LF_1 (0.3-0.5 s vs 0.5-1.3 s).
- Field: |B| 64 vs 61 nT, but split differently (LF_1 30 transverse / 56 z, 250 uL 42 / 45 nT); a different sample
  position or a trade-off between the components (the fit gives no uncertainty for the field here).

Volume (owner: LF_1 is 2 mL). Expected: a larger sample spans a wider range of the fringe-field gradient, so the
field-sensitive lines broaden more (inhomogeneous, about proportional to gradient x size and to each line's
dnu/dB), while the field-insensitive J line keeps the intrinsic decay. Observed: J line alike (2.6 / 2.9 s), J
sidebands and weak outer 2J lines faster in the 2 mL sample (0.2-0.5 s vs 0.4-1.3 s): consistent. The central 2J
triplet decays alike (about 0.4 s in both): not explained by the volume alone (smaller dnu/dB of these lines, or
widths from unresolved lines). Signal per scan 15 x the 250 uL one (J line 0.5 vs 0.034 in the averages) for 8 x
the volume: the rest from geometry (filling closer to the sensor) or other conditions, not checked.

Check without the exponential window (owner: "remove the added decay and fit again"; then a longer record): series
runs/series/acn_2ml_c03_a0_r16 and acn250_c03_a0_r16 (`--crop 0.3 --record 16 --apodization 0`), same fit command
with 2 workers (runs/processed/acn_c03_a0_fits.sh). The window never entered the fitted rates (the model is
rendered through the same processing; data and model |J line| widths agree, 0.285 / 0.283 Hz), and the refit
confirms it:

| line (Hz) | 2 mL, 8 s + 0.3 1/s window | 2 mL, 16 s, no window | 250 uL, 8 s + window | 250 uL, 16 s, no window |
|---|---|---|---|---|
| 136.3 (J line) | 2.57 | 2.56 | 2.85 | 2.69 |
| 134.6 (J sideband) | 0.48 | 0.46 | 1.29 | 1.17 |
| 138.0 (J sideband) | 0.31 | 0.30 | 0.52 | 0.52 |
| 270.4 (weak) | 0.33 | 0.35 | 0.59 | 0.58 |
| 271.4 / 271.6 | 0.43 | 0.42 | 0.54 | 0.53 |
| 272.6 | 0.43 | 0.43 | 0.45 | 0.44 |
| 273.6 / 273.7 | 0.39 | 0.39 | 0.38 | 0.38 |
| 274.7 (weak) | 0.16 | 0.14 | 0.37 | 0.35 |

1J 136.283 / 136.282 Hz, field 30.8 / 55.6 nT (2 mL) and 42.4 / 45.1 nT (250 uL), unchanged. Scores (0.020 / 0.400)
are not comparable with the 8 s fits (twice the points, the added half mostly noise, no window).

## Conclusion

Conditional numerical result. With the record starting at 0.3 s (no early transient), LF_1 and 250 uL give the same
1J (136.28 Hz), a similar field size (64 / 61 nT) and the same decay of the J line (2.6 / 2.9 s) and of the 2J
triplet (0.4 / 0.4-0.5 s); the J sidebands decay faster in LF_1 (0.3-0.5 s vs 0.5-1.3 s). The crop-0.1 fit of LF_1
(first section) is not valid: the fit imitated the transient with a 117 nT field.

## Open points

- The crop 0.1 s of every earlier whole-grid fit (ethanol, the overnight batch) can carry the same transient bias;
  compare with crop 0.3 s for the high-SNR data sets. A model term for the transient would keep the early signal.
- Field uncertainty (transverse / z split) not estimated.
- A second field region (two copies with their own field) was implemented and tested before the diagnosis; it is
  kept unmerged on branch xuehan/second-field (not needed here).
- Position of the 2 mL sample relative to the 250 uL one; dnu/dB of each line to test the broadening quantitatively.

Commits: d3cf6a5 (first section); this commit (crop 0.3 s).
