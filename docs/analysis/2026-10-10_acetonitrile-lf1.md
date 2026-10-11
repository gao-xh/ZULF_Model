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

Field sensitivity of every line vs the volume effect (scripts/field_sensitivity.py on the 16 s window-free fits;
figure runs/figures/acn_rate_vs_field_sensitivity.png, numbers in runs/figures/acn_rate_vs_field_sensitivity.json).
dnu/d|B| from the fitted spin system at each fit's field (|B| scaled by 1.001, direction fixed). Expected:
rate difference = pi |dnu/dB| dB (dB: the extra field spread of the 2 mL sample).

| line (Hz) | dnu/d\|B\| (Hz/nT) | rate 2 mL (1/s) | rate 250 uL (1/s) | difference (1/s) | dB (nT) |
|---|---|---|---|---|---|
| 134.6 (1J, Delta m = -1) | -0.026 | 2.18 | 0.86 | +1.32 | 15.9 |
| 136.3 (1J, Delta m = 0) | +0.000 | 0.39 | 0.37 | +0.02 | - |
| 138.0 (1J, Delta m = +1) | +0.027 | 3.39 | 1.93 | +1.46 | 17.2 |
| 270.4 (2J) | -0.034 | 2.90 | 1.74 | +1.16 | 10.7 |
| 271.4 / 271.6 (2J, one shared rate) | -0.019 / -0.016 | 2.38 | 1.89 | +0.49 | 8.4 / 9.9 |
| 272.6 (2J, center) | +0.000 | 2.31 | 2.30 | +0.01 | - |
| 273.6 / 273.8 (2J, one shared rate) | +0.016 / +0.019 | 2.58 | 2.67 | -0.09 | -1.8 / -1.5 |
| 274.8 (2J, weak) | +0.035 | 7.36 | 2.84 | +4.52 | 41.3 |

The slope-zero lines (1J 136.3 and 2J 272.6 Hz) do not change; the 1J sidebands give a consistent dB of 16-17 nT;
the 2J lines follow the trend but scatter (8-11 nT for the low-frequency side, about 0 for 273.6/273.8, 41 nT for
the weak 274.8 Hz line). Proportionality is supported qualitatively; a quantitative test needs uncertainties of the
rates (not estimated yet).

Lower-noise refit of both (record 8 s, window 1.0 1/s, crop 0.3 s; runs/series/acn_2ml_c03_w1, acn250_c03_w1; fits
runs/processed/<series>_fit): 2 mL 1J 136.281 Hz, B 28.6 / 56.6 nT; 250 uL 1J 136.271 Hz, B 39.9 / 45.0 nT. Lines of the
two fits are now matched as the same transition by dnu/dB (scripts/field_sensitivity.py; the field differs, so
the outer lines move by about 0.1 Hz). Rate difference 2 mL - 250 uL (1/s): 136.3: +0.09, 272.6: +0.18,
271.4/271.6: +0.60, 273.6/273.8: +0.02, 134.6: +1.34, 138.0: +1.10, 270.4: +1.81; 274.8 Hz: the 2 mL rate went to the
bound (40 1/s, the weak line suppressed), not usable. Fit through the origin without 274.8 Hz: dB difference
13.7 nT, rms 0.35 1/s (runs/figures/acn_delta_rate_vs_slope_w1both.png). The rates of the same sample change by
0.1-0.3 1/s between processings (e.g. 250 uL 134.6 Hz: 0.86 with 16 s / no window, 0.54 with 8 s / 1.0 1/s): an
estimate of their uncertainty. Next: tie the rates of mirror lines and fit R = R0 + pi |dnu/dB| dB directly.

Mirror-pair asymmetry checks (8 s, window 1.0 1/s fits):
- One rate free, everything else held (scratch scripts on the fits' own objective): 250 uL 273.5/273.7 started at
  1.7 or 2.3 returns to 2.48, 271.4/271.6 started at 2.3 returns to 1.70; 2 mL cross-starts return to 2.31 / 2.50.
  Profiles (runs/figures/acn_w1_rate_profiles_both.png) have clear minima (about +-0.3 1/s for +1 % in 250 uL,
  +-0.1 in 2 mL). Within this model the high-frequency member is faster in both samples.
- Magnitude spectrum (2J band, the two pair rates free): 2 mL 2.19 / 2.40, 250 uL 1.43 / 2.26 1/s; one shared rate
  costs +64 % / +117 %: not a phase effect (runs/figures/acn_2j_magnitude_fit.png).
- No fixed narrow line at 273-275 Hz in any other data set (4000 / 8333 Hz); acetonitrile is the only 2000 Hz data,
  so an aliased instrument line is not excluded. Residuals (runs/figures/acn_residual_mirror_pairs.png): no isolated
  extra peak at 273.6; 2 mL shows band-wide ripples in the 2J residual (line-shape mismatch).
- Second-order Zeeman shifts make the high-frequency member slightly more field sensitive (1-2 %), the right sign
  but 0.01-0.03 1/s instead of 0.2-0.8 1/s.
- Voigt (new fit_joint_series --fit-sigma: one Gaussian width per isotopologue on top of the per-line Lorentzian
  rates; runs/processed/<series>_voigt_fit, same family edges): 2 mL sigma 0.008 Hz (none; score unchanged), 250 uL
  sigma 0.069 Hz (score -1.3 %, all rates about 0.2 1/s lower); the asymmetry is unchanged. A uniform Gaussian
  cannot describe field broadening, which scales with each line's dnu/dB; next: sigma_k = |dnu/dB|_k sigma_B.
- Model-free envelope decay (demodulate, Gaussian low-pass): the J line gives 0.39 / 0.38 1/s as the fits, but the
  sidebands reach the noise floor within the window; inconclusive.

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
