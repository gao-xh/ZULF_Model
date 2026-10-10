# acetonitrile fringe field LF_1 vs 250 uL: decay times per line (2026-10-10)

- Data (drive only): `original/2026-09-23-acetonitrile-fringe-field/ACN_fringe_field_ms1_LF_1/`, 2000 scans,
  2026-09-23 to 09-25, 2000 Hz, no volume in the name. z5 average (--exclude-z 5): 1744 scans,
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

## Conclusion

Conditional numerical result. Intrinsic (field-insensitive) decay: T 2.6 s (LF_1) vs 2.26 s (250 uL). The
field-sensitive lines of LF_1 are not described by a uniform-field model; the 250 uL decay times (triplet 0.37-0.48 s,
J sidebands 0.3-1.1 s) have no valid LF_1 counterpart from this fit.

## Open points

- Fit LF_1 with two field components (two copies of the acetonitrile model with separate fields and amplitudes,
  or a field distribution) to get the decay times of the triplet and of the wings separately.
- Sample volume and position of LF_1 (collaborator).

Commit: (this commit)
