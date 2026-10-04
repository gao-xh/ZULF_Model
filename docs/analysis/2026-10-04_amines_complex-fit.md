# Ethylenediamine, triethylamine, N-ethylmethylamine: complex known-structure fits (2026-10-04)

- Data: confirmed samples fde3fbb2 (ethylenediamine, H2N-CH2-CH2-NH2), 4322bdfc
  (triethylamine, N(CH2CH3)3) and e66a4b08 (N-ethylmethylamine,
  CH3-CH2-NH-CH3); NMRduino average FIDs, 4000 Hz, 65516 points (uploads, read
  from ZULF_DATA_DIR, not committed). Identities given by the user after the
  blind analyses (docs/ANALYSIS_LOG.md); these are known-structure fits.
- Question (user): fit the other amines with the current fitter, as for
  isopropylamine (docs/analysis/2026-10-03_isopropylamine_complex-fit.md):
  fit quality, couplings, NH / 15N traces, and whether localized misfits like
  the isopropylamine 121.69 / 122.28 Hz pair also appear.
- Code: aa07a76 at the start; run logs runs/processed/{eda,tea,nema}_fast*/RUN_LOG.md;
  scratchpad amines/make_series.py, auto_edges.py, run_stage1.sh, run_stage2.sh.

## Processing (same as isopropylamine)

Crop 0.1-8.1 s (start sample 400), exponential 0.3 1/s, zero fill 3 (0.0417 Hz
grid), phase from the instrument calibration (configs/confirmed_samples.json:
phase0 176.3 deg at the switching edge - 0.033 ms); edges 3.470 / 3.439 / 3.477
ms. The complex fit refines a residual zero-order phase in the gain. Fit ranges
(mains harmonics n x 60.06 Hz +-0.4 Hz out; the regions where 15N lines would
lie kept even without signal; the broad low-frequency features below 55 Hz,
which no model line explains, left out):

| sample | ranges (Hz) |
|---|---|
| ethylenediamine | 85-119.6, 120.4-179.6, 180.5-237.5 |
| triethylamine | 85-119.6, 120.4-142.1, 180.5-215.1, 231.6-239.6, 240.5-263.2 |
| N-ethylmethylamine | 55-59.6, 60.5-77.3, 102.2-119.6, 120.4-153.5, 177.9-179.6, 180.5-215.9, 231.6-239.6, 240.5-276 |

Model lines (line table, motif defaults): ethylenediamine 13C@C1 (x2) at
190-206 Hz, slow exchange 15N lines at 98-110 Hz; triethylamine 13C@C1 (CH2)
185-196 Hz, 13C@C2 (CH3) 123-127 and 236-258 Hz; N-ethylmethylamine 13C@C1
(CH3 of ethyl) 123-127 / 237-258, 13C@C2 (CH2) 192-206, 13C@C3 (N-CH3)
128-134 / 259-266 Hz, slow exchange 15N lines at 62-67 Hz.

## Stage 1: fast exchange, one rate per isotopologue

fit_joint_series --real-only false --shape free --exchange fast, signal
threshold 2.5 / taper 4, peak penalty 5 (smooth 0.03, min sigma 3), 32 starts
(spread 2 Hz) from the earlier regression couplings (runs/regression/final;
2J(C2,H1) of N-ethylmethylamine started at -4.5 instead of -11.9 Hz):

| sample | objective | relative residual | time |
|---|---|---|---|
| ethylenediamine (eda_fast) | 0.1015 | 0.302 | 266 s |
| triethylamine (tea_fast) | 0.3115 | 0.544 | 710 s |
| N-ethylmethylamine (nema_fast) | 0.3129 | 0.548 | 700 s |

(isopropylamine at this stage: 0.146.) Couplings (Hz):

- ethylenediamine: 1J(C,H) 132.40, 2J(C,H) -2.93, 3J(H,H) 5.63.
- triethylamine: 1J(CH2) 130.71, 1J(CH3) 124.97, 2J(C1,H2) -4.61, 2J(C2,H1)
  -3.18, 3J(H1,H2) 7.22, 3J(C1,HX) (through N) -2.21.
- N-ethylmethylamine: 1J(CH3) 124.91, 1J(CH2) 138.33, 1J(N-CH3) 131.76,
  2J(C1,H2) -3.49, 2J(C2,H1) -5.91, 3J(H1,H2) 7.25, 3J(C2,H3) 11.45 (through
  N, implausibly large), 3J(C3,H2) 3.73.

Misfits (runs/processed/*_fast/spectra.png): one rate per isotopologue cannot
carry both the narrow and the broad lines (triethylamine: 125 Hz lines too
weak, 250 Hz too strong; N-ethylmethylamine: 195-205 Hz CH2 lines too weak).

## Stage 2: rate families

Running (run_stage2.sh): edges from the stage-1 lines (auto_edges.py: lines of
relative amplitude >= 0.05 clustered at gaps > 0.5 Hz, at most 18 families),
remote couplings at 0 in the motif freed (triethylamine J(H1,HX), J(H2,HX),
J(C2,HX); N-ethylmethylamine J(C1,HC3), J(HC1,HC3), J(HC2,HC3), J(C3,HC1);
seeds 0 and +-0.8 Hz), rate bounds 0.2-15 1/s, 4 starts.

Ethylenediamine stage 2 (eda_fast_fam, 4 families 195.3 / 199.5 / 202.4 Hz,
123 s): objective 0.0719, relative residual 0.253; 1J 132.06, 2J -2.86, 3J(H,H)
6.35 Hz. The data show sharp pairs (190.6 / 192.1, 197.8 / 198.8, 200.2 Hz,
widths 0.3-0.5 Hz); the model has five lines (190.6, 198.09, 198.55, 199.5,
204.6 Hz) and covers the pairs with broad lines (rates 7.9 / 5.2 / 11.0 / 2.5
1/s, FWHM 2.5 / 1.6 / 3.5 / 0.8 Hz): broadening is too cheap where positions
are missing. The fitted delay is off: +2.3 ms (stage 1) and the -10 ms bound
(stage 2), against -3.6 ms for isopropylamine, triethylamine and
N-ethylmethylamine (near 200 Hz a delay is ambiguous by about 1 / f inside one
band). New option `--phase-delay-bounds lo,hi` (ms) in fit_joint_series
(instrument prior -4.6,-2.6); rerun eda_fast_d / eda_fast_fam_d running.
Structural candidates next: two vicinal couplings J, J' of the AA'BB' CH2-CH2
unit (line table: J 3 / J' 9.7 Hz moves the 199.5 Hz line to 200.14 Hz, the
data line is at 200.2 Hz; no line at 192.1 Hz yet) and slow NH2 exchange.
