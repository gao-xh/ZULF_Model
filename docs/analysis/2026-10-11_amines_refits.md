# Amines: refits with background off, 1J from bands, noise filter (2026-10-11)

- Question (Xuehan): improve the amines overlay with the improvements found on the metabolites.
- Command: `python scripts/run_batch.py configs/batch_2026-10-10_v13.json --parallel 3 --workers 2 --tag _v13
  --table docs/analysis/2026-10-10_batch_results_v13.md` (all entries: `background_order -1`,
  `one_bond_from_data`, 1J starts from the previous best pass; window 1.5 1/s and record 4 s for diethylamine,
  n,n-dimethylethylenediamine, sec-butylamine). Applied to every entry as a probe; the overlay switched only
  where better.
- sec-butylamine stage 2 (6.4 h, 24 evaluations) and dipropylamine stage 1 (7 h) were stopped by PID;
  sec-butylamine fit.json is the best monitor point (objective 0.0781, `stopped_early`); dipropylamine has no v13.
- Figures: runs/figures/batch/amines_v13_compare.png (previous vs v13; runs/figures/batch/compare_versions.py),
  runs/figures/overlay_amines.png (version printed in every panel).

## Results

| entry | previous | v13 objective | verdict |
|---|---|---|---|
| tert-butylamine | v1 s2 (0.121) | 0.169 | no gain, keep v1 |
| isopropylamine | v2 s2 (0.063) | 0.133 | no gain, keep v2 |
| triethylamine | v1 s2 (0.090) | 0.090 | same, keep v1 |
| ethylenediamine | v1 s2 (0.062) | 0.078 | small lines 190-193 Hz no better, keep v1 |
| diethylamine | v2 s2 (0.559, almost no model signal) | 0.855 (other data vector) | **v13**: 125, 251, 255 Hz lines fitted; model lines near 200 Hz not seen in the data |
| n,n-dimethylethylenediamine | v4 s2 (0.124) | 0.074 | **v13**: all bands follow the data |
| sec-butylamine | v2 s2 (0.070) | 0.078 (stopped early) | **v13**: bands follow the data with less noise |
| dipropylamine | v2 s2 (0.066) | - | keep v2 |

Objectives of passes with the noise filter are not comparable with the 0.3 1/s passes.

## Open points

- The 2J (XA3) band near 245-255 Hz is too weak in the model for isopropylamine, triethylamine and (less)
  tert-butylamine, in every pass. It belongs to the same isotopologue as the fitted J band, so free amplitudes
  cannot help; candidates: a faster 2J-band decay limited by the rate families, field sensitivity of the 2J
  lines, or the early switching transient (compare crop 0.3 s for these slow-decaying spectra).
- diethylamine: spurious model lines near 200 Hz.

## Slow N-H exchange and the 15N isotopologue (v14, v15; owner's suggestion)

With fast exchange the N-H protons are decoupled and the 15N isotopologue has only small couplings to
carbon-bound protons, so its lines fall below the fitted range; with slow exchange the N-H protons stay in the
spin system and both the 13C and the 15N isotopologues gain lines (1J(N,H)). Triethylamine has no N-H: neither
applies. Command: `python scripts/run_batch.py configs/batch_2026-10-10_v14.json --parallel 2 --workers 2 --tag
_v14` (exchange slow, static; start from the chosen pass of each entry). v15 (one fitted exchange rate k,
`nh_exchange` 50 1/s) did not finish one evaluation of stage 1 in 46 min (exact Liouville model with the
component search) and was stopped.

| entry | previous (same data) | v14 objective | 1J(N,H) (Hz) | verdict |
|---|---|---|---|---|
| tert-butylamine | v1 s2 0.121 | 0.110 | -75.76 +- 0.05 | **v14**, 2J band slightly better |
| isopropylamine | v2 s2 0.063 | 0.116 (stage 2 stopped) | -75.80 (stage 1) | worse, keep v2 |
| ethylenediamine | v1 s2 0.062 | **0.033** | -65.00 +- 0.31 | **v14**: the small lines at 190-193 Hz now fitted |
| diethylamine | v13 s2 0.855 | **0.535** | -75.89 +- 0.15 | **v14**: 251 / 255 Hz lines better, spurious 200 Hz lines weaker; 2J(C2,HN1) -11.6 Hz implausibly large |

Figure runs/figures/batch/amines_v14_compare.png; the overlay draws v14 for these three.
Isopropylamine and triethylamine keep a too-weak 2J band: not explained by N-H coupling.

Conditional numerical results of the given structures.

Commit: (this commit)
