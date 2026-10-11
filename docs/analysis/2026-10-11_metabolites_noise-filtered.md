# Metabolites: noise-filtered refits of lactic acid, alanine, glycerol, serine (2026-10-11)

- Question (Xuehan): the noisiest panels of the metabolites overlay (lactic acid, alanine, glycerol, serine) should
  be refitted with less noise.
- Data: the z5 averages of the batch (docs/analysis/2026-10-10_overnight_plan.md); new whole-grid spectra
  runs/series/batch_<id>_apod1p5/ (exponential window 1.5 1/s instead of 0.3 1/s, record 4 s instead of 8 s,
  crop 0.1 s). The model goes through the same window, so the fitted rates and couplings are window-free;
  objectives are not comparable with the 0.3 1/s passes (other data vector, less noise).
- Commands: `python scripts/run_batch.py configs/batch_2026-10-10_v9.json --parallel 2 --workers 2 --tag _v9
  --table docs/analysis/2026-10-10_batch_results_v9.md` (1J starts from the previous best pass of each entry);
  crop test `configs/batch_2026-10-10_v10.json` (tag _v10, crop 0.3 s, glycerol and serine).
- Figures: runs/figures/batch/denoise_compare.png (before / after; runs/figures/batch/denoise_compare.py),
  runs/figures/overlay_metabolites.png (now v9; the 0.3 1/s version kept as overlay_metabolites_window0p3.png),
  runs/figures/batch/glycerol_rate_floor.png.

## Results

| entry | pass | objective | 1J C1 / C2 (Hz) | 2J(C1,H2) | 3J(H1,H2) | 2J(C2,H1) |
|---|---|---|---|---|---|---|
| lactic acid | v1 (0.3 1/s) | 0.223 | 147.07 / 129.13 | -4.26 | 5.98 | -3.87 |
| lactic acid | v9 (1.5 1/s) | 0.077 | 147.36 / 129.24 | -4.04 | 6.24 | -3.82 |
| alanine | v1 | 0.383 | 145.46 / 129.89 | -4.19 | 6.68 | -4.54 |
| alanine | v9 | 0.282 | 145.51 / 129.94 | -4.11 | 6.64 | -4.43 |
| glycerol | v2 | 0.759 | 145.41 / 142.68 | -12.0 | 6.13 | -7.84 |
| glycerol | v9 | 0.832 | 142.73 / 143.06 | -5.37 | 5.56 | -2.36 |
| serine | v3 | 0.663 | 146.42 / 144.95 | -3.67 | 5.37 | -3.46 |
| serine | v9 | 0.616 | 146.59 / 144.48 | -3.59 | 3.91 | -1.01 |

- Lactic acid and alanine: the lines stand clearly above the noise and the fit follows them; couplings agree
  with the earlier pass within 0.3 Hz. The fine structure (127.33 / 149.58 Hz lactic acid lines) is smoothed
  by the window (each line about 0.5 Hz wider).
- Glycerol and serine: the small couplings move by 1-6 Hz between passes (glycerol 2J(C1,H2) -12.0 -> -5.4 Hz,
  serine 2J(C2,H1) -3.5 -> -1.0 Hz): they are not determined by these data; the 1J values are stable within
  0.5 Hz (glycerol C1 within 2.7 Hz).
- Crop 0.3 s (v10) is worse for glycerol and serine: their lines decay at 5-10 1/s, and the later start loses
  most of the signal (spectra become ringing). Rejected.
- Glycerol v9 model shows a ripple (period about 10 Hz, 160-245 Hz, broken at the mains gaps 180 / 240 Hz) that
  the data do not show. Not the 0-14 Hz lines outside the grid: a stage-2 refit with rates >= 2 1/s
  (runs/processed/batch_glycerol_v9_rate2_s2; their rate 0.51 -> 3.17 1/s) leaves it unchanged (objective
  0.834). Cause open.

## Conclusion

The noise-filtered fits replace the earlier ones in the metabolites overlay (scripts/overlay_panels.py CHOSEN
"_v9"). 1J values are stable; the long-range couplings of glycerol and serine are not determined. Conditional
numerical results of the given structures.

## Open points

- Origin of the 10 Hz ripple in the glycerol model (crop wings of the strong broad J-band lines with crop 0.1 s
  are the candidate; compare with the model rendered without the SG drift filter).
- Glycerol and serine long-range couplings: a prior or a scan, as for methyl 2-bromopropionate.

Commit: (this commit)
