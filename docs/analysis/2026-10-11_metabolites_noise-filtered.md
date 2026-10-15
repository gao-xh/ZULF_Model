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
  0.834). Cause: the per-band background (below).

## Glycerol ripple: cause and fix (later the same night)

- Decomposition of the v9 model (forward predict: component spectra and background): the components are below
  0.0005 in 225-240 Hz, the model 0.0023 there. The ripple is the per-band linear background
  (`default_fit_base` background_order 1), which is smooth in the record's frame and multiplied by the
  first-order phase correction, i.e. a ripple of period 1/(crop + delay) = 10 Hz, broken at the band gaps.
  With signal weighting the empty regions cost little, so the fit used it to absorb broad structure near the
  lines. Not the AsLS display baseline (raw model, no AsLS) and not broad lines (rate caps 8 and 12 1/s:
  unchanged).
- New option `fit_joint_series.py --background-order` (default unchanged: 1) and run_batch entry key
  `background_order`.

| version | change | objective | 1J C1 / C2 (Hz) | result |
|---|---|---|---|---|
| glycerol_v9_s2 | background linear | 0.832 | 142.73 / 143.06 | ripple 160-245 Hz |
| glycerol_v9_bg0_s2 | background constant, stage 2 only | 0.869 | 142.76 / 143.03 | ripple smaller |
| glycerol_v9_bgoff_s2 | background off, stage 2 only | 1.011 | 142.79 / 143.08 | no ripple, CH2 band shifted |
| glycerol_v11_s2 | background off, both stages | 0.989 | 144.05 / 141.08 | no ripple; half of the middle CH2 doublet (213.3 Hz) and the 147 Hz CH line missed |
| **glycerol_v12_s2** | as v11, 1J starts from the bands | 0.962 | 142.32 / 143.13 | both bands fitted, no ripple |
| serine_v11_s2 | background off | 0.857 | 146.55 / 144.52 | as serine_v9 |

Objectives with and without background are on the same data, but the background lowers the objective by
absorbing structure; v12 is chosen by the bands (runs/figures/batch/glycerol_zoom.png, raw real part, model
lines marked; runs/figures/batch/glycerol_background.png). glycerol_v12 couplings: 2J(C1,H2) -1.98,
3J(C1,H3) 2.93, 3J(H1,H2) 5.11, 2J(C2,H1) -1.98 Hz (not determined: they moved by several Hz between passes).
The overlay now draws glycerol_v12 and serine_v11, with the version printed in every panel.

## Conclusion

The noise-filtered fits replace the earlier ones in the metabolites overlay (scripts/overlay_panels.py CHOSEN
"_v9"). 1J values are stable; the long-range couplings of glycerol and serine are not determined. Conditional
numerical results of the given structures.

## Open points

- Other fits with the default linear background may carry the same ripple: check the model away from the lines.
- Glycerol and serine long-range couplings: a prior or a scan, as for methyl 2-bromopropionate.

Commit: (this commit)
