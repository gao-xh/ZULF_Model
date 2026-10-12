# Pyridine dilution series and pyridinium: refits with the toolbox (2026-10-12)

- Data: the batch entries pyridine-02/10/33/50/66/75/100-mol and pyridinium-6m-acidic (z5 averages; scans in
  original/ and raw/2025-12-09-pyridine-dilutions-and-pyridinium/ on the drive).
- Question (Xuehan): run the pyridine series with the new methods.
- Command: `python scripts/run_batch.py configs/batch_2026-10-10_v28.json --parallel 2 --workers 2 --tag _v28
  --table docs/analysis/2026-10-10_batch_results_v28.md`.

## Results so far


configs/batch_2026-10-10_v28.json: free isotopologue amplitudes, background off, 1J from the bands, line groups
<= 0.4 Hz; window 1.5 1/s for 2 mol %; slow N-H exchange for pyridinium. Finished so far (stage 2; previous v1):
2 mol % 0.584 (window 1.5, not comparable; smoother, owner keeps the previous version), 10 mol % 0.070 (0.064,
similar picture), 33 mol % 0.132 (0.099; the 160 Hz doublet better resolved, 182 Hz still high). Figure
runs/figures/batch/pyridine_v28_compare.png. The overlay keeps v1 unless a v28 panel is clearly better.
Pyridine overlay now coloured by 13C site (C2 / C3 / C4), with a tick at each line of each site
(scripts/overlay_panels.py SITE_STYLES).

Commit: (this commit)

## Stopped (16:40)

50 mol % stage 2 0.045 (v1 0.016), 66 mol % 0.027 (v1 0.016) when stopped: the toolbox makes the pyridine fits
worse (v1 had no slow exchange, field or merged-line problem to fix; background off and free amplitudes cost
here). The rest of v28 was stopped; the overlay keeps v1 for every pyridine. Only pyridinium (N-H+, slow exchange
physically motivated) is run alone with the v28 entry.
