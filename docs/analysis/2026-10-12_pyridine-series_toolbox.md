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

## Coupling trend over the dilution series (v1 stage-2 fits)

`python scripts/coupling_trend.py runs/figures/pyridine_J_trend.png --x 2,10,33,50,66,75,100 --runs
runs/processed/batch_pyridine-{02,10,33,50,66,75,100}-mol_s2 --xlabel "pyridine (mol %)" --csv
runs/figures/pyridine_J_trend.csv` (new script).

- 1J(C2,H2) falls almost monotonically from 180.1 (2 mol %) to 176.9 Hz (100 mol %), 66 mol % below the line
  (177.5); 1J(C4,H4) falls from 164.3 to 162.3 Hz, monotonic; 1J(C3,H3) scatters between 161.5 and 163.5 Hz.
- Long-range J(C,H) and J(H,H) jump by several hertz between neighbouring concentrations (e.g. 3J(H2,H3) 4.6 at
  50 mol %, 9.9 at 66 mol %, 5.5 at 75 mol %) with no trend: not determined by these fits. The linearised error
  bars are far smaller than this scatter (known 10-30x underestimate). A joint fit of the series with smooth
  (monotone) couplings would tie them down.
