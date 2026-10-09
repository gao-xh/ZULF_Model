# ethanol (70 % rubbing alcohol): visually best fit, run and shown in Studio (2026-10-09)

- Data: `runs/series/ethanol_full/series.json` (whole grid 20-380 Hz, 8101 points; the averaged FID of the ethanol
  run in ~/research/zulf/data/raw, see docs/analysis/2026-10-08_ethanol_field-fit.md).
- Question (Xuehan): continue from the best ethanol fits towards the fit that looks best against the data (the
  slide overlay: 120.7 Hz and 204 Hz peaks too low in the simulation, broad bumps at 196-200, 225-228 and
  245-247 Hz missing, 251 Hz low), and show it in Studio.
- Code: 3c36182; run logs: runs/studio/fits/20261009-134414_ethanol/RUN_LOG.md

## Setup

Started from Studio (JSON API, so the window follows the fit) at the applied
ethanol_spinsystem2 result (objective 0.1210: 1J 125.447 / 141.525 Hz, 2J -2.400 / -4.796 Hz, 3J 7.132 Hz,
B transverse 2.6 nT, B z 64.6 nT). Changes against that fit, aimed at the visible misfits:

- decay-rate bounds 0.1-40 1/s (were 0.2-15; four 13C@C2 rates sat at 15 1/s, so broad features could not form);
- free amplitude per isotopologue (spin_system fixed_weights false): the measured 13C@C1 / 13C@C2 intensities
  need not follow the abundance ratio (relaxation and polarisation transfer differ by site);
- one rate family per data peak (`--family-edges peaks`), field fitted, guard rows default, 4 starts on 4 workers,
  300 evaluations per start, no component search (the start is the previous optimum).

## Results

### Stage 1: wider rate bounds, free isotopologue amplitudes (runs/studio/fits/20261009-134414_ethanol, 1383 s)

| start | objective |
|---|---|
| 0 | **0.11822** |
| 1 | 0.12269 |
| 2 | 0.17834 |
| 3 | 0.12403 |

Final score 0.1183 (was 0.1210), residual 0.348 (was 0.352). Couplings (J_std from the linearised covariance):
1J(C1,H1) 125.448 +- 0.002, 1J(C2,H2) 141.512 +- 0.006, 2J(C1,H2) -2.403 +- 0.008, 2J(C2,H1) -4.775 +- 0.013,
3J(H,H) 7.153 +- 0.007 Hz: all within 0.03 Hz of the previous fit, so the couplings were already settled.
Field: B transverse at its lower bound 0, B z 66.6 nT; delay -3.6 ms. 8 of 76 decay rates at a bound (0.1 or
40 1/s).

Seen in the overlay (runs/figures/ethanol_overlay_visual.png, scripts/overlay_figure.py with the run's options):
the 125.5 Hz line now reaches the data (was about 8 % low), 204, 217 and 256 Hz closer; still missing: the
120.7 Hz line (simulation about 90 against 215), the 124 Hz shoulder and the broad bumps at 196-200, 225-228 and
245-247 Hz.

### Stage 2: residual-peak rows (runs/studio/fits/20261009-140825_ethanol)

From the stage-1 result, `--residual-peaks 5` (data peaks the model leaves unexplained weighted in a refit of
the best candidates), 2 starts, 250 evaluations each, otherwise as stage 1.

Stopped after 38 min (Xuehan): the residual-peak rows act only after the starts, and the starts only repeated
stage 1 from its optimum (best 0.11836 and 0.11942 after about 120 evaluations each, against 0.11822); about
1.5 h more before the residual-peak stage. (Studio's Stop ended only the main process; its two workers were
stopped by PID, and Stop now ends the workers too, d0a76ac.)

### Why 120.7 Hz is missing

The model has a 13C@C1 line at 120.84-120.87 Hz (relative 0.33 of the strongest line); in stage 1 it shared a
4.6 1/s family (FWHM 1.5 Hz), so it came out broad and low (about 90 against 215 in the overlay). `peaks` finds
families in the data's smoothed magnitude: there the line is about 6 sigma and was not taken as a peak (the
peaks found near it were 115.38, 118.71 and 123.69 Hz). Not mains: the 120 Hz harmonic is 0.7-1.3 Hz away (and
excluded +-0.4 Hz), and acetonitrile and isopropanol show nothing at 112-130 Hz; it is weak, though (6.7 sigma in
the even half, not separated in the odd half). New option `--family-edges lines` (one family per model line,
D59 amendment, 03f77d2): 38 families, 120.87 Hz on its own (123.96 / 124.12 Hz, 0.16 Hz apart, share one).

### Stage 3: one decay rate per model line (runs/studio/fits/20261009-144726_ethanol)

From the stage-1 result: `--family-edges lines`, `--residual-peaks 5`, rate bounds 0.1-40 1/s, free
isotopologue amplitudes, field fitted (transverse start 20 nT, as Studio starts a zero component), 4 starts,
300 evaluations each, no component search.

| start | objective |
|---|---|
| 0 | **0.11067** |
| 1 | 0.13214 |
| 2 | 0.51996 |
| 3 | 0.27314 |

Residual-peak stage (1-2 assignable residual peaks of 16-18 per candidate, so it changed little): final score
**0.1117** (the candidate ranked best on the common residual-peak windows; start 0 alone 0.1107), residual 0.338
(stage 1 0.348, ethanol_spinsystem2 0.352); 3247 s. Couplings (J_std from the linearised covariance):
1J(C1,H1) 125.410 +- 0.002, 1J(C2,H2) 141.389 +- 0.008, 2J(C1,H2) -2.424 +- 0.012, 2J(C2,H1) -4.821 +- 0.029,
3J(H,H) 7.043 +- 0.011 Hz. Field B z 63.6 nT, B transverse 0 (at its bound); delay -3.7 ms; 2 of 76 rates at a
bound. Against stage 1 the couplings moved by 0.04-0.12 Hz, more than their linearised errors: the decay-rate
model (families) changes them at that level, so the model spread, not J_std, is their uncertainty.

Visually (runs/figures/ethanol_overlay_lines.png, start 0 point): 211 Hz now reaches the data (about 620 against
650), the 251 Hz multiplet shows its splitting, 217, 256 and 128 Hz shapes closer. The 118-125 Hz region still
looked low.

### The 118-125 Hz "misfit" was the display baseline (D62)

Overlays with four display baselines (zulf_processing.display_baseline method none / separate / shared /
residual) at 114-131 Hz: without a baseline data and model agree there (the crop's rolling baseline is in both);
"separate" (the earlier figures) estimated one baseline on the data and another on the model and left that
roll in the data only, which looked like a missing line. One baseline for both curves keeps the residual of the
fit; estimated on the noise-free model ("model", now the default of paper_figure.py, overlay_figure.py and
Studio) it adds no noise ripple to the model, as "shared" (estimated on the data) did. Final overlay:
runs/figures/ethanol_overlay_final.{png,svg,pdf} and ethanol_overlay_final_120.png (final fit.json, baseline
"model").

### Halves (even / odd scans) with the stage-3 model

Whole-grid series of the two halves made with the same command as ethanol_full
(`make_series_entry.py --fid .../z5/average_{even,odd}.npy --id ethanol --record 7.5 --exclude 81.5,86`; the
spectrum of average_fid.npy made this way equals ethanol_full point for point) and the fit ranges of
ethanol_full copied into them (runs/series/ethanol_full_{even,odd}). Fits: the stage-3 command (same family
edges, rate bounds, free amplitudes), started at the final stage-3 couplings and field, 2 starts each, no
residual-peak stage (runs/processed/ethanol_lines_{even,odd}).

(running)

## Conclusion

Conditional numerical result for the ethyl spin system with a uniform static field: the visually best fit so
far is stage 3 (one decay rate per model line, free isotopologue amplitudes, rates 0.1-40 1/s), score 0.1117
(from 0.1210); with one display baseline for both curves the overlay shows no systematic misfit at the strong
lines. The couplings are stable to about 0.1 Hz across the decay-rate models (1J 125.41-125.45 / 141.39-141.53,
2J -2.40 to -2.42 / -4.78 to -4.88, 3J 7.04-7.15 Hz); the field direction fit (transverse at 0) differs from the
earlier fits (2.6-12 nT transverse) and is not determined by these data at this level.

## Open points

- Halves (even / odd scans) and more starts for stage 3; profile 3J(H,H).
- Small features still above the model: 126.9 Hz, 196-200 Hz, 245-247 Hz (no or weak model lines).
- The 1-2 assignable residual peaks per candidate: which transitions (fit.json residual_peaks).

Commit: (this commit)
