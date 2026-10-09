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

(running)

## Conclusion

(pending)

## Open points

- (pending)

Commit: (pending)
