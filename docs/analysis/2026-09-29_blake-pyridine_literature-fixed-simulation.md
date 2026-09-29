# Blake pyridine series: literature couplings held fixed (2026-09-29)

- Data: Blake `6_3_26_pyridine_*` processed real spectra (140-200 Hz).
- Question: how well do the literature couplings alone describe the spectra?
- Code: f3d1ba9 (scratchpad script lit_sim.py; results in
  runs/processed/literature_fixed_residuals.json).

## Setup

All 19 couplings fixed at the literature (1J 177.63 / 163.04 / 162.41 Hz, small
couplings as in the literature-start fits, both 4J signs); fitted per spectrum:
isotopologue decay rates, phase delay, gains and shared phase (ratios variant,
signal weighting 2.5 sigma / 4 Hz taper), one start.

## Results

| x | 0.02 | 0.33 | 0.50 | 0.66 | 0.75 | 1.00 |
|---|---|---|---|---|---|---|
| literature, 4J > 0 | 0.620 | 0.645 | 0.577 | 0.531 | 0.512 | 0.468 |
| literature, 4J < 0 | 0.628 | 0.669 | 0.674 | 0.613 | 0.559 | 0.414 |
| couplings fitted | 0.228 | 0.229 | 0.225 | 0.254 | 0.231 | 0.171 |

Figure: runs/processed/literature_fixed_spectra.png (4J > 0).

## Conclusion

At x = 1.00 (neat, nearest to the literature conditions) the main line
positions agree roughly; towards water the literature pattern does not match
(the fit broadens the lines instead). Solvent effects on the couplings are
real, and even at x = 1.00 the literature set leaves 2.4x the fitted residual.

Commit: affebe9
