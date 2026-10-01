# Blake pyridine series: line widths, peak weighting and a missing-peak penalty (2026-10-01)

- Data: Blake `6_3_26_pyridine_{02,33,50,66,75,100}` processed real spectra
  (140-200 Hz), x = pyridine mole fraction 0.02 ... 1.00 (not committed).
- Question: the best data-only solution of the 48-more run
  (2026-10-01_blake-pyridine_ms48-continuation.md) loses the small peak about
  1.5 Hz above the strongest line (169.9 / 169.5 / 169.1 Hz at x 0.33 / 0.50 /
  0.66) that the earlier 48-start best reproduced. Is it the single line width
  per isotopologue, the weighting, or a trade with other peaks; and can the
  objective keep such peaks without giving up the tall ones?
- Code: scripts/fit_joint_series.py (`--peak-penalty`, `--peak-smooth`),
  zulf_core/solver/forward.py (`height_factor`, `signal_height_power`);
  scratchpad scripts ratefam.py, rescore_height.py, region_trade.py,
  peakpen_check.py, smooth_conv.py, smooth_eval.py.

## A. Own decay rate for the lines near the small peak

Couplings held at a stored solution; per spectrum the decay rates and delay
refitted with (a) one rate per isotopologue, (b) an own rate for lines in a
window main + 0.8 ... main + 2.3 Hz (rates outside tied), (c) an own rate for
the whole main cluster (main +- 2.5 Hz).

| couplings | x | peak data | (a) model | (b) model | (b) score | (b) 13C2 rate in window / outside (1/s) |
|---|---|---|---|---|---|---|
| 48-start best | 0.33 | 0.350 | 0.285 | 0.356 | -10.0 % | 0.81 / 3.35 |
| 48-start best | 0.50 | 0.309 | 0.253 | 0.329 | -10.7 % | 0.92 / 3.14 |
| 48 more best | 0.33 | 0.350 | 0.225 | 0.255 | -3.5 % | 0.88 / 2.26 |
| 48 more best | 0.50 | 0.309 | 0.191 | 0.186 | -1.3 % | 1.17 / 2.11 |

With the 48-start couplings the 13C2 lines there are narrower (T2* about
1.1 s against 0.3 s) and the peak is reproduced; with the 48-more couplings the
lines are not in the right place and a width cannot help. At x >= 0.66 the
window rate runs to its bound (20 1/s), i.e. a free width is used to remove
lines: per-line widths are not safe without a constraint.

## B1. Peak-height weighting (dropped)

Option `signal_height_power` p (default 0): weight times (height of the nearest
data core / tallest)^-p, height floored at 0.1. A first draft that used the
local running maximum raised empty gaps to the floor (x20) and was replaced.
Rescoring the stored best solutions (p = 1): the 48-start best moves from
+8.3 % to +17 % behind the 48-more best. Per 2 Hz band the two solutions trade
peaks: the 48-start best is better at 168-170 Hz and 184-186 Hz, the 48-more
best at 152-162 Hz. Raising all small peaks raises both sides of the trade, and
every reweighting lowers the share of the tall peaks. Not pursued.

## B2. Missing-peak penalty

The ordinary residual and its weights are unchanged. For every data peak top
(prominence >= 12 % of the maximum and >= 2 noise sigma; 9-15 tops per
spectrum, the small peak found at x 0.33, 0.50 and 0.66) one extra row
lambda * min(max of r within +-0.15 Hz, 0): zero while the model reaches the
data top, growing when it stays below. Jacobian rows are combinations of the
ordinary rows (checked against finite differences, relative error 9e-5).

Stored solutions rescored, lambda = 5 (couplings and spectrum parameters as
stored):

| solution | data | penalty | total | largest misses |
|---|---|---|---|---|
| 48-start best | 0.2881 | 0.0388 | 0.3268 | 164.7, 150.3, 160.2 Hz |
| 48 more best | 0.2661 | 0.0638 | 0.3298 | 164.1, 169.5 / 169.9, 151.0 Hz |
| staged | 0.3005 | 0.0748 | 0.3753 | |
| literature start, 4J < 0 | 0.3040 | 0.0762 | 0.3802 | |

The hard rows have kinks; a 16-start run with them was stopped after 50 min
(each start ran to the 200-evaluation cap). Smooth rows for the optimiser
(`--peak-smooth 0.03`): soft max (log-sum-exp) over the window and soft hinge,
width 3 % of each peak's height in residual units; solutions are rescored and
ranked with the hard rows.

One start from the 48-start best, 200 evaluations (both reach the cap):

| objective | data | hard penalty | small peak model / data, x 0.33 / 0.50 / 0.66 |
|---|---|---|---|
| smooth penalty, lambda 5 | 0.2906 | 0.0174 | 0.323 / 0.350, 0.300 / 0.309, 0.172 / 0.195 |
| no penalty | 0.2867 | 0.0368 | 0.291 / 0.350, 0.254 / 0.309, 0.166 / 0.195 |

The penalty costs 1.4 % on the data score and keeps the small peak.

## 48 starts with the smooth penalty

    python scripts/fit_joint_series.py --series series.json \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162.5, "A4": 161.5}}' \
      --couplings "$(cat pyr_lit.json)" --prior-sigma-hh 0.5 --prior-sigma-ch 1.5 --prior-weight 10 \
      --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 \
      --from-joint runs/processed/joint_new_ms48/fit.json --seeds seeds_peakpen.json \
      --starts 40 --spread 1.5 --prior-starts 8 --seed 23 --workers 4 \
      --out runs/processed/joint_peakpen_smooth_ms48

Seeds (10): the two single-start solutions above, the 48-start best, the
48-more best (total and data-only) and the 5-solution near-best set of the
48-start run. Started 2026-10-01 20:06 UTC (one BLAS thread per worker).

48 starts in 2096 s (run log runs/processed/joint_peakpen_smooth_ms48/RUN_LOG.md).
Best: start 11 (a perturbed start around seed 1), hard total 0.2911; next
0.3248 (+11.6 %), 1 solution within 3 %. Rescored with the same hard rows
(scratchpad pp_result.py):

| solution | data | penalty | prior | total |
|---|---|---|---|---|
| penalty run best | 0.2557 | 0.0146 | 0.0209 | 0.2911 |
| 48-start best | 0.2879 | 0.0390 | 0.0205 | 0.3474 |
| 48 more best | 0.2659 | 0.0643 | 0.0261 | 0.3562 |

The new best is lower on the data alone than every earlier solution (best
data-only before: 0.2634) and misses fewer peaks. Relative residuals 0.221 /
0.209 / 0.202 / 0.201 / 0.236 / 0.201 (x 0.02 ... 1.00). Main peak height
(model, data = 1): 1.02 / 0.96 / 0.98 / 0.99 / 1.02 / 0.95. Small peak (model /
data): 0.207 / 0.350, 0.203 / 0.309, 0.142 / 0.195 at x 0.33 / 0.50 / 0.66;
better than the 48-more best (0.158, 0.127, 0.157), below the 48-start best
(0.285, 0.252, 0.166) and below the single smooth start (0.323, 0.300, 0.172).
Largest remaining misses: x 0.66 150.5 Hz, x 0.75 164.1 Hz, x 0.33 169.9 Hz,
x 0.50 169.5 Hz.

1J(C2,H2) 178.79 -> 177.32, 1J(C3,H3) 164.22 -> 162.73, 1J(C4,H4)
162.09 -> 156.93 Hz (x 0.02 -> 1.00), as in the earlier solutions. Small
couplings still change by several hertz across the series, e.g. J(C4,H3)
-3.42 -> 6.42 Hz, J(C2,H3) 6.45 -> -4.72 Hz.

Figures: runs/processed/joint_peakpen_smooth_ms48/spectra.png,
couplings_vs_x.png, J_trends_compare.png (four solutions, literature), compare_zoom.png (x 0.33 / 0.50 / 0.66, three solutions,
peak tops marked).

## Conclusion

The smooth missing-peak penalty finds a solution that is better on the data
alone and on missed peaks than every earlier one, so the penalty does not cost
the tall peaks. The small peak 1.5 Hz above the strongest line is still only
about 60 % reproduced: in the data it is narrower than any model line nearby
(compare_zoom.png), which agrees with test A (own, smaller decay rate for
those 13C2 lines). The best solution again comes from a perturbed start and
stands 11.6 % below the next: the search has not converged.

## Open points

- Narrower lines near the small peak: a constrained per-family rate (bounded
  ratio to the isotopologue rate) in the joint fit, instead of a free width.
- More starts around the new best (search not converged).
- Reliability classes on the penalised objective (reliability_series.py does
  not include the penalty rows yet).

Commit: see git log (this file).
