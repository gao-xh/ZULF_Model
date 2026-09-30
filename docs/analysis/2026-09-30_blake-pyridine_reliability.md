# Blake pyridine series: which couplings are reliable (2026-09-30)

- Data: Blake `6_3_26_pyridine_*` processed real spectra (140-200 Hz); fits
  already done (runs/processed/joint_new_ms48, joint_real_staged,
  joint_hj73_neg, joint_hj73_pos; same spectra, range, variant and signal
  weighting: threshold 2.5 sigma, taper 4 Hz).
- Question: which couplings of the joint monotone fits are determined, given
  the other solutions that fit the data about as well?

## Revised criterion (data-only score)

The fit score of fit_joint_series.py is the weighted data residual sum of
squares plus the prior term. The prior centres differ between these fits
(48-start and staged: the earlier centres, e.g. J(A3,HA2) 3.5 Hz; the two
literature starts: Hansen & Jakobsen), so their scores are not comparable.
Every stored solution (all starts of all four fits, 96) was therefore scored
on the data alone: score minus its prior term, recomputed from the solution's
couplings and its fit's prior centres, sigmas and weight. For the best
solution of each fit the data score was also recomputed directly from its
couplings and spectrum parameters; both agree to 1e-6.

    python scripts/reliability_series.py --series series.json \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162.5, "A4": 161.5}}' \
      --couplings "$(cat pyr_lit.json)" --signal-threshold 2.5 --signal-taper 4.0 \
      --fit 48-start=runs/processed/joint_new_ms48/fit.json \
      --fit staged=runs/processed/joint_real_staged/fit.json \
      --fit "lit 4J<0"=runs/processed/joint_hj73_neg/fit.json \
      --fit "lit 4J>0"=runs/processed/joint_hj73_pos/fit.json \
      --literature configs/literature/pyridine_couplings.json --out runs/processed/reliability_data

| fit | solutions | best score = data + prior | best data score | vs best data |
|---|---|---|---|---|
| 48-start | 48 | 0.3084 = 0.2879 + 0.0205 | 0.2788 | 0 |
| staged | 24 | 0.3358 = 0.3004 + 0.0354 | 0.3004 | +7.7 % |
| literature start, 4J < 0 | 12 | 0.3665 = 0.3039 + 0.0626 | 0.3039 | +9.0 % |
| literature start, 4J > 0 | 12 | 0.3648 = 0.3509 + 0.0139 | 0.3509 | +25.8 % |

The best data score belongs to a different 48-start solution than the best
total score (the old reference solution, total 0.308, is +3.3 % on data).
Near-equivalent set: all solutions within 3 % of the best data score, pooled
over the four fits: 5 solutions, all from the 48-start fit (data 0.2788 -
0.2859; total scores 0.312 - 0.323). Spread = largest (max - min) over the
set at any concentration. Reliable: spread <= 1.6 Hz. Trend only: every
solution changes the same way by more than 1 Hz from x 0.02 to 1.00.
Otherwise not determined.

| coupling | class | spread (Hz) | best-data solution, x 0.02 -> 1.00 | literature |
|---|---|---|---|---|
| 1J(C3,H3) | reliable | 0.43 | 164.21 -> 162.79 | 163.04 |
| 1J(C2,H2) | reliable | 0.66 | 178.71 -> 177.29 | 177.63 |
| J(H2,H5) | reliable | 1.26 | 2.92 -> 0.79 | 0.9 |
| J(H2,H4) | reliable | 1.34 | 3.66 -> 5.08 | 1.8 |
| 1J(C4,H4) | trend only | 4.82 | 162.59 -> 157.15 | 162.41 |
| 10 small couplings | trend only | 1.7 - 7.2 | | |
| J(C3,H2), J(H2,H3), J(C2,H4), J(C3,H4) | not determined | 1.9 - 2.9 | | |

Figure: runs/processed/reliability_data/J_trends_reliability.png; numbers
runs/processed/reliability_data/reliability.json; run log RUN_LOG.md there.

## Caveats

- The set holds 5 solutions of one multi-start run whose objective had the
  earlier prior centres, so it describes one region of parameter space. The
  literature-start fits (12 starts each) found no solution within 3 % on the
  data: a search limit or the literature values not fitting these spectra
  with the present model; the data score alone cannot tell.
- With 5 solutions the trend test is weak (5 curves agreeing in direction);
  "trend only" for couplings with spreads of 4-7 Hz means little.
- Reliable means consistent across the set, not accurate: all solutions
  share one forward model, so a missing term shifts them together.
- 1J(C2,H2) and 1J(C3,H3) are reliable under every set tried (also the first
  version below); J(H2,H5) and J(H2,H4) only under this one.

## First version (superseded)

Eleven solutions: the 8 of the 48-start fit within 3 % of its best total
score plus the best of the staged and the two literature-start fits, ranked
by total score (prior included). Reliable: 1J(C2,H2) 1.5 Hz, 1J(C3,H3)
1.5 Hz; trend only: 1J(C4,H4) 4.4 Hz; the other 16 not determined (3.0 -
18.7 Hz). Superseded because the total scores mix different prior centres
and the three added fits are 7.7-25.8 % worse on the data alone. Figure:
runs/processed/J_trends_reliability.png.

## Conclusion

1J(C2,H2) and 1J(C3,H3) are determined (spread < 0.7 Hz, both decreasing by
about 1.4 Hz from x 0.02 to 1.00, within 1.1 Hz of the literature).
1J(C4,H4) decreases in every solution but its value is not fixed. Among the
small couplings only J(H2,H5) and J(H2,H4) agree within 1.6 Hz, and only
within one multi-start run; the rest are not determined by these processed
spectra with the present model.

Commit: see git log (this file).
