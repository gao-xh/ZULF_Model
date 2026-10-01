# Blake pyridine series: 48 more starts from the near-best set (2026-10-01)

- Data: Blake `6_3_26_pyridine_{02,33,50,66,75,100}` processed real spectra
  (140-200 Hz), x = pyridine mole fraction 0.02 ... 1.00 (not committed).
- Question: the near-equivalent set of the reliability analysis
  (2026-09-30_blake-pyridine_reliability.md) holds only 5 solutions, all from
  one 48-start run. Do more starts around them find better data minima, and
  which couplings stay within 1.6 Hz in a larger set?
- Code: fit_joint_series.py with the new `--seeds` option (several start
  tables; perturbed starts spread over them in turn); run log
  runs/processed/joint_ms48_cont/RUN_LOG.md.

## Setup

Same model and settings as the 48-start fit (joint_new_ms48): pyridine ring,
13C2 / 13C3 / 13C4 isotopologues, monotone J(x), earlier prior centres
(sigma 0.5 Hz H-H, 1.5 Hz C-H, weight 10), signal weighting threshold 2.5,
taper 4 Hz. Seeds: the 5 solutions within 3 % of the best data score (data
scores 0.2788-0.2859); spectrum parameters from the 48-start best.

    python scripts/fit_joint_series.py --series series.json \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162.5, "A4": 161.5}}' \
      --couplings "$(cat pyr_lit.json)" --prior-sigma-hh 0.5 --prior-sigma-ch 1.5 --prior-weight 10 \
      --signal-threshold 2.5 --signal-taper 4.0 --from-joint runs/processed/joint_new_ms48/fit.json \
      --seeds seeds_ms48_set.json --starts 40 --spread 1.0 --prior-starts 8 --seed 7 --workers 4 \
      --out runs/processed/joint_ms48_cont

48 starts: the 5 seeds, 35 perturbed (7 per seed, level spread 1 Hz), 8 drawn
from the priors. Started 2026-10-01 00:07 UTC, finished 00:44 UTC.

## Results

48 starts in 2228 s (4 workers). Best fit score 0.2920 (data 0.2659 + prior
0.0261), against 0.3084 for the first 48-start run; next total scores
0.3011, 0.3014, 0.3092.

Pooled with the earlier fits on the data-only score
(scripts/reliability_series.py, 5 fits, 144 solutions; out
runs/processed/reliability_data_v2):

| fit | best data score | vs best |
|---|---|---|
| 48 more (this run) | 0.2634 | 0 |
| 48-start | 0.2788 | +5.9 % |
| staged | 0.3004 | +14.0 % |
| literature start, 4J < 0 | 0.3039 | +15.4 % |
| literature start, 4J > 0 | 0.3509 | +33.2 % |

Set within 3 %: 4 solutions, all from this run (data 0.2634 - 0.2693); the
earlier 5-solution set is now 5.9 % above the best.

| coupling | class | spread (Hz) | best-data solution, x 0.02 -> 1.00 | literature |
|---|---|---|---|---|
| 1J(C3,H3) | reliable | 0.22 | 164.24 -> 162.79 | 163.04 |
| 1J(C2,H2) | reliable | 0.73 | 178.73 -> 177.26 | 177.63 |
| J(H2,H5) | reliable | 0.95 | 3.43 -> 0.79 | 0.9 |
| J(C3,H6) | reliable | 1.58 | 2.42 -> 4.64 | 1.65 |
| 1J(C4,H4) | trend only | 3.02 | 162.21 -> 157.15 | 162.41 |
| 8 other small couplings | trend only | 2.2 - 5.8 | | |
| J(C3,H2), J(H2,H4), J(H2,H3), J(H3,H4), J(C3,H4), J(C2,H4) | not determined | 2.1 - 4.4 | | |

Best data-score solution, spectra refitted (decay rates and delays, couplings
held; data score 0.26343 = stored): relative residuals 0.211 / 0.220 /
0.205 / 0.234 / 0.230 / 0.177 (x 0.02 ... 1.00).

Figures: runs/processed/reliability_data_v2/J_trends_reliability.png,
runs/processed/reliability_data_v2/best_data_spectra.png,
runs/processed/joint_ms48_cont/spectra.png.

## Conclusion

More starts keep finding lower data minima (best data score 0.2879 ->
0.2788 -> 0.2634 over three runs), so the search has not converged and every
"near-equivalent set" so far is the neighbourhood of the latest best. 1J(C2,H2)
and 1J(C3,H3) are again reliable (spread 0.7 and 0.2 Hz, within 0.4 Hz of the
literature at x = 1.00) and match every earlier set; J(H2,H5) stays within
1 Hz (0.79 Hz at x = 1.00, literature 0.9). J(H2,H4) dropped out and J(C3,H6)
entered at 1.58 Hz, just under the threshold: the small-coupling classes are
not stable between rounds. Many small couplings change by several hertz in one
step between x 0.66 and 0.75, sometimes through zero (e.g. J(C2,H3) 5.5 ->
-3.2 Hz); such jumps are not plausible chemistry and point to the small
couplings absorbing something the model lacks (line shape, processing,
field-switch timing).

## Follow-up: unfitted small peaks and the field-switch sequence

The best data-only solution leaves peaks of 10-26 % of the maximum unfitted
(runs/processed/reliability_data_v2/missing_peaks.png): one about 1.5 Hz above
the strongest line that moves with it (169.9 / 169.5 / 169.1 Hz at x 0.33 /
0.50 / 0.66; the model has only a weak 13C2 line there, 0.09 of the maximum,
0.3 Hz away), one near 151 Hz (model lines 0.04), and features at 163-165 Hz.
Model lines exist near most of them; their intensities are too low.

Diagnostic (scratchpad hold_diag.py, couplings held at the best data-only
solution, decay rates and delay refitted per case): initial state from the
Ajoy-lab sequence (80 uT along x, GF1 ramped to 0 over 30 ms with 1 uT along
-y, hold, sudden switch-off; brute-force propagation as in
scripts/field_switch_check.py) instead of the sudden drop. Data score change
per spectrum (x 0.33 / 0.50 / 0.66): ramp without hold +5 / +10 / +5 %;
hold 1 ms +42 / +62 / +38 %; 3 ms +4 / +2 / +3 %; 10 ms +18 / +21 / +14 %;
30 ms +83 / +106 / +77 %; 100 ms +95 / +125 / +69 %; long hold dephased at
1 uT -2.4 / -3.5 / +1.9 %. No hold time raises the 169.5 Hz peak without
spoiling the rest. With these couplings the field-switch sequence does not
explain the missing peaks (caveat: the couplings were fitted with the sudden
drop; a joint refit with the hold model is not done).

Model peak height near 169.5 Hz (fraction of the maximum; scratchpad
peak_compare.py, each fit's stored best, windows 169.5-170.3 / 169.1-169.9 /
168.7-169.5 Hz at x 0.33 / 0.50 / 0.66):

| | data | 48-start best (total 0.308) | 48 more best (data-only) | staged | lit 4J<0 | lit 4J>0 |
|---|---|---|---|---|---|---|
| x 0.33 | 0.350 | 0.285 | 0.158 | 0.242 | 0.121 | 0.182 |
| x 0.50 | 0.309 | 0.252 | 0.127 | 0.195 | 0.092 | 0.142 |
| x 0.66 | 0.195 | 0.166 | 0.157 | 0.162 | 0.070 | 0.152 |

The earlier 48-start best (lowest total score) reproduces this peak much
better; the lower data score of the new best is bought elsewhere at the cost
of this peak. The data-only ranking alone is therefore not a safe choice of
reference solution; the peak should be checked in any candidate.

## Open points

- The set holds 4 solutions of one run; the trend test stays weak.
- Inspect the x 0.66 / 0.75 spectra for a change in acquisition or processing.
- Unfitted peak 1.5 Hz above the strongest line: test a free extra line (width, height, x dependence) and per-line widths; 14N effects on 13C2.
- Further rounds would likely lower the score again; a global search or the
  missing model terms are more useful than more local starts.

Commit: (pending)
