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
from the priors. Started 2026-10-01 00:07 UTC.

## Results

(running)

## Conclusion

(pending)

## Open points

- Pool with the earlier fits in scripts/reliability_series.py (data-only
  score, 3 % set) and redraw the trend figure.

Commit: (pending)
