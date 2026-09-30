# Blake pyridine series: which couplings are reliable (2026-09-30)

- Data: Blake `6_3_26_pyridine_*` processed real spectra; fits already done
  (runs/processed/joint_new_ms48, joint_real_staged, joint_hj73_neg, joint_hj73_pos).
- Question: which couplings of the best joint fit (48-start, score 0.308) are
  determined, given the other near-equivalent solutions?

## Setup

Eleven near-equivalent solutions: the 8 solutions of the 48-start fit within 3 %
of its best score, plus the staged (0.336) and the two literature-start fits
(0.366, 0.365). Spread = largest range over the six concentrations of the 11
solutions. Classes: reliable, spread <= 1.6 Hz; trend only, every solution
changes in the same direction by more than 1 Hz from x 0.02 to 1.00; otherwise
not determined. (Linearised standard errors of the best fit are 0.01-0.4 Hz and
describe one minimum only.)

## Results

| coupling | class | spread (Hz) | 48-start best, x 0.02 -> 1.00 |
|---|---|---|---|
| 1J(C2,H2) | reliable | 1.5 | 178.75 -> 177.33 (literature 177.63) |
| 1J(C3,H3) | reliable | 1.5 | 164.24 -> 162.80 (literature 163.04) |
| 1J(C4,H4) | trend only (decreasing) | 4.4 | 162.49 -> 156.99 (literature 162.41) |
| other 16 (H-H and long-range C-H) | not determined | 3.0-18.7 | |

Figure: runs/processed/J_trends_reliability.png.

## Conclusion

Only 1J(C2,H2) and 1J(C3,H3) (values within ~0.75 Hz, both decreasing by ~1.4 Hz
from x 0.02 to 1.00) are supported by every near-equivalent solution; 1J(C4,H4)
decreases in every solution but its value is not fixed (and is low at x = 1.00
against the literature). The small couplings are not determined by these
processed spectra with the present model.

Commit: see git log (this file).
