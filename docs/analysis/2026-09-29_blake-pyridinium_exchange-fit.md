# Blake pyridinium (6 M HCl): fitted N-H exchange rate (2026-09-29)

- Data: Blake `6_14_26_pyridinium_6M` processed real spectrum (140-200 Hz).
- Question: does an intermediate N-H exchange rate describe the spectrum better
  than the decoupled (fast) limit?
- Code: cc92e12 (exchange model, D45); run 2026-09-29 19:26-21:05 UTC.

## Setup

Protonated pyridine ring (N-H proton HA1 kept), HA1 exchanging with the solvent
at a fitted rate (one k tied across the three isotopologues); small couplings
started at the pyridine literature with priors sigma 1 Hz H-H / 2 Hz C-H,
weight 10 (no pyridinium literature at hand); 1J started at the earlier
pyridinium fit 187.5 / 169.0 / 171.2 Hz; k starts 10, 1, 100, 3000 1/s.

    python scripts/fit_staged.py --freq .../6_14_26_pyridinium_6M_frequency.npy \
      --values .../6_14_26_pyridinium_6M_amplitude.npy --id pyridinium_6M \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 187.5, "A3": 169.0, "A4": 171.2}}' \
      --protonate --couplings "$(cat pyr_lit_hj73_neg.json)" --prior-sigma-hh 1.0 --prior-sigma-ch 2.0 \
      --prior-weight 10 --exchange HA1 --kex-start 10 --starts 2 --spread 0.5 \
      --out runs/processed/pyridinium_6M_exchange

## Results

| model | k | data-core residual | delay |
|---|---|---|---|
| N-H decoupled (earlier, runs/processed/pyridinium_6M_staged) | 25 | 0.195 | 3.7 ms |
| N-H exchanging, fitted k | 27 | 0.171 (stage 1: 0.172) | 10 ms (bound) |

k_ex = 78 1/s; 1J 190.7 / 163.7 / 171.3 Hz (C3 moved from 169.0);
J(A2,HA1) 4.9, J(HA2,HA1) 4.1 Hz. Boundary hits: phase delay (10 ms), C2
isotopologue rate (0.05 1/s); flags optimizer_not_converged,
search_budget_exhausted. The model shows broad lines near 154 and 198 Hz that
the data do not.

Figures: runs/processed/pyridinium_6M_exchange/fit.png,
runs/processed/pyridinium_exchange_vs_decoupled.png.

## Conclusion

The exchange model fits better (182 Hz peak reproduced), but the value of k is
conditional: two parameters sit at bounds. Expected physically (estimate, not a
measured literature value): dissociation to water k ~ k_r Ka ~ 1e10 x 10^-5.2
~ 6e4 1/s (1e4 in 6 M HCl), plus 14N scalar relaxation, i.e. the fast limit;
78 1/s more likely absorbs other model errors.

## Open points

- Same settings with wider delay bounds: k free vs k fixed at 1e4 1/s vs N-H
  dropped.

Commit: e9bc834
