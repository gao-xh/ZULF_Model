# Overnight batch 2026-10-10: results

Conditional numerical results of the given structures (scripts/run_batch.py, configs/batch_2026-10-10.json);
stage 1: one rate per isotopologue; stage 2: one rate per model line. J in Hz (J_std: linearised).

| id | status | stage 1 score | stage 2 score / residual | field (nT, transverse / z) | couplings (stage 2) |
|---|---|---|---|---|---|
| tert-butylamine | done | 0.3420 | 0.1101 / 0.326 | zero field | J(C1,HC1) 124.47, J(C1,HC3) 4.40, J(C1,HN1) 5.29, J(N1,HC1) 1.66, J(N1,HN1) -75.76 |
| isopropylamine | failed: stage 2 | 0.2603 | - | zero field | J(C1,HC1) 123.63, J(C1,HC2) -5.61, J(C1,HC3) 3.96, J(C1,HN1) 4.98, J(HC1,HC2) 6.10, J(HC2,HN1) 7.46, J(C2,HC1) -4.13, J(C2,HC2) 133.46, J(C2,HN1) -7.01, J(N1,HC1) 0.19, J(N1,HC2) 7.05, J(N1,HN1) -75.80 |
| ethylenediamine | done | 0.1335 | 0.0331 / 0.181 | zero field | J(C1,HC1) 131.20, J(C1,HC2) -3.63, J(C1,HN1) -4.58, J(C1,HN2) 3.15, J(HC1,HC2) 8.56, J(HC1,HN1) 5.60, J(N1,HC1) 7.00, J(N1,HC2) 4.92, J(N1,HN1) -65.00 |
| diethylamine | done | 0.9192 | 0.5349 / 0.649 | zero field | J(C1,HC1) 124.81, J(C1,HC2) -2.69, J(C1,HN1) 3.80, J(HC1,HC2) 8.40, J(HC2,HN1) 0.19, J(C2,HC1) 1.52, J(C2,HC2) 132.29, J(C2,HN1) -11.58, J(C2,HC3) 11.57, J(N1,HC1) 0.20, J(N1,HC2) 3.99, J(N1,HN1) -75.89 |
