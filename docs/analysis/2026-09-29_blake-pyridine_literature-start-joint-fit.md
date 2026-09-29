# Blake pyridine series: joint fits started from the literature (2026-09-29)

- Data: Blake `6_3_26_pyridine_{02,33,50,66,75,100}` processed real spectra
  (140-200 Hz), x = pyridine mole fraction 0.02 ... 1.00.
- Question: starting from (and with priors at) the literature couplings, do the
  joint monotone fits approach the literature, and does 1J(C4) stay low?
- Code: af60f59 (literature file), runs from 2026-09-29 18:36-19:26 UTC.

## Setup

Literature: Hansen & Jakobsen, JMR 10, 74 (1973), via Bax, JMR 52, 330 (1983),
Table 1 (absolute values, 70 % v/v acetone-d6/TMS); H-H from Castellano et al.,
JCP 46, 327 (1967). Stored in `configs/literature/pyridine_couplings.json`.
The earlier prior centres had J(A3,HA2) = 3.5 Hz (literature 8.47).
Two sign variants for 4J(A2,HA5), 4J(A3,HA6): negative, positive.

    python scripts/fit_joint_series.py --series series.json \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 177.6, "A3": 163.0, "A4": 162.4}}' \
      --couplings "$(cat pyr_lit_hj73_{neg,pos}.json)" --prior-sigma-hh 0.5 --prior-sigma-ch 1.5 \
      --prior-weight 10 --signal-threshold 2.5 --signal-taper 4.0 --starts 6 --spread 0.5 \
      --prior-starts 6 --seed 17 --workers 2 --hold-small-first --out runs/processed/joint_hj73_{neg,pos}

## Results

| fit | score | data-core residuals (x 0.02 ... 1.00) |
|---|---|---|
| literature start, 4J < 0 | 0.366 | 0.250 0.240 0.220 0.239 0.222 0.206 |
| literature start, 4J > 0 | 0.365 | 0.279 0.248 0.244 0.275 0.240 0.214 |
| earlier 48-start (old priors) | 0.308 | 0.228 0.229 0.225 0.254 0.231 0.171 |

(scores include the prior term, whose centres changed.)

| coupling | literature | 4J<0, x 0.02 -> 1.00 | 4J>0, x 0.02 -> 1.00 |
|---|---|---|---|
| J(A2,HA2) | 177.63 | 178.74 -> 176.59 | 178.73 -> 178.07 |
| J(A3,HA3) | 163.04 | 163.49 -> 161.52 | 165.00 -> 162.54 |
| J(A4,HA4) | 162.41 | 160.63 -> 157.79 | 161.42 -> 157.98 |
| J(A3,HA2) | 8.47 | 6.25 -> 6.22 | 6.32 -> 5.42 |
| J(A3,HA4) | 0.84 | 5.59 -> 7.91 | 3.21 -> -3.70 |
| J(A4,HA3) | 0.70 | 6.36 -> 3.86 | 3.65 -> 1.33 |
| J(A4,HA2) | 6.34 | 3.90 -> -2.57 | 4.75 -> 1.19 |
| J(A2,HA6) | 11.16 | 14.38 -> 16.75 | 11.02 -> 12.77 |

Figures: runs/processed/joint_hj73_{neg,pos}/spectra.png,
runs/processed/J_literature_start.png.

## Conclusion

1J(C2,H2) and 1J(C3,H3) agree with the literature and decrease with the pyridine
fraction. 1J(C4,H4) at x = 1.00 is ~158 Hz in every fit (literature 162.2-162.4):
not caused by the priors. Several small couplings move Hz away from the
literature whatever the start; the sign of 4J is not decided. A missing model
term (field-switch hold, 14N broadening of C2, unknown processing) is more
likely than a search problem.

## Open points

- Fitted hold time at 1 uT; raw FIDs and sequence timings from Blake.

Commit: f3d1ba9
