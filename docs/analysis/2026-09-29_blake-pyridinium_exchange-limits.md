# Blake pyridinium (6 M HCl): exchange rate vs its limits (2026-09-29)

- Data: Blake `6_14_26_pyridinium_6M` processed real spectrum (140-200 Hz).
- Question: is the fitted N-H exchange rate (78 1/s, previous analysis
  2026-09-29_blake-pyridinium_exchange-fit.md, delay at its 10 ms bound)
  determined by the data, or does it absorb other model errors? Physical
  expectation (estimate): k ~ 1e4-6e4 1/s, i.e. the fast limit.
- Code: see commit below; run logs runs/processed/pyridinium_6M_{kfree,k1e4,nh_dropped}/RUN_LOG.md.

## Setup

Identical settings for three fits: small couplings started at the pyridine
literature (sigma 1 Hz H-H / 2 Hz C-H, weight 10), 1J started at 187.5 / 169.0 /
171.2 Hz, staged (H-H held first), 2 stage-2 starts, spread 0.5 Hz, delay bounds
widened to -30 ... +30 ms.

- A (kfree): protonated, HA1 exchanging, k free (starts 10, 1, 100, 3000 1/s).
- B (k1e4): protonated, HA1 exchanging, k fixed at 1e4 1/s.
- C (nh_dropped): N-H proton not in the spin system (fast limit).

## Results

(running)

## Conclusion

(pending)
