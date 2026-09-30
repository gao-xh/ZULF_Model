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

- C (N-H dropped, delay -30 ... +30 ms): data-core residual 0.181
  (stage 1: 0.235), k 25, delay -1.35 ms, no boundary hits,
  141 s; 1J 190.1 / 167.5 / 171.1 Hz.
  Same model as the earlier decoupled fit (0.195, delay 3.7 ms, bounds +-10 ms): the wider delay bounds alone
  lower the residual. The earlier exchange fit (0.171) had its delay at the +10 ms bound.
- A (k free) and B (k = 1e4 1/s): first attempt lost in a container restart, second one stopped with its
  background task after ~30 min (time limit); third attempt started detached (setsid nohup) at
  2026-09-30 00:27 UTC. Lesson: long fits must run detached, and fit_staged writes nothing before the end.

## Conclusion

(pending)
