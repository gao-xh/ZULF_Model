# Checks before believing a fit

Residual scale
- no-lines residual: background polynomial only, same window, range, mode.
- noise floor: residual expected from noise alone, estimated from
  signal-free parts of the range (rough; it can be wrong where the range
  contains structured background, e.g. below 100 Hz).
- report the explained fraction (null - r) / (null - floor).

Red flags
- Decay rates at the lower bound: very narrow lines chasing background.
- Several parameters at bounds on one component: that component is absorbing
  something else.
- Amplitude ratios that contradict natural abundance or carbon counts.
- Couplings that must be equal differ (symmetry partners, shared protons).
- Implausible sizes: H-H between separate methyl groups of several Hz; 1J and
  large N-H on the same proton; H-H of 20 Hz or more.
- A fit that improves the residual by fitting the edges of a band or the
  low-frequency background instead of the molecular lines.
- A component whose weight collapses to zero: the hypothesis did not need it.

Comparisons
- Same window, range, zero filling, mode (complex or real-only).
- F-type test for nested models: F = (RSS1 - RSS2)/dk / (RSS2/(N - k2)),
  with N reduced for zero filling (correlated points).
- Prefer the model that explains every clear line with consistent physics
  over one with a marginally lower residual.

Model candidates
- Refine them like hypotheses; they usually end near the no-lines residual
  when the models do not resolve the couplings. Report that plainly.
