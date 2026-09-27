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
  In a free variant of one molecule this means the fit describes a
  different set of isotopologues (check `collapsed_component`, < 0.2 of the
  strongest per unit abundance); the search ranks such fits last (D42).
  4322bdfc: a CH-CH3 fit with its CH at 0.061 beat the correct ethyl.
- A component that turns into background: much larger amplitude than its
  partners together with a high decay rate (e3d282da: C-alpha 7.5x at
  7.7 1/s). Its broad lines can cut holes (139 Hz) that cost little outside
  the cores. Fix with natural ratios and a tied rate.
- Decay rates that differ strongly between isotopologues of one molecule
  (b683220d H7: 4.7 vs 1.9 1/s): the faster one hides unresolved couplings.
- Couplings that differ between equally good starts are not determined;
  report them as such (e3d282da small couplings, b683220d extra spin).

- An extension (a nested model) that scores worse than its parent: suspect
  the start, not the physics. A parent that compensated for missing
  couplings (broadened lines, shifted small couplings) traps a warm start;
  refine extensions from both the parent's optimum and the proposal values.

Comparisons
- Same window, range, zero filling, mode (complex or real-only).
- With signal weighting compare `data_region_residual` (data cores only);
  `signal_region_residual` includes each model's own cores and differs
  between models.
- Plot the weights actually used (data cores and model cores) under the
  fits; check a weighting claim on the fit itself before stating it.
- F-type test for nested models: F = (RSS1 - RSS2)/dk / (RSS2/(N - k2)),
  with N reduced for zero filling (correlated points).
- Prefer the model that explains every clear line with consistent physics
  over one with a marginally lower residual.

Model candidates
- Refine them like hypotheses; they usually end near the no-lines residual
  when the models do not resolve the couplings. Report that plainly.
