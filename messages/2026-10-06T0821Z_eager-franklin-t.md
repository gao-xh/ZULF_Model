From: eager-franklin-t
To: xuehan, all
Re: 2026-10-06T0812Z_xuehan.md
Subject: fitter features done on claude/eager-franklin-tp94gf (ff0d903); interface changes before the merge into main

All on claude/eager-franklin-tp94gf, merged with main e222fb1, tests/test_j_tuner.py 18 OK.

Interface changes in scripts/fit_joint_series.py:
- fit.json gains "family_edges_hz" (the edges used).
- --from-joint and j_tuner.load_fit remap the decay rates when the previous fit.json has other edges: each
  new family takes the rate of the old family containing its centre. Old fit.json files without edges load
  as before (cbf7655).
- --tie-rates REGEX on component labels: every rate family of a matching component follows its family 0.
  For ethylenediamine species 2: --tie-rates '^P2:'. The start value is family 0's rate from --from-joint,
  which for eda_2sp_fam198 is the broadened 9/s family. Give a start file with a sensible P2 rate (or
  accept the refit from there) (8caf49b).
- --family-edges auto: edges from the start vector's line clusters and the sharp data peaks; the numbers
  are written into RUN_LOG / monitor argv and fit.json (8f5a09a).
- --component-search-hold-gains / JointSeries.local_fit(hold_gains=True): window cost = window part of the
  whole objective (04b9df6).
- new scripts/snapshot_fit.py RUN OUT.json: best monitor vector (rescored on the plain objective) as a
  --from-joint file (64bcc2e).

Suggested real-data checks for the local line once this is on main:
1. ethylenediamine: eda_2sp_fam198 settings plus --tie-rates '^P2:' (one rate for species 2). Does the
   -18 Hz geminal survive?
2. --family-edges auto against the hand edges of tea_fam121, nema_fam121_final and eda_2sp_fam198 (same
   starts): objective and couplings.

Next on my line: porting plot_runs (and a components figure) into scripts/, then 14N in zulf_core/physics.
Merge into main: waiting for Xuehan's go-ahead in my session chat.
