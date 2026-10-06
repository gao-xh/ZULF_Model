# eager-franklin-t

- Kind: Claude Code cloud session (no FIDs, no persistent disk)
- Branch: claude/eager-franklin-tp94gf
- Task (assigned by Xuehan 2026-10-06): 1. fitter features, 2. then the 14N spin-1 model (PLAN Phase 3e, claimed)
- Area: scripts/fit_joint_series.py (rate ties, automatic family edges, rate remap in --from-joint, local fit
  with gains held), new scripts/snapshot_fit.py, tests/test_j_tuner.py; later zulf_core/physics (14N) and tests
- State: working (fitter features done, waiting for merge go-ahead; porting plot helpers)
- Last update: 2026-10-06T0821Z
- Next step: port plot_runs and a components figure into scripts/, then 14N
- Questions: none
