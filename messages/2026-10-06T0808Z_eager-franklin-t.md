From: eager-franklin-t
To: xuehan
Re: -
Subject: claimed: fitter features, then 14N model (area: scripts/fit_joint_series.py, later zulf_core/physics)

Xuehan assigned both to me in the session chat. Claimed in docs/PLAN.md on claude/eager-franklin-tp94gf.
Area until further notice: scripts/fit_joint_series.py, new scripts/snapshot_fit.py, tests/test_j_tuner.py;
afterwards zulf_core/physics for 14N. Please avoid those files on xuehan/local, or post here first.

Interface changes (announced before merging into main):
- fit.json gains a "family_edges_hz" key;
- --from-joint remaps decay rates when the previous fit used other family edges;
- new option --tie-rates REGEX: one decay rate for every rate family of the matching components.
Old fit.json files (no edges recorded) load as before.
