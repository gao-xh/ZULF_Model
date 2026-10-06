From: eager-franklin-t
To: xuehan, all
Re: 2026-10-06T0815Z_eager-franklin-t.md
Subject: merged into main at dfda00a: fitter features and run helpers (full test suite passed)

main fast-forwarded from 68a505f to dfda00a, on Xuehan's request after a full test run.
- Full suite on 70ac144 (the fitter work): all 25 modules OK, test_hypothesis 1070 s.
- dfda00a merges that commit with main 68a505f (your data commits; no file overlap). Rechecked:
  check_ascii, test_diagnostics, test_processing, test_j_tuner and test_joint_series, 48 tests OK.

New on main (interfaces as announced in my previous message):
- fit.json "family_edges_hz"; rates remapped across edges in --from-joint and j_tuner.load_fit;
- --tie-rates REGEX; --family-edges auto; --component-search-hold-gains;
- scripts/snapshot_fit.py, scripts/run_problem.py, scripts/plot_runs.py, scripts/plot_components.py.

For the local line: --tie-rates is available now for the ethylenediamine species-2 test, e.g.
  ... eda_2sp_fam198 options ... --tie-rates '^P2:' --from-joint <eda_2sp_fam198 fit.json or a rebuilt start>
Figures: python scripts/plot_components.py runs/processed/NAME OUT.png --band 176,216

Not on main yet: 14N quadrupolar relaxation in zulf_core/physics/exchange.py (f71cdab, on my branch, its
tests pass). The model and fitter integration comes next; I will announce it before merging.
