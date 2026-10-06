From: eager-franklin-t
To: all
Re: 2026-10-06T0750Z_eager-franklin-t.md
Subject: handoff from the original developer to the local-machine line (Xuehan's desktop Claude)

Branch: claude/eager-franklin-tp94gf was deleted on 2026-10-06 and recreated on the user's request from
main at 8ab5c84. It is again the working branch of this cloud session (name: eager-franklin-t). Nothing
new on it yet beyond main.

Who I am: the session that wrote most of this repository up to e9a1580:
- the solver and physics;
- the exchange model (D45);
- fit_joint_series with rate families, component search, several molecules and model-line passes (D46);
- band_diagnosis, the J tuner, the fit monitor, the paper figures and J-to-structure (Phase 7);
- the amine analyses.

Where to start: docs/HANDOFF.md (2026-10-05). It has:
- the state of the amines, with run settings, rate-family edges and best couplings;
- how to fit a known structure;
- the open work, in order.

Then read docs/WORKFLOW.md and docs/analysis/2026-10-05_amines_overnight.md.

Things that are NOT in git, which the local line should know:
1. runs/ and my scratchpad were in a cloud container and are gone or will be. Rebuild a result from the
   commands and couplings in the analysis logs. Rebuild the series files with:
   python scripts/make_series_entry.py --fid DATA/<id>-average_fid.npy --id NAME --out runs/series/NAME
   Defaults: crop 0.1 s, record 8 s, exponential 0.3 1/s, zero fill 3. Fit ranges leave out 179.6-180.5 Hz.
2. Helper scripts that lived only in my scratchpad (not in scripts/):
   - snapshot_fit (best monitor vector -> --from-joint file);
   - remap_families (old rates onto new family edges);
   - plot_runs (overlay of runs);
   - eda_species_split (two-species figure);
   - flip_2j (sign test).
   I can port the generic ones into scripts/ with tests if you want them. Tell me here.

Suggested split by area, to keep conflicts rare:
- local line (has the FIDs and a persistent disk): real-data fits under runs/, figures, the experiments
  in HANDOFF "Open work" item 2 (ethylenediamine species 2: tied-rate refit, high-field NMR, CO2 test);
- eager-franklin-t (no persistent data): code in zulf_core/physics (14N spin-1 model, approximate
  intermediate N-H exchange) and fitter features in scripts/fit_joint_series.py, tested on synthetic
  data.
Before touching either area, check status/ and post here.

Conventions: author Xuehan Gao <gao.xh@berkeley.edu> only (AGENTS.md); one branch per line; merge into
main by fast-forward after merging origin/main (docs/DEVELOPMENT.md); claim PLAN items with [~].

Note for the cloud lines: this container's stop hook flags commits whose committer email is not
noreply@anthropic.com as "Unverified" and asks to amend them. That contradicts the author rule and would
rewrite pushed history, so I did not do it; the user is deciding the rule.
