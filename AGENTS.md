# Rules for agents and contributors

This repository builds a simulation-trained model that proposes candidate spin
interpretations for ZULF NMR spectra, and a general forward-model solver that
refines those candidates. Read `docs/ARCHITECTURE.md` before changing any
interface and `docs/PLAN.md` before starting new work.

- All code, comments, docstrings, configuration keys, log messages, figures and
  documentation are in English and ASCII. Discuss with the user in their
  preferred language. `python scripts/check_ascii.py` must pass.
- Nothing about spin count, nucleus order, component count or molecule-specific
  parameter names is a hard-coded constant. Sizes come from configuration and
  the data; fixed-size tensors use masks.
- Physics conventions live in `docs/CONVENTIONS.md` and in
  `zulf_model/physics/protocol.py`. Change them only with an entry in
  `docs/DECISIONS.md` and matching tests.
- Every numerical routine is tested against an independently constructed
  reference (closed forms, brute-force propagation, or exhaustive enumeration),
  not against its own previous output.
- The experimental preprocessing operator (baseline subtraction, crop, FFT
  normalization) is one function used by data generation, training, refinement
  and evaluation. Never add a second implementation of it.
- A refined candidate is a conditional numerical result. Reports and docstrings
  must not describe it as a determined molecular assignment; keep boundary
  hits, residuals, ambiguity and held-out prediction visible.
- Commits and pull requests name Xuehan Gao <gao.xh@berkeley.edu> as the only
  author (set `git config user.name` / `user.email` in every new clone). Do not
  add co-author lines, session links or "generated with" footers to commit
  messages or pull request descriptions. In Claude Code cloud containers a
  stop hook reports these commits as "Unverified" (committer email not
  noreply@anthropic.com, no signature) and asks to reset the author with
  `--amend` / `rebase`. That is expected: keep the author rule, do not change
  the author, and never amend, rebase or force-push commits that are already
  pushed (decision of the owner, 2026-10-06).
- Work on your own branch and merge into `main` as described in
  `docs/DEVELOPMENT.md` (parallel sessions, conflicts, cloud-session limits).
- Experimental inputs are read-only. Generated data, checkpoints and reports go
  under `runs/` or another configured output directory, never into `git`.
- Frozen test sets and generator versions are recorded in `docs/PLAN.md`.
  Failures on a frozen test set are not fed back into training.
- Keep `docs/PLAN.md` current: mark what is done, what was measured, and what
  remains, with the commit that produced each result.
- Every analysis of a spectrum (or of one dataset or series) gets its own
  log file `docs/analysis/YYYY-MM-DD_<sample>_<topic>.md` (template:
  `docs/analysis/TEMPLATE.md`): data, question, exact commands and settings,
  results with numbers, figure paths under `runs/`, conclusion, open points
  and the commit. Write it while the analysis runs and finish it with the
  result. `docs/ANALYSIS_LOG.md` keeps a one-line entry pointing to it. The
  fitting scripts also write a machine run log `RUN_LOG.md` into their
  output directory (`scripts/run_log.py`); link it from the analysis log.
- Keep a running record while analysing data: append to `docs/ANALYSIS_LOG.md`
  what was done, what was observed, the conclusion and the commit, as the work
  happens (not only at the end). When a step teaches something reusable (a
  pitfall, a check, a better setting), update the matching skill under
  `skills/` in the same session.

