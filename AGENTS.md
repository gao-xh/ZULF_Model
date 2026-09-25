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
- Experimental inputs are read-only. Generated data, checkpoints and reports go
  under `runs/` or another configured output directory, never into `git`.
- Frozen test sets and generator versions are recorded in `docs/PLAN.md`.
  Failures on a frozen test set are not fed back into training.
- Keep `docs/PLAN.md` current: mark what is done, what was measured, and what
  remains, with the commit that produced each result.
