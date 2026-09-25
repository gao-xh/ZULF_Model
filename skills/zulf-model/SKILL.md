---
name: zulf-model
description: Propose and refine spin interpretations (isotopologues, magnetic-equivalence groups, J matrices) for zero- to ultralow-field NMR spectra with the ZULF_Model toolkit. Use for simulating ZULF transitions and spectra, diagnosing and processing averaged experimental FIDs, generating simulated training data, training or running the candidate models, refining candidates against data with held-out validation, and identifiability checks.
---

# ZULF_Model agent guide

The toolkit is a backend of typed tools. Call them; do not re-implement physics
or processing in chat, and never invent numbers that a tool did not return.

## Access

- JSON CLI: `zulf-model tool NAME --json '{...}'`, background jobs with
  `zulf-model submit NAME --json '{...}'` and `zulf-model tool get_job --json '{"job_id": "..."}'`.
- MCP (Claude Code, Claude Desktop, Codex): `python -m zulf_model.agent.mcp_server`.
- API tool definitions: `zulf-model tools export --format anthropic` or `--format openai`.
- Outputs go to `$ZULF_MODEL_WORKSPACE` (default `./runs`). Experimental folders are read-only.

## Standard workflow for an experimental dataset

1. `describe_project` once per session.
2. `diagnose_fid` on the averaged FID. Report the first-point anomaly, saturation
   plateau, ringing end and baseline fits. The candidate recipes are hypotheses.
3. `process_fid` with two or three candidate recipes; inspect the figures.
4. Obtain candidates: `propose_candidates` with a trained checkpoint, or an
   explicit hypothesis written as an interpretation.
5. `refine_candidates` (use `submit_job` for long records) with held-out group
   averages whenever they exist. Keep all branches. Use `settings.ties` for
   couplings shared across isotopologues and `settings.policy.nuisance` for
   explicit instrument baseline terms.
6. Compare recipes: a parameter that moves between recipes beyond its tolerance is
   processing-sensitive and must be reported as such.
7. Optional: `identifiability` for the chosen candidate.

## What to report

- Every result is a conditional numerical candidate. Quote flags
  (`search_boundary`, `optimizer_not_converged`, `search_budget_exhausted`,
  `continuation_unavailable_without_fid`, `provisional_ranking`), band residuals
  and frozen held-out residuals next to any J value.
- A decay rate at its bound, or a J at its bound, is not a measurement.
- A global sign flip of all J is unobservable; relative signs are observable.
- Do not describe a candidate as the molecular structure or assignment.

## Training

`train_model` with a run configuration (see `configs/run_*.json`). Use
`submit_job` and poll. Device selection is automatic (CUDA, Apple MPS, CPU).
Never feed frozen test failures back into training.
