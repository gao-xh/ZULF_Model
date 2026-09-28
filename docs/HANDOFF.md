# Handoff (2026-09-28)

State of the project at handoff to the team workspace, what is in flight,
and how to continue. Read `AGENTS.md` (rules), `docs/ARCHITECTURE.md`
(interfaces) and `docs/PLAN.md` (ledger) first; `docs/ANALYSIS_LOG.md` is
the running lab notebook, `docs/DECISIONS.md` (D1-D44) the reasons.

## Repository

- GitHub `gao-xh/ZULF_Model`, branch `claude/eager-franklin-tp94gf`. This is
  currently the only branch on the remote (no `main` yet); create one from it
  or keep working on it.
- Packages:
  - `zulf_core`: NumPy/SciPy physics, rendering, the single processing
    operator (`render/acquisition.py`), solver.
  - `zulf_model`: PyTorch models, generator, training, agent tools.
  - `zulf_hypothesis`: programmatic hypotheses, blind search,
    known-structure fits, reports.
  - `zulf_processing`: experimental processing per dataset (D44).
- Skills for agents under `skills/` (blind analysis, FID processing,
  phasing, hypothesis refinement, spectrum interpretation, model).

## Setup and checks

```bash
pip install -e ".[dev]"
python scripts/check_ascii.py                      # must pass (code and docs ASCII English)
python -m unittest discover -s tests               # 184 tests, about 46 min on 4 cores
python -m unittest tests.test_processing           # processing only, about 1 min
```

## Data (not in git)

Experimental FIDs are read-only and never committed (`*.npy` and `runs/` are
ignored). The six confirmed samples (`configs/confirmed_samples.json`) are
the regression set; put the files in one folder and set `ZULF_DATA_DIR`:

| id | file | compound |
|---|---|---|
| e3d282da | e3d282da-average_of_two_batch_fids.npy | L-alanine |
| b683220d | b683220d-average_fid.npy | lactic acid |
| 7ad4aafb | 7ad4aafb-average_fid.npy | pyridine |
| fde3fbb2 | fde3fbb2-average_fid.npy | ethylenediamine |
| 4322bdfc | 4322bdfc-average_fid.npy | triethylamine |
| e66a4b08 | e66a4b08-average_fid.npy | N-ethylmethylamine |

NMRduino, 4 kHz, 65516 points (16.4 s), sequence
standard_zf_4000Hz_no_dead.seq. Source inventory: PLAN "Experimental data
inventory".

## Reference results to regenerate (runs/ is not in git)

`scripts/validate_processing.py` reads the known-structure fits from
`runs/regression/final/summary.json`. Regenerate before using it:

```bash
ZULF_DATA_DIR=... python scripts/regression_confirmed.py --mode known --out runs/regression/final   # about 1 h, 4 workers
ZULF_DATA_DIR=... python scripts/regression_confirmed.py --mode blind --out runs/regression/blind    # about 1 h
ZULF_DATA_DIR=... python scripts/validate_processing.py                                             # about 7 s
```

Expected (ANALYSIS_LOG, "Final regression"): blind search right skeleton
at rank 1 on 6 of 6 (two extension rounds); known-structure reduced chi2
alanine 17.8, lactic acid 18.8, pyridine 70.6, ethylenediamine 680,
triethylamine 115, N-ethylmethylamine 196.

## New sample

```bash
python scripts/analyze_sample.py FID.npy --id ID                    # blind: processing, hypotheses, search, report
python scripts/analyze_sample.py FID.npy --id ID --structure '{...}' # known structure: every exchange regime and variant
```

Outputs in `runs/blind/ID`: processing.{png,json} (plan, edge, phase),
overview.png, blind or structure .{json,md,png} (ranked table, J matrices,
phased fit figure). Procedure and pitfalls: `skills/zulf-blind-analysis`.

## Processing (zulf_processing): state

- Per-dataset plan with reasons: crop start after this dataset's
  ringing; the other parameters are defaults (SG 201/2, window 1 s,
  apodization 0.6 1/s, zero fill 4).
- Delay: measured per dataset from the switching edge in the raw FID. The
  edge is at 3.41-3.51 ms; the delay is minus the edge time.
- Phase0: instrument calibration by default (`phase_calibration` in the
  config, 176.3 deg; delay = edge - 0.033 ms). Leave-one-out it is within
  14 deg of each sample's own complex fit.
- Model-free criteria (entropy, lines; global search plus fine-tune) err
  2-86 deg on J-spectra; they are kept for comparison and validated by
  `validate_processing.py`.
- Complex fits fit their own phase and are not affected; the processed
  phase is used by the phased route and figures.
- If the pulse sequence, hardware or processing defaults change, refresh
  the calibration (skill `zulf-phasing`).

## In flight

1. **Window length per dataset.** The records are 16.4 s, but the fixed
   window uses only 1 s.
   - `signal_extent` gives each dataset's signal end: 1.25-4.75 s. The
     1 s window misses 11-14 % of the signal energy on four samples.
   - Plan `window_mode="signal_extent"` sets the window from it, and
     zero fill keeps 4 points per Hz. The default is still "fixed".
   - Known fits with the long window (`regression_confirmed.py
     --processing dataset --window-mode signal_extent`):
     - Lactic acid and triethylamine: 1J stable, fitted delay at the
       edge, phased-route chi2 halved; the triethylamine best model became
       the protonated, slow-exchange form.
     - N-ethylmethylamine: mixed (delay moved 0.15 ms off the edge, one 1J
       moved 1.35 Hz).
   - Cost is 4-6x (render and SG scale with record length).
   - Numbers: ANALYSIS_LOG "Window length per dataset".
2. **Soft crop and weighted SG** (user proposal, not started). Replace
   the hard crop with a window that is flat in the middle and Gaussian at
   both ends, and give the SG kernel the same shape (weighted local
   polynomial). Design agreed so far:
   - Put it in `zulf_core.render.acquisition` (the single operator) as new
     `Acquisition` fields, for example `taper_start_s`, `taper_stop_s` and
     `sg_taper`. Default off, so existing results are unchanged.
   - `savgol_baseline` and the renderer must use the same (weighted) SG
     coefficients.
   - The analytic renderer folds exponential apodization into the line
     damping (closed form). Keep the closed form for the flat part and sum
     the Gaussian edges (about 3 sigma long) numerically.
   - Tests against independent references: pointwise multiplication for
     the window, brute-force weighted least squares for the SG.
   - Expected use: the stop taper replaces exponential apodization, which
     broadens every line by 0.19 Hz. A short start taper (a few ms)
     suppresses residual ringing and the hard-edge step, and may allow an
     earlier crop start.
   - The model applies the same operator, so these choices do not bias
     parameters. Choose them by fit quality: J stability across settings,
     uncertainties, c_hat, residual whiteness.
3. **Renderer cost for long records**, needed before the long window can
   become the default.

## Next steps (suggested order)

1. Renderer / SG speed for long records.
2. Soft crop and weighted SG in the operator (item 2 above), with tests.
3. Re-run all six known fits and the blind regression with signal-extent
   windows and candidate tapers. Change the defaults only if J values stay
   stable and the fits improve; record the results in ANALYSIS_LOG.
4. SG window per dataset from the lowest signal band; instrument-line
   detection per dataset.
5. Open items in PLAN:
   - zero-order phase from the data alone (e.g. a reference signal
     recorded with every acquisition);
   - known-structure budget for 11-spin exchange forms;
   - the physics double diagonalisation;
   - intermediate exchange;
   - the low-frequency range (a 48-51 Hz feature is not modelled).

## Working rules learned the hard way

- Never kill processes by command-line pattern (`pkill -f`). It can kill
  the calling shell; use PIDs.
- Long runs (regressions, the full test suite) belong in the background
  with logs under `runs/`. Do not run two heavy jobs on 4 cores if the
  timings matter.
- A processing step that takes more than seconds is a bug. The unused
  exponential fits of `diagnose_fid` once made every dataset take 546 s.
- Every numerical routine gets a test against an independent reference.
  Record what was done and observed in ANALYSIS_LOG as you go, and update
  the matching skill when a step teaches something reusable.
- Results are conditional numerical candidates: keep flags, residuals,
  boundary hits and alternatives visible in reports.
