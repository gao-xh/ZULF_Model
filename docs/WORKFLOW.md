# Fitting workflow and algorithms (known structure, processed or raw ZULF spectra)

This is the working record of how a ZULF spectrum (or a concentration series) is turned into fitted couplings
with honest uncertainties, as practised on the Blake pyridine series and on isopropylamine (October 2026). It
names every step, the algorithm behind it, the code that does it, the settings that worked, and the pitfalls that
were found. Each analysis keeps its own log under docs/analysis/; this file is the method. Update it when a step
changes (AGENTS.md: keep the record current).

Overview:

    raw FID(s)
      1. diagnostics and processing (per data set: crop, window, SG, zero fill)
      2. phase (switching edge + instrument calibration or data-only symmetry)
      3. structure -> spin-system model (motif, isotopologues, exchange regime, couplings to fit)
      4. forward model through the same processing; complex comparison
      5. objective (weighted residual, missing-peak rows, priors, optional residual-peak rows)
      6. parameterization (one spectrum, monotone or free series, shared parameters, rate families)
      7. search (multi-start, seeds, nested starts, smoothing, coordinate scan)
      8. fine structure (line table, peak sources, residual peaks, targeted local fit, then global)
      9. reliability and uncertainty budget
     10. report (figures with display-only baseline, trace viewer, analysis log, commit)

## 1. Diagnostics and processing (per data set)

Code: zulf_processing (process_dataset, plan_for_dataset, diagnostics.switching_edge, diagnose_raw);
skill skills/zulf-fid-processing (checklist "Tune the parameters for every new FID set").

- Sampling rate from the data or the instrument settings (ini), never by habit.
- Diagnostics: switching edge (half height), ADC plateau, ringing end (about 500 Hz, ends near 50 ms on the
  NMRduino), first-point anomaly.
- Crop start: scan (e.g. 0.03 ... 0.5 s) and look at the main-line FWHM and a valley depth; take the earliest
  start where both stop changing. Pyridine 0.1 s, isopropylamine 0.1 s (lines 0.53 Hz from 50 to 100 ms).
- Record length: the signal extent (lines decay at 1-7 1/s), e.g. 0.1-8.1 s; longer adds noise only, shorter
  costs resolution.
- Window: compare 0, 0.3, 0.6, 1 1/s (noise/max against width); 0.3 1/s is the usual choice.
- SG drift removal on the full record before the crop (process_record); zero fill 3 (0.042 Hz grid at 8 s).
- One processing rule for every spectrum of a series (per-spectrum hand baselines created an outlier band).
- Do not remove what the processing creates: a crop at t_c gives every line oscillating wings (first lobe about
  -0.21 of the height at 0.1 s). Fit with the model rendered through the same processing instead of flattening
  them (AsLS on the data biased 1J(C4,H4) of pyridine by -3.5 Hz).

## 2. Phase

Code: zulf_processing.phase (calibrated_phase), zulf_core.render.phasing; skill skills/zulf-phasing.

- First order: the delay at this FID's own switching edge (plus the instrument offset).
- Zero order: the instrument calibration of the sequence (configs/confirmed_samples.json, "phase_calibration")
  or a data-only estimate by peak symmetry with a model-spectrum bias check.
- In a complex fit the shared-phase gain refines a residual zero-order phase, and a model delay (per spectrum or
  shared) is fitted; with the data phased at the edge the fitted delay should stay near the edge (pyridine:
  -1.9 to -4.5 ms against -3.2 ms; a free per-spectrum delay that jumps signals another problem).

## 3. Structure to spin-system model

Code: zulf_hypothesis (motifs, fragments, builder.build_model, fit.exchange_variants); scripts/regression_confirmed
(structure_for); scripts/fit_processed_spectrum.override_couplings.

- A motif (e.g. "pyridine ring", "(CH3)2CH-NH2") gives the fragment: sites, protons, symmetry, default couplings.
- Natural-abundance isotopologues (13C at each inequivalent site, 15N) become components with fixed or free
  abundance ratios (variant "ratios").
- Exchange regime of N-H / O-H protons: fast (protons dropped: decoupled) or slow (kept); the slow model
  contains the fast one (couplings to the exchangeable protons at 0). Intermediate exchange: Liouville model
  with a fitted rate (D45, fit_staged --exchange).
- Couplings the template leaves unspecified are fixed at 0. List them with their |df/dJ| (line table, section 8)
  and free those that move lines by more than a line width (isopropylamine: methyl-methyl 4J(H,H), fitted
  +0.19 Hz).

## 4. Forward model

Code: zulf_core.physics (Hamiltonian, sectors, transitions), zulf_core.render (Renderer, acquisition operator,
NUFFT), zulf_core.solver.forward.MixtureForward.

- Zero-field Hamiltonian H = sum_k J_k 2 pi I_a . I_b, block diagonal in total-spin sectors; transitions
  (frequency, amplitude) from the eigen-decomposition with the sudden-drop protocol (gamma-weighted preparation
  and detection; all amplitudes positive).
- Each transition is a damped complex exponential exp(2 pi i f (tau + t) - R t) on the acquisition times,
  synthesised by a NUFFT and pushed through the same acquisition operator as the data (SG, crop, mean, window,
  zero fill, phase reference), then evaluated on the observed grid: the model carries the crop wings, the window
  and the record length exactly.
- Linear parameters by variable projection: per component a complex gain (or a shared phase with real
  amplitudes), a per-band background polynomial and nuisance terms are solved by least squares inside every
  evaluation; the optimiser sees only the nonlinear parameters (couplings, decay rates, delay).
- Complex comparison (real and imaginary part) when the processing is linear; real part only for supplied
  processed spectra.
- Speed switches (for local work): `MixtureForward.line_band_hz` renders and differentiates only transitions in a
  band; `MixtureForward.jacobian_only` computes only chosen Jacobian columns.

## 5. Objective

Code: zulf_core.solver.forward (weighting), scripts/fit_joint_series.JointSeries (residual, rows).

- Weighted residual W (model - data) / norm. Signal weighting (`--signal-threshold 2.5 --signal-taper 4`):
  data-driven peak cores full weight, a Gaussian fall-off, a floor elsewhere.
- Missing-peak rows (`--peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 3`): for every data peak top
  (prominence >= 12 % of the maximum and >= k noise sigma), one row lambda * min(max over +-0.15 Hz of r, 0): zero
  while the model reaches the data top, growing when it stays below. Smooth version (log-sum-exp soft max and
  soft hinge) for the optimiser; solutions are rescored and ranked with the hard rows. It keeps small peaks
  without trading the tall ones (peak-height weighting did trade them and was dropped).
- Gaussian priors on couplings (series average around literature values, `--prior-sigma-hh/-ch --prior-weight`).
  Caution: on the pyridine series the series-average priors split the series between two coupling branches and
  produced a jump; check against each spectrum alone (two-branch test).
- Residual-peak rows (`--residual-peaks S`, section 8): a spread-out residual hardly affects the couplings, a
  localized one marks a structural mismatch; such windows are weighted in an outer loop.

## 6. Parameterization

Code: scripts/fit_joint_series.JointSeries, zulf_core.solver.parameterization.

- Series: every coupling monotone in the concentration (`--shape monotone`: J = v + A c(w), c a softmax
  cumulative profile, direction and shape free) or independent at every concentration (`--shape free`, also the
  mode for one spectrum).
- Shared spectrum parameters (`--shared phase_delay`: one model delay for the series).
- Decay-rate families (`--family-edges`): every isotopologue's transitions split by frequency into families, one
  rate each, with bounds (`--rate-bounds 0.2,15` so that no line is switched off by an extreme width). Different
  transitions decay differently (rates depend on the eigenstates of each coherence; residual-field broadening
  scales with each line's g factor); frequency families are an empirical stand-in for a relaxation model.
  Isopropylamine: one rate per isotopologue 0.146, three per isotopologue 0.107, one per line cluster 0.050,
  couplings within 0.16 Hz; the 2J-band lines decay about twice as fast as the J band.
- Per-entry fit ranges in a series file (e.g. to leave out mains harmonics at 120 and 240 Hz).

## 7. Search

Code: scripts/fit_joint_series (main, _solve_start, _coordinate_scan).

- Multi-start: seeds (earlier solutions, literature, other groups' values), perturbed starts (`--spread`), starts
  drawn from the priors; `--workers` processes (one BLAS thread each).
- Nested starts for nested models: start the slow-exchange model from the fast solution with the exchangeable
  couplings at 0 (isopropylamine: generic slow starts ended 26 % above the fast fit, the nested start 13 % below).
- Branches: a series can sit in two near-equivalent coupling branches; start a run inside each branch
  (pyridine: branch A started from the single-spectrum fits gave a 22 % lower total than every earlier run).
- Optional coarse-to-fine smoothing schedule and coordinate grid scans for barriers along single couplings.
- Ranking always on the hard objective (hard missing-peak rows; residual-peak rows on one common window set).

## 8. Fine structure: from a misfit to the couplings that cause it

Code: zulf_core.physics.lines (line_table), scripts/line_table.py, JointSeries.find_residual_peaks,
JointSeries.peak_sources, JointSeries.local_fit, the residual-peak stage in fit_joint_series.

Not every step: as an occasional diagnostic and as a final fine-tuning stage.

1. Line table (analytic). The Hamiltonian is linear and homogeneous of degree 1 in the couplings, so every line
   frequency splits exactly into coupling contributions, f = sum_k J_k df/dJ_k (Euler), with
   df/dJ_k = <a|H_k|a> - <b|H_k|b> (Hellmann-Feynman). `line_table` lists every transition with its amplitude,
   df/dJ for every coupling (zero and tied ones included) and the contributions; `--second-order` adds
   d2f/dJk dJl (cross terms), the nearest other line, a trust step and a near-degenerate flag. Second-order terms
   are for the report only; the fits always recompute the spectrum exactly. Near-degenerate lines (isopropylamine
   122.13 / 122.15 Hz, 0.026 Hz apart) mix strongly: a Taylor step of 1 Hz is off by 0.15 Hz even at second
   order.
2. Residual peaks. |model - data| (complex residual, both signs) with a prominence of k sigma over a robust local
   noise level (1.4826 MAD within +-5 Hz), narrower than 1.5 Hz: ripples on a broad misfit are not peaks. The
   sign of the real residual is reported (model above data: a filled dip or a too-tall line; below: a missing
   line). A peak is assignable when a model transition lies within 0.6 Hz; otherwise it is reported only
   (impurity, instrument line, missing species) and never penalised.
3. Peak sources. For an assignable peak: the transitions under it (component, frequency, amplitude) and per
   coupling its local effect, its selectivity (local / whole-spectrum effect), `shift` (amplitude-weighted mean
   df/dJ: moves the feature) and `split` (their spread: changes the spacing inside the feature, e.g. opens a
   dip). A one-bond coupling shifts a cluster as a whole and cannot open a dip; the split-type couplings can.
4. Targeted local fit. Only the window (1-3 Hz), only the chosen couplings and the rate family of those lines
   free, everything else held; a forward model of the window points alone that renders only transitions within
   the window +- 2 Hz and computes only the needed Jacobian columns (isopropylamine window: Jacobian 0.21 s
   against 5.4 s for the full model). Multi-start inside the window is cheap.
   Caveats (isopropylamine NH2 test): the window forward solves gains, phase and background again on the window,
   so its cost overstates what survives in the full model (local 0.0198, the same couplings 0.122 in the full
   model window); to do: hold gains and phase at the global values. A coupling started at 0 between an
   otherwise uncoupled group and the rest (NH2 protons added to a fast-exchange solution) is a stationary
   point (symmetric splitting, zero gradient): give it nonzero starts.
   Width before position: if narrowing the lines of the feature (one rate family) gives one peak where the data
   have two, the misfit is a line position, not a width.
5. Back to global. Every local candidate becomes a global start; it is accepted only if the global objective
   improves (and, with residual-peak rows, on one common window set). A feature whose local optimum the rest of
   the spectrum rejects points to a model error (missing coupling, species or line-shape physics), to be reported.
5b. Sharp structure rows (`--peak-penalty S --dip-penalty D --peak-max-width W`): rows at every sharp data line
   (model below it: a sharp negative residual peak) and every sharp data valley (model filling it); together they
   stop one broad line from covering a resolved splitting (the fit otherwise widens a rate family over a
   doublet), while a broad data line refined into close model lines costs nothing extra.
6. Residual-peak stage (`--residual-peaks S --residual-peak-rounds 3 --residual-peak-candidates 3`): outer loop
   detect -> weight windows -> refit. Known issues (isopropylamine, to fix): assignability can flip between rounds
   when the lines move away (keep a window once found, judged on the starting transitions), and the ranking chose
   a solution 45 % worse on the plain objective for a 9 % better weighted one (cap the plain loss).
7. Interactive tuning (`scripts/j_tuner.py`, the fit's own command line plus `--fit RUN/fit.json`; open
   http://127.0.0.1:8765). The problem is built by `fit_joint_series.build_problem`, so every number is the
   fitter's objective. Per coupling a coarse slider (+-6 Hz, 1J +-3 Hz around the baseline), a fine slider
   (+-0.3 Hz around the value at the start of the drag) and a number field; on every change (about 0.3 s for
   isopropylamine) the objective with the missing-peak rows, the cost in the view, the plain and relative
   residuals and their change against the baseline, data / model / residual of the view, the model lines and the
   residual peaks (both signs, assignable or not). On demand: the sources of a residual peak (couplings ranked
   by local effect with shift / split, "free only" ticks them for the refinement) and suggested steps (gradient
   and the one-coupling Gauss-Newton step of every coupling; 2 s). "Start refinement": least squares on the full
   objective with the ticked couplings (plus the rate families in view, or every spectrum parameter) free and
   everything else held, streamed to the page, stoppable (the best point is kept), undoable. "Accept" makes the
   state the baseline; "Save" writes OUT/tuned.json, which fit_joint_series reads with --from-joint (a manual
   point as the start of a full multi-start fit). Use it to test a hypothesis by hand (move one coupling, watch
   the window and the global cost together) before writing a scripted local fit.
8. Watching a fit (`scripts/fit_monitor.py runs/processed` -> http://127.0.0.1:8770, or `--text` in a terminal).
   fit_joint_series records every run in OUT/monitor (on by default, `--monitor off`): the objective of every
   residual evaluation of every start with its stage (start / smoothing / fit / coordinate scan / residual-peak
   stage), the best point so far (couplings and the full vector, at most once a second), the start's own starting
   couplings, the phase of the run and the console output. The hook only reads the objective the fit has already
   computed; a test checks that a start gives bit-identical results with and without it. The page shows the
   objective curves of all starts, a table of the starts, the couplings of a chosen start against its start, the
   console, and on request data and model at that start's best point (rebuilt in the viewer process from the
   recorded command; the fit is not touched).

9. Band diagnosis (`scripts/band_diagnosis.py <fit options> --fit RUN/fit.json [--band lo,hi] --figure OUT.png`).
   - Misfit share per isotopologue band.
   - Per bad band, every free parameter as a lever: wanted step, band gain, cost elsewhere, selectivity.
   - Verdict per lever:
     - "free knob": refit it;
     - "conflict": the band and the rest want different values of a shared parameter;
     - "weak": no real gain.
   - When every lever is weak, no free parameter acts on the misfit and the model lacks something there.
10. Component search (`--component-search start|end|both`, default both).
    - Every isotopologue on the window it dominates, from random starts of its own couplings and rate families.
    - Accepted by the whole objective only. Window fits re-solve gains and phase, so they look better than they
      are.
    - The start pass becomes the first centre of the multi-start; the end pass is followed by a global refit,
      kept only if it wins.
    - Record: fit.json["component_search"].
11. When the model lacks something (checklist in skills/zulf-hypothesis-refinement), cheapest first:
    1. **Line list against data peaks** (`TuningSession.lines`):
       - a data peak with no model line is a position problem, which rates cannot fix;
       - a model line where the data show none, sharp because it shares a rate family with sharp lines, needs a
         family edge between them.
    2. **Narrow rate families** around a sharp line that shares a family with broad ones. Remap the old rates onto
       the new edges.
    3. **A second species with a free ratio** (`--structure '[{...}, {...}]'`). Seed it away from species 1: a copy
       has zero gain and no gradient. If seeds stay at the old objective, grid its couplings with the gains solved.
    4. **Model-line passes** (`--model-line-passes N`): the model's own line envelope joins the signal-weighting
       cores and the best solution is refit. Compare runs by
       fit.json["model_line_passes"]["objective_original_weights"].
    Tools for these steps:
    - `--family-edges auto`: edges from the start vector's line clusters (gaps > 1.5 Hz) and from the sharp data
      peaks (each line under one gets its own family). The edges used are printed and recorded in fit.json.
    - `--from-joint` remaps the decay rates when the previous fit used other edges, so changing edges does not
      restart the rates.
    - `--tie-rates REGEX` gives the matching components one decay rate for all their families, so a species
      cannot hide a predicted line by broadening it (for example `'^P2:'` for a second species).
    - `--component-search-hold-gains` (`local_fit(hold_gains=True)`): window costs are the window part of the
      whole objective.
    - `scripts/snapshot_fit.py RUN OUT.json`: the best vector of a running or killed run as a `--from-joint`
      start.
12. Tests of a fitted parameter's sign or of a regime:
    - flip it (for example every 2J(C,H)) and refit everything, then compare objectives;
    - for N-H exchange, check the predicted natural-abundance 15N-H lines against the data before accepting a
      slow-exchange gain.

## 9. Reliability and uncertainty budget

Code: scripts/reliability_series.py; scratchpad budget scripts of the pyridine analysis
(docs/analysis/2026-10-02_blake-pyridine_fid-joint-fit.md).

- Linearised errors of fit.json are 10-30 x too small (correlated residual, model misfit, several minima).
- Reliability classes: all stored solutions of comparable fits pooled on the data-only score (priors removed,
  hard missing-peak rows kept, settings checked equal), near-equivalent set within 3 %; reliable if the spread
  is <= 1.6 Hz, trend only if all agree in direction by > 1 Hz, otherwise not determined.
- Budget per coupling and concentration, combined in quadrature: (1) noise: linearised error x sqrt(residual
  correlation length); (2) near-equivalent solutions: half range of the 3 % set; (3) processing: local refits of
  processing variants (crop start, window, record length; the zero-order phase is absorbed by the gain);
  (4) recovery: synthetic series = model(truth = reference + constant offsets) + the real complex residual,
  refitted from the reference (check that the recovered solution scores at or below the truth). The recovery
  part dominates; it does not cover a model error.

## 10. Report

- Figures of complex fits: real part (and imaginary), with the same data-derived AsLS baseline subtracted from
  data and model for display only (the residual is unchanged).
- Publication figure: `scripts/paper_figure.py` (skills/zulf-figures): whole range and detail bands, experiment and
  simulation overlaid, isotopologue bars, RDKit structure insets, display processing of the whole record and the
  two-step display baseline (anchor spline, then AsLS under the line clusters only); caption file with every
  setting. Couplings on the molecule: `scripts/coupling_diagram.py` (width ~ sqrt|J|, sign colour, reliability
  line style, J +- sigma from the noise term and the spread over model variants).
- J-trend figures annotated with the n-bond type and ring positions, literature values and reliability classes;
  error bars from the budget.
- Fit trace: `fit_joint_series --trace N` writes N frames of the best start's path (and the residual-peak
  stage): OUT/trace.npz, OUT/trace.json; `scripts/trace_view.py OUT` makes OUT/trace.html with a slider and a
  play button (data, model, residual, objective along the path, couplings of the frame).
- Every analysis: its log in docs/analysis/ (data, question, commands, results with numbers, figures, conclusion,
  open points, commit), an index row in docs/analysis/README.md, a line in docs/ANALYSIS_LOG.md, PLAN updated,
  reusable lessons in the matching skill.

## Typical command (one spectrum, complex, rate families, fine structure)

    python scripts/fit_joint_series.py --series series.json --real-only false --shape free --exchange fast \
      --range 85,265 --structure '{"motif": "(CH3)2CH-NH2", "one_bond": {"C1": 133, "C2": 125, "N1": -65}}' \
      --couplings '{"J(HC2,HC3)": 0.0}' --family-edges "50,118.5,121,...,253.5" --rate-bounds 0.2,15 \
      --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 3 \
      --seeds seeds.json --starts 16 --spread 1.0 --workers 4 --trace 60 \
      [--component-search both] [--model-line-passes 1] [--residual-peaks 3] --out runs/processed/NAME
    python scripts/band_diagnosis.py <the fit_joint_series options> --fit runs/processed/NAME/fit.json \
      --bands 2 --local-fit 3 --figure runs/processed/NAME/bands.png
    python scripts/line_table.py --structure '...' --exchange fast --fit runs/processed/NAME/fit.json \
      --band 120,124 --second-order
    python scripts/trace_view.py runs/processed/NAME --band 112,140
    python scripts/make_series_entry.py --fid DATA/average_fid.npy --id NAME --out runs/series/NAME
    python scripts/j_tuner.py <the fit_joint_series options of the fit> --fit runs/processed/NAME/fit.json
    python scripts/fit_monitor.py runs/processed            # live view of running fits (or --text)
    python scripts/plot_runs.py OUT.png runs/processed/A runs/processed/B --zooms "186,196;196,204"   # compare runs
    python scripts/plot_components.py runs/processed/NAME OUT.png --band 176,216   # every component alone
    python scripts/snapshot_fit.py runs/processed/NAME snap.json   # best monitor vector -> --from-joint file
