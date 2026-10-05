# Handoff (2026-10-05)

State of the project, the work in flight and how to continue. Read in this order:
1. `AGENTS.md` (rules);
2. `docs/ARCHITECTURE.md` (interfaces);
3. `docs/PLAN.md` (ledger);
4. `docs/WORKFLOW.md` (how a known-structure fit is done, step by step).

`docs/ANALYSIS_LOG.md` is the running lab notebook, with one log per analysis in `docs/analysis/` (index:
`docs/analysis/README.md`). `docs/DECISIONS.md` (D1-D46) gives the reasons. The previous handoff (2026-09-28)
covered processing; its open items are kept under "Open work" below.

## Repository

- GitHub `gao-xh/ZULF_Model` (the clone's remote `gao-xh/zulf_model` redirects there), branch
  `claude/eager-franklin-tp94gf`. This is the
  only branch on the remote; there is no `main` yet.
- Packages:
  - `zulf_core`: NumPy/SciPy physics (including chemical exchange), rendering, the single processing operator
    (`render/acquisition.py`), solver.
  - `zulf_model`: PyTorch models, generator, training, agent tools, J-to-structure route B
    (`zulf_model.structure`).
  - `zulf_hypothesis`: programmatic hypotheses, motifs, blind search, known-structure fits, reports, J-to-structure
    route A (`j_structure`).
  - `zulf_processing`: experimental processing per dataset (D44).
- Main scripts for experimental fits (all in `scripts/`):

| Script | Purpose |
|---|---|
| `fit_joint_series.py` | complex joint / single-spectrum fitter |
| `make_series_entry.py` | FID -> processed spectrum + series file |
| `band_diagnosis.py` | which band is bad and which parameters act on it |
| `j_tuner.py` | interactive sliders |
| `fit_monitor.py` | live view of running fits |
| `paper_figure.py` | publication figure |
| `line_table.py` | every line with its df/dJ |
| `j_structure.py` | J network -> heavy-atom structure |

- Skills for agents are under `skills/`:
  - blind analysis, FID processing, phasing, spectrum interpretation, model;
  - hypothesis refinement, which holds the fine-structure workflow and the "model lacks something" checklist;
  - figures.

## Setup and checks

```bash
pip install -e ".[dev]"
python scripts/check_ascii.py                      # must pass (code and docs ASCII English)
python -m unittest discover -s tests               # 254 tests; several hours on 4 cores (test_hypothesis alone > 30 min)
python -m unittest tests.test_j_tuner              # fitter tools (tuner, monitor, band diagnosis, component search), 15 s
python -m unittest tests.test_solver.SignalWeightingTests   # signal weighting incl. add_signal_cores, 2.5 min
python -m unittest tests.test_exchange tests.test_render     # exchange physics and renderer
python -m unittest tests.test_processing           # processing only, about 1 min
```

## Data (not in git)

Experimental FIDs are read-only and never committed (`*.npy` and `runs/` are ignored).

**Confirmed samples** (`configs/confirmed_samples.json`, the regression set). Put the files in one folder and set
`ZULF_DATA_DIR`:

| id | file | compound |
|---|---|---|
| e3d282da | e3d282da-average_of_two_batch_fids.npy | L-alanine |
| b683220d | b683220d-average_fid.npy | lactic acid |
| 7ad4aafb | 7ad4aafb-average_fid.npy | pyridine |
| fde3fbb2 | fde3fbb2-average_fid.npy | ethylenediamine |
| 4322bdfc | 4322bdfc-average_fid.npy | triethylamine |
| e66a4b08 | e66a4b08-average_fid.npy | N-ethylmethylamine |

**Other data:**
- Isopropylamine: the 10000-scan average (development molecule, never a blind test).
- Blake pyridine series: processed spectra and raw FIDs (Phase 3d).

All are from the NMRduino: 4 kHz, 65516 points (16.4 s), sequence standard_zf_4000Hz_no_dead.seq. Source
inventory: PLAN "Experimental data inventory".

**Container caveat.** In the cloud session, the uploaded files, `runs/` and the scratchpad are lost when the
container is reclaimed. Everything needed to rebuild a result is in git:
- the commands, in the analysis logs;
- the couplings, in the analysis logs and the tables below;
- the series files, rebuilt from the FIDs with `make_series_entry.py`.

## Current state 1: known-structure complex fits of the amines (2026-10-03 .. 10-05)

Logs:
- docs/analysis/2026-10-03_isopropylamine_complex-fit.md
- docs/analysis/2026-10-04_amines_complex-fit.md
- docs/analysis/2026-10-05_amines_overnight.md

**Common settings:**
- Processing (make_series_entry.py defaults): crop 0.1-8.1 s (start sample 400), exponential 0.3 1/s, zero fill
  3, phase from the instrument calibration.
- Fit options: `--real-only false --shape free --exchange fast --rate-bounds 0.2,15 --signal-threshold 2.5
  --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 3`.
- The 180 Hz instrument line is left out of the ranges.

| Sample | Best run | Objective: stage 1 -> now | Data residual | Model |
|---|---|---|---|---|
| ethylenediamine | eda_2sp_fam198 | 0.1015 -> **0.0051** | 0.072 | two AA'BB' units, free ratio |
| triethylamine | tea_fam121 | 0.3115 -> **0.0383** | 0.196 | N-ethyl (Et3N) |
| N-ethylmethylamine | nema_fam121_final | 0.3129 -> **0.0458** | 0.215 | chain C1-C2-N1-C3 |
| isopropylamine | ipa_fast_fam2_4j | 0.146 -> 0.0475 | 0.217 | (CH3)2CH-NH2, 4J(H,H) free |

**Rate-family edges of the best runs (Hz)** (needed to reproduce):
- eda_2sp_fam198: `185.0,195.3,198.3,199.5,202.4`; range 85-240.
- tea_fam121: `120.8,121.2,121.6,122.6,124.4,126.2,153.3,182.0,185.3,190.6,192.2,193.8,194.8,199.2,200.1,205.5,
  222.8,241.1,247.9,252.7`; range 85-265.
- nema_fam121_final: `121.0,121.65,122.8,126.1,129.3,131.0,133.0,156.8,183.3,189.5,192.6,205.2,211.9,223.4,227.7,
  234.2,241.1,247.8,252.7,257.7`; range 55-280.
- ipa_fast_fam2_4j: `50,118.5,121,123.5,126,128.8,131,132.7,134.2,136,137,138,190,241,245.5,249.8,253.5`;
  range 85-265.

**Best couplings (Hz)** (start values for a rebuild with `--seeds` or a `--from-joint` file):

| Sample | Couplings |
|---|---|
| ethylenediamine, species 1 | 1J 131.34, 2J -1.99, geminal -11.78, J 1.93, J' 14.69 |
| ethylenediamine, species 2 | 1J 126.68, 2J -3.23, geminal -18.01, J 5.85, J' 15.46; gain ratio 0.26 |
| triethylamine | 1J(C1,H1) 130.69, 1J(C2,H2) 125.00, 2J(C1,H2) -4.86, 2J(C2,H1) -3.07, 3J(H1,H2) 7.17, J(C1,HX) -1.72, J(C2,HX) -0.46, J(H1,HX) -0.21, J(H2,HX) 0.38 |
| N-ethylmethylamine | 1J(C1) 124.93, 1J(C2) 131.57, 1J(C3) 131.75, 2J(C1,HC2) -3.37, 2J(C2,HC1) -4.50, 3J(C2,HC3) 5.80, 3J(C3,HC2) 4.27, 3J(HC1,HC2) 7.24, small: J(C1,HC3) -0.36, J(C3,HC1) -0.35, J(HC1,HC3) 0.30, J(HC2,HC3) -0.22 |
| isopropylamine | 1J(C1) 133.50, 1J(C2) 124.43, 2J(C1,HC2) -4.07, 2J(C2,HC1) -1.57, 3J(HC1,HC2) 5.94, 3J(C2,HC3) 5.07, 4J(HC2,HC3) 0.19 |

**Conclusions so far** (all conditional numerical results; see the logs for boundary hits and alternatives):
- **N-H exchange is fast** in all four amines:
  - no natural-abundance 15N-H lines in the data;
  - static N-H models are worse (isopropylamine 0.074 against 0.0475; ethylenediamine 0.0090 against 0.0051);
  - their N-H couplings fit to within 1.3 Hz of zero;
  - the fitted N-ethylmethylamine exchange rate goes fast.
- **2J(C,H) is negative** (relative to 1J > 0). Flipping its sign and refitting is 14x (triethylamine) and 40x
  (ethylenediamine) worse. The literature ethanol value is -4.6 Hz (docs/REFERENCES.md).
- **Ethylenediamine needs a second ethylene unit**, whose identity is open:
  - its lines are 190.54, 192.00, 193.71 and 200.25 Hz, with a strong 177.56 Hz line broadened by its own rate
    family;
  - its geminal (-18 Hz) is outside the usual range, and two of its rate families are at the bound;
  - with one shared rate it still gives 0.0103 against 0.0449 for one species;
  - candidates (carbamate from CO2, a model stand-in for 14N) and decisive experiments are in the overnight log;
  - figure: runs/processed/amine_overview/ethylenediamine_two_species.png.
- **Common open misfit:** a methyl 13C line near 1J that the model does not have:
  - isopropylamine 121.67 Hz (60 % of its cost);
  - N-ethylmethylamine 121.33 Hz and its N-CH3 band;
  - triethylamine 235 / 237.5 Hz small lines.

## Current state 2: Blake pyridine series (Phase 3d)

Logs: docs/analysis/2026-09-29_* .. 2026-10-02_blake-pyridine_fid-joint-fit.md.
- Joint complex fit of the raw FIDs (7 concentrations): 1J(C2,H2), 1J(C3,H3) and 1J(C4,H4) are determined.
- The small couplings are at an information limit.
- An uncertainty budget exists around the branch-A reference.

Open items: PLAN Phase 3d (H-H held at literature, per-spectrum priors, pyridinium, line-shape tests).

## Current state 3: J to structure (Phase 7, docs/J_TO_STRUCTURE.md)

- Route A: rule-based ranking with bond orders and a hybridization-dependent 1J.
- Route B: learned bond-count and hybridization likelihood.
- Synthetic top-1 0.91 / top-3 1.00. All four amine networks (`configs/j_networks/`) rank first.
- Open: a GNN route B, per-coupling sigma from the fits, and the readout on a structure-free network.

## How to fit a known structure now (short form of docs/WORKFLOW.md)

```bash
python scripts/make_series_entry.py --fid DATA/average_fid.npy --id NAME --out runs/series/NAME   # check OUT/phased.png
python scripts/fit_joint_series.py --series runs/series/NAME/series.json <common options> \
    --range LO,HI --structure '{...}' --family-edges "..." --seeds seeds.json --starts 16 --workers 2 \
    --max-nfev 300 --model-line-passes 1 --out runs/processed/NAME      # component search runs by default
python scripts/fit_monitor.py runs/processed --text                    # watch
python scripts/band_diagnosis.py <same options> --fit runs/processed/NAME/fit.json --bands 2 --figure OUT.png
python scripts/paper_figure.py <same options> --fit runs/processed/NAME/fit.json --fid DATA/average_fid.npy \
    --segments "lo,hi;..." --figure runs/processed/NAME/paper_figure.png
```

When a band stays bad, follow the checklist in skills/zulf-hypothesis-refinement ("Model lacks something"):
1. line list against data peaks;
2. narrow rate families around sharp lines;
3. a second species with a free ratio (`--structure` as a JSON list);
4. model-line passes.

Compare runs with model-line passes by `fit.json["model_line_passes"]["objective_original_weights"]`.

## Open work (suggested order)

1. **Amines, missing methyl lines** (PLAN Phase 3e):
   - 14N as a spin-1 nucleus with quadrupolar relaxation for the N-bonded carbons;
   - a common minor species (CO2 carbamate / ammonium) for all four samples.
2. **Ethylenediamine species 2:**
   - a refit with its rates tied (one rate per species), to check whether the -18 Hz geminal and the line
     intensities survive;
   - experiments: high-field 1H/13C of the same sample, inert vs air-exposed sample, added CO2.
3. **Intermediate N-H exchange.** The exact Liouville model is too large for 9 spins without equivalence groups
   (35 GiB for ethylenediamine), so it needs an approximation (for example averaging over the N-H spin states).
4. **Independent J values:** compare with DFT-predicted zero-field multiplets (Andrews et al., arXiv 2604.26071).
5. **Fitter items** (PLAN Phase 3e):
   - local fit with gains held at the global values;
   - sticky residual-peak windows;
   - synthetic validation of close splittings;
   - an automatic placer for rate-family edges.
6. **Processing items still open from the 2026-09-28 handoff** (PLAN Phase 3c):
   - renderer / SG cost for long records, needed before the signal-extent window becomes the default;
   - soft crop window and weighted SG inside the single operator, default off. Design: put the new fields in
     `Acquisition`; give the renderer and `savgol_baseline` the same weighted SG; keep the closed form for the
     flat part and sum the Gaussian edges; test against pointwise multiplication and brute-force weighted
     least squares;
   - per-dataset SG window and instrument-line detection;
   - zero-order phase from the data alone.
7. **Pyridine series:** PLAN Phase 3d open items.

## Working rules learned the hard way

- **Process handling.**
  - Never kill processes by command-line pattern (`pkill -f`, and also `pgrep -f` in a script): it can match and
    kill the calling shell. Use explicit PIDs from `ps -eo pid,cmd`.
- **Cloud container.**
  - The container is reclaimed when the session is idle, and background jobs die with it.
  - For long fits: keep the session active (poll in under-10-minute waits) and snapshot the best monitor vector
    (`OUT/monitor/start_*.jsonl`) into a `--from-joint` file before it is lost.
  - Use at most 4 top-level fits with 2 workers each on 4 cores.
- **Starts and windows.**
  - A second species that starts as a copy of the first has zero gain and no gradient: seed it away from species 1,
    or grid its couplings with the gains solved.
  - Couplings that start at 0 between an uncoupled group and the rest are a stationary point: give them nonzero
    starts.
  - Window fits re-solve gains and phase on the window. Judge every local candidate by the whole objective.
- **Speed.** A processing step that takes more than seconds is a bug.
- **Checks and records.**
  - Every numerical routine gets a test against an independent reference.
  - Record what was done and observed in the analysis log as you go, and update the matching skill when a step
    teaches something reusable.
- **Reporting.** Results are conditional numerical candidates: keep flags, residuals, boundary hits and
  alternatives visible in reports.
