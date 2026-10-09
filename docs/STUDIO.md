# ZULF Studio: user guide

ZULF Studio is a desktop program for zero-field NMR spectra: live simulation with sliders, fits, figures, and
an interface through which AI models can drive it. Design: D48 and D49 (docs/DECISIONS.md); code in
`zulf_studio/`; environment: docs/ENVIRONMENT.md.

## Start

    conda activate zulf
    python scripts/run_studio.py                                       # acetonitrile, zero field
    python scripts/run_studio.py --series runs/series/acn/series.json --fit runs/processed/acn_field_fam2
    python scripts/run_studio.py --structure '{"motif": "N-ethyl (Et3N)", "one_bond": {"C1": 131, "C2": 125}}'
    python scripts/run_studio.py --no-gui --series runs/series/acn/series.json     # API only, no window

Options: `--series` a processed spectrum (`scripts/make_series_entry.py` writes it), `--fit` a fit run to apply
(couplings, field, decay rate, trace), `--api-port` (default 8766; -1 switches the API off; a port in use falls
back to a free one, shown in the AI API tab), `--workspace` (default `runs/studio`).

## Modes (D58)

The header has four task modes (View menu, Cmd+1 to Cmd+4); the model (structure, couplings, field, line width)
is the same in all of them, and the session file records the mode. Studio starts in Fit when it is given data,
else in Simulate (`--mode` chooses).
- **Simulate**: the model alone. Left: structure, couplings, field and line width. Plot: the simulation (the
  weighted sum, or every component and the sum); the loaded spectrum only on request. Right: display, the model
  summary, export of the simulation, "Fit this model to data".
- **Process**: scans or an FID to a spectrum, no model (PLAN 8b). Left: the averaged FID (FID ..., or Scan
  folder ... which averages first) and the recipe: crop start and record (in seconds, or in points: "crop and
  record in" switches; the kept samples are shown), window (exponential, 1/s), drift filter
  (Savitzky-Golay length, 0 = default), phase0 and delay (from the calibration and the switching edge, or by
  hand), zero fill, frequency grid, bands not fitted. Every change recomputes the preview after 150 ms
  (`zulf_processing.series_spectrum`, the same steps a saved series uses). Plot: the FID start with the crop and
  the switching edge, and the spectrum of the recipe. Drawer **Scans**: the per-scan deviation and late noise of
  the run (from scans.json) as a plot and a table with keep boxes; keep z <= ..., keep a range (`0-5989`), all;
  **Average the kept scans** runs `average_scans.py --keep` as a job into `runs/studio/averages/` (the run folder
  stays read only) and processes the new average. Right: **Save spectrum and recipe** (`make_series_entry.py`
  with the recipe; series.json, the spectrum and recipe.json; loaded for Fit when done), other sources, then Fit
  or Blind analysis.
- **Fit**: the data, the model and the residual; right: Fit and Figure.
- **Blind analysis**: the data alone; right: the blind analysis and its report.

## Window

Three columns in the order of the work (2026-10-08 layout; which cards and panels show depends on the mode):
- **Left, the model in workflow order**, each a card that folds away (remembered):
  1. **Data**: the loaded spectrum (points, range, fit ranges, series file); Open, Import FID, Scan folder; display
     phase (deg) and delay (ms), which rotate the shown data only; **Auto phase** "match the simulation" (phase
     and delay that best match the current simulation, delay within +-3 ms; needs the lines roughly in place) or
     "data only" (model-free, delay within +-0.5 ms); "also delay" off keeps the delay; Reset sets both to zero.
  2. **Model**: a summary (name, kind, the isotopologues in their plot colours), the model drawn (click: large view;
     see Molecule below) and **Edit model ...**, which opens the model editor window (owner, 2026-10-09: building
     the model is a set-up step, so it no longer fills the left column). The editor has **From molecule ...**,
     **Attach molecule ...**, JSON, and two sources (switch at the top):
     - **From structure**: a registered motif, or the structure JSON (button JSON shows the editor and Build);
       the natural-abundance 13C isotopologues; N-H/O-H fast (decoupled) or slow (kept).
     - **Spin system** (PLAN 8c): components, each with a name, a weight and its isotopes (`13C 1H 1H 1H`;
       Set resizes the matrix), and the J matrix as an upper-triangle table (the lower triangle mirrors it) or as
       JSON text. A cell holds Hz or a variable name; one name in several cells ties those couplings (the
       Couplings card then has one slider per variable and per numeric coupling, key `component:J(i,j)`). Apply
       builds the model (new variables start at 0 Hz). **From structure** turns the current structure into a spin
       system (each isotopologue a component with its abundance weight; equal couplings share J1, J2, ...; rename
       a cell to split two that are equal by chance); **Load** reads a model JSON or a ZULF_NMR_Suite molecule
       folder (structure.csv); **Save** writes the model JSON. Exact diagonalisation (zulf_core.SpinSystem), no
       MATLAB. Fit works as for a structure (D60): the variables and numeric couplings are the fit parameters;
       **fit weights** (on for typed-in systems) fits every component's amplitude, off holds the weight ratios
       (a converted structure holds its abundance ratios).
  - **The model drawing** (in the Model card; D61), drawn with matplotlib in the plot style, RDKit giving only the
    2D coordinates. A structure model shows its heavy-atom skeleton with
    the protons (CH3, CH2), the site labels and the labelled site of every isotopologue in its plot colour; bond
    orders appear only for a model built from a molecule (the structure itself has none). A spin system shows a
    network: one node per group of equivalent spins, one edge per nonzero coupling with its value (width ~ |J|,
    dashed when negative), a component chooser for several components. The labelled site of each isotopologue sits in
    a pill of its plot colour, written 13C. Protons are drawn as their own atoms ("H atoms"; exchangeable ones on O and
    S faint, they are not in the spin model) or grouped on their atom ("CH3 groups"); the font follows the size of
    the view. A spin system converted from a structure
    keeps the molecule; for any other one **Attach molecule ...** (SMILES or mol file, atoms C1, C2, O1, ... in order)
    adds it for drawing only (model, couplings and an applied fit unchanged); a switch shows molecule or network. **From molecule ...** builds the model
    from a SMILES string or mol-file text (zulf_hypothesis.molecule.structure_from_molecule; API tool
    structure_from_molecule); **Large** opens the drawing in its own window. The plot legend and the field badge
    name the applied fit's run, so the main plot and the Monitor tab can be told apart when they show different
    runs.
  - **Display (phase, delay)** (Fit and Blind analysis; folded by default): the display phase and delay of the
    data and Auto phase, for display only (fits fit their own phase and delay).
  3. **Couplings (Hz)**: a slider and a spin box per coupling; "Add / set" adds a coupling by key `J(a,b)` (for
     example `J(C2,HC1)` makes the nitrile 13C isotopologue of acetonitrile visible). **fine sliders** shows a
     second, fine slider under every coupling, field and phase slider.
  4. **Field and line width**: B transverse and B z in nT, |B|, Zero field, decay rate in 1/s (log scale) and the
     Lorentzian FWHM (= rate / pi).
- **Centre, the plot**: data, simulation and residual, line sticks coloured by isotopologue; one bar with the
  view range, part (real, imaginary, magnitude), lines, fit trace, **lock scale** (freezes the display scale of
  the simulation; automatic: least squares on the magnitudes) and the zoom tools. Below it a drawer (drag to
  resize or close): **Lines**, **Jobs**, **Log**, **AI assistant**.
- **Right, what to run**: **Fit**, **Analysis**, **Figure**.

**baseline** (plot bar, Simulate and Fit): display only, the same baseline correction for data and model as the
publication figures (`zulf_processing.display_baseline`: a spline through anchor points away from the lines,
then AsLS under the line clusters), with the display scale taken on the corrected curves; fits never subtract a
baseline (WORKFLOW W2). The drawer shows the pages of the mode: Lines (Simulate, Fit), Scans (Process), Jobs,
Log and AI assistant.

The static field of the model is written in the lower-left corner of the spectrum plot ("(applied fit)" when the
parameters are exactly an applied fit's).

The simulation is a quick look: complex Lorentzian lines with one decay rate, sign convention of the processed
spectra (numpy FFT). It is not the fit model. In Fit mode with an applied fit whose parameters are unchanged, the
plot draws the **fit model (applied fit)** instead: the fit's own forward model (phase, delay, the decay rate of
every family, amplitudes) at its final parameter vector (`z_final` in fit.json; older fits: the best point of the
monitor record), rebuilt in the background from the run's recorded command, in the data's units (no display
scale) and with the same display phase as the data. Its objective equals the fit's reported score. Any change of
couplings, field or line width goes back to the quick look.

Pages:
- **Lines**: every transition (isotopologue, frequency, relative = |amplitude| / strongest line of all
  isotopologues, amplitude, decay time). **decay time (s)** = 1 / decay rate: with an applied, unchanged fit from
  the fitted decay rate of the line's rate family (found from fit.json family_edges_hz as the forward model does);
  otherwise from the one decay rate of the quick look. Tooltip: the rate and FWHM = rate / pi Hz. The lines.csv
  export has both (decay_per_s, decay_time_s).
- **Fit**: runs `scripts/fit_joint_series.py` in the background from the current couplings and field (starts,
  workers, evaluations, trace frames, coupling precision (0.01 Hz default, D55), field on/off, rate families, rate bounds, extra options). Apply result
  copies couplings, field and decay rate into the sliders. The field start matters: try several (D47).
  **Following a fit** (owner, 2026-10-09): when a fit starts, the main view follows it. Every second the best
  point so far of the best start is applied (couplings, field, decay rates, the exact fit model in the plot),
  rebuilt from the run's record with fit_monitor.point_result, the same way fit.json is written; the sliders are
  read-only and a bar at the top of the left column names the run, start, evaluation and objective. A start or
  point chosen in the Monitor tab (click, slider, Starts table) is applied the same way, so the main plot, the
  sliders and the Monitor always show one point; "Latest best" goes back to following. **Stop following** frees
  the sliders. When the followed fit ends, its final result (fit.json) is applied; a fit that was not followed
  waits for Apply result. (The old trace-frame slider is gone; the Monitor records every evaluation.)
- **Figure**: the publication figure of `scripts/paper_figure.py` (PNG at the chosen dpi, PDF, SVG and
  caption.txt), preview, Export to a folder, Open folder. If the sliders are exactly an applied fit, that fit is
  drawn; otherwise the current parameters are drawn and labelled "manual parameters (not a fit)".
- **Analysis**: blind analysis of an averaged FID (`scripts/analyze_sample.py`, W5): the FID of the loaded series
  (its `source_fid`) or one chosen with Browse; mode blind (hypotheses and search) or known structure (the current
  structure, with a labelling); workers; Run / Stop / Open report. The report (blind.md or structure.md) is shown
  next to the form when it is written. The ranking is a list of conditional candidates.
- **Jobs**: every background job of the session (fits, blind analyses, imports, figures) with state, starts
  finished / total, stage, elapsed time, best objective and output folder; Stop, Live monitor, Open folder, Apply
  fit result, Open report. **Live monitor** (or a double-click on a job) opens the job in the **Monitor** tab.
- **Monitor** (drawer, Fit and Blind analysis modes; View > Monitor, Cmd+Shift+M; shown by itself when a fit
  starts): the live view of a fit run inside Studio (WORKFLOW step 8), read from OUT/monitor every 2 s while the
  tab is visible: objective of every start (thin: each evaluation, thick: best so far), data, model and residual
  at the best point so far of the selected start (rendered in a background thread from the run's recorded
  command; the fit is not touched; "follow" re-renders when the best point improves). **Any point of any start**: click
  near a curve of the objective plot (that start, the nearest evaluation; a ring marks it) or drag the "point"
  slider of the selected start; the spectrum and the Couplings tab then show that evaluation's own parameter
  vector (records since 2026-10-08 store it for every evaluation; older records give the best point up to it,
  and the title says so). "Latest best" goes back to following the run. Side tabs Couplings
  (value and change from the start), Starts and Console. A run that has not written its record yet shows a
  waiting line. While the tab is shown the drawer takes 65 % of the middle column; the previous split returns
  with another tab. **Large window** (button in the tab; Run > Large fit monitor window, Cmd+Shift+L)
  opens the same monitor in its own window of 85 % of the screen, which can be resized or maximized. Run > "Live
  fit monitor in the browser" opens the same view as the web page of
  `scripts/fit_monitor.py`, served by Studio on 127.0.0.1.
- **Log**: everything the session did (sources session, api, ai, fit, figure, terminal); also written to
  `runs/studio/studio.log`.
- **AI assistant**: the chat with a language model that operates the session (below); its status line shows the
  provider, the model and whether a key is set.

**Working indicator**: while a job runs, the status bar (on every tab) shows a spinner, the job, its stage
(component search, multi-start fit, global refit, writing results; processing for an analysis), starts finished
/ total, the elapsed time and the best objective so far, a progress bar (filled by finished starts; moving when
the job has no start count) and Stop; the header shows the same as a badge. When nothing runs it shows the last
job and how it ended.

**Cores and load**: the status bar also shows the cores (Apple silicon: performance + efficiency), the load
average and the analysis workers running on this machine from any program (fit_joint_series, analyze_sample,
...; parsed from `ps`, forked workers counted with their run); red when they exceed the cores (tooltip: each
run). Fits and analyses started from Studio use one BLAS thread per worker; the default worker count is what is
free, and starting more than the free cores asks first. Several fits with many threads each made every run slow
(load 95 on 10 cores, 2026-10-08).

**Files** (File menu; also drag files or folders onto the window):
- **Open** (Cmd+O) by content: a session (`.zulfstudio`), a `series.json`, a fit run (folder or its fit.json; the
  fit's structure is taken over), a folder with series.json; an averaged FID (`.npy`) or an instrument scan folder
  opens the import dialog. **Open recent** keeps the last 12; the dialogs remember their folders.
- **Import** (Cmd+I for an FID, File > Import scan folder, or drop it on the window): first an inspection
  (`StudioSession.inspect_path`, nothing written): scans (numbers), sampling rate and where it was read, points,
  record, sequence, first and last scan time, size, a cached average of another tool (`halp_compiled.npy`),
  switching edge, warnings (scans without .ini, unequal scan sizes, no sampling rate), and a preview of the FID
  start and its spectrum. Then the source (average the scans with `average_scans.py`, optionally leaving out scans
  above a robust z, into `runs/studio/averages/<folder>` or a chosen folder; or the cached average) and the
  processing (`make_series_entry.py`: record, crop, frequency grid, bands left out of the fit ranges, sampling
  rate override, name). The last settings are kept. The run folder is never written to; averages kept for other
  work belong in `~/research/<project>/data/processed/<measurement>/`. The import runs as a job (Jobs tab) and the
  spectrum is loaded when it is done.
- **Save session** (Cmd+S) / **Save session as** (Shift+Cmd+S): structure, couplings, field, line width, view,
  display, and the paths of the loaded series and applied fit (JSON; data and fits stay where they are). Opening
  ignores unknown keys and reports a series or fit that no longer exists. The title shows the session name and a
  dot while there are unsaved changes; closing or opening another file asks first, and closing while jobs run
  asks too (they keep running).
- **Export** (Cmd+E): a dialog for one folder (`<name>_<time>`, place remembered) with any of: the spectrum
  (data, simulation, residual; whole spectrum or the view; `spectrum.csv` with the header `frequency_hz,
  data_real, data_imaginary, data_magnitude, simulation_...` and/or `spectrum.npz` with complex arrays), the
  FID (`fid.csv`: time_s, signal), the parameters (`parameters.json`, `couplings.csv`, `lines.csv`, the session
  file), the applied fit (`applied_fit.json`, `applied_J_table.csv`), the plot (PNG, PDF, SVG), and always
  `information.json` (sources with sha256, time, code commit, settings). File > Export also has the plot alone
  and the publication figure.

**View** menu: Cmd+1 ... Cmd+7 bring a page forward (Fit, Analysis, Figure, Lines, Jobs, Log, AI assistant). **Run** menu: Start fit (Cmd+Return), Blind analysis (Cmd+B),
Stop the running job (Cmd+.), Jobs (Cmd+J), Live fit monitor in the browser; View > Monitor (Cmd+Shift+M). The window size and splitters are remembered.

Window size on any screen: the window fits the screen it is on (restored or moved to another monitor: shrunk to
that screen's free area and moved onto it). Tool rows (the plot bar, the Lines and Jobs buttons, the Monitor
options) wrap onto a second row instead of widening the window, and the left and right columns scroll on short or
narrow screens, so the smallest window is about 1040 x 470 px in every mode.

**Settings** (File > Settings, on macOS ZULF Studio > Settings, Cmd+,): the AI assistant configuration
(provider, model, steps, API key, Keychain, setup guide), the AI API address and tool list, and the appearance
(theme: follow the system, light, dark).

**Tools** menu (separate window): **Terminal** (shell commands in the repository directory, the Studio's Python
first on PATH; Ctrl+Shift+T) and **Python** console with `session` (StudioSession), `api` (StudioAPI) and `np`
(Ctrl+Shift+P).

## JSON API (for scripts and AI agents)

While Studio runs, it serves its tools on `http://127.0.0.1:8766` (only this computer can reach it). A call
changes the same session the window shows; the window follows.

| Request | Purpose |
|---|---|
| `GET /api/tools` | all tools with their parameters; `?format=anthropic` or `?format=openai` gives tool definitions for those APIs |
| `GET /api/state` | structure, couplings, field, line width, data, fit, figure |
| `POST /api/<tool>` | call a tool; the request body is its arguments (JSON) |
| `POST /api/call` | the same as `{"tool": NAME, "args": {...}}` |

Answers are `{"result": ...}`; errors are HTTP 400 with `{"error": ...}`.

`curl` (a command-line program that sends web requests; macOS includes it) works from any terminal on this
computer (iTerm2, the Studio Terminal tab, or Claude Code), for example:

    curl -s http://127.0.0.1:8766/api/state
    curl -s -X POST http://127.0.0.1:8766/api/set_field -d '{"transverse_nt": 37, "z_nt": 45}'
    curl -s -X POST http://127.0.0.1:8766/api/auto_phase -d '{"method": "model"}'
    curl -s -X POST http://127.0.0.1:8766/api/lines -d '{"min_relative": 0.05}'

Main tools: get_state, list_motifs, set_structure, set_couplings, remove_coupling, set_exchange, set_field,
set_linewidth, set_view, set_display, auto_phase, lock_scale, lines, simulate, load_spectrum, fit_command,
start_fit, fit_status, stop_fit, apply_fit, load_trace, trace_frame, figure_command, make_figure, figure_status,
export_figure, export, read_log. `simulate` returns full arrays; prefer `lines`, `get_state` and `fit_status`
when a model reads the results.

Claude Code in this repository can also drive Studio: ask in the chat (for example "use the Studio API to scan
J(C1,HC1) from 136.2 to 136.4 Hz"); it calls the API with curl and the window shows each step.

## AI assistant (Anthropic or OpenAI)

The AI assistant tab sends a request to a model (configured in Settings), which then calls the Studio tools in a loop
(`zulf_studio/assistant.py`, executed in-process on the session). Every call and its result go to the transcript
and to the log (source "ai"). Large arrays are summarised before they reach the model. The model is told that a
fitted result is a conditional numerical result.

| Provider | Library | Credentials | Model |
|---|---|---|---|
| Anthropic (Claude) | `anthropic` | `ANTHROPIC_API_KEY`, `ANTHROPIC_AUTH_TOKEN`, or a profile from `ant auth login` | `claude-opus-5-5` (default); the server-side refusal fallback is on |
| OpenAI (e.g. Codex) | `openai` | `OPENAI_API_KEY` | required: a model your account offers (field "model" or `OPENAI_MODEL`) |

Usage: choose the provider and model in Settings, write the request in the AI assistant tab (Ctrl+Enter sends), watch the tool calls; Stop ends the
loop after the current step; New conversation forgets the history; "max steps" limits the tool calls of one
request. The same assistant works from Python:

    from zulf_studio import StudioSession
    from zulf_studio.assistant import StudioAssistant
    s = StudioSession(); s.load_spectrum(series="runs/series/acn/series.json")
    print(StudioAssistant(s, "anthropic").ask("auto-phase the data and report the rms residual"))

### Getting API access

API use is billed per token, separately from a Claude or ChatGPT subscription. An organisation (lab,
university) may already have an API account; ask its administrator first.

- Anthropic: create an account on platform.claude.com, add a payment method under Billing, create a key under
  API Keys (`sk-ant-...`); or run `ant auth login`. Prices per million tokens (input / output): Claude Opus 5.5
  $4 / $20, Claude Sonnet 5.5 $2 / $10, Claude Haiku 4.5 $1 / $5. A Studio session of a few dozen tool calls
  uses a few hundred thousand tokens, about 1-3 USD with Opus 5.5.
- OpenAI: create an API key in the OpenAI platform and choose a model available to the account.

Giving Studio the key (Settings > AI assistant shows the same guide):
- paste it into the key field, tick "Remember (macOS Keychain)" and press "Use": the key is stored encrypted in
  the macOS Keychain (service "zulf-studio", one item per variable; zulf_studio/credentials.py) and loaded every
  time Studio starts, so it is entered once. Without Remember it stays in memory until Studio closes. Never
  written to a file or the log. "Forget" removes it from Studio and the Keychain (or delete it in the Keychain
  Access app); "Check" re-reads the credentials. A key already set in the environment wins over the Keychain;
- or set it in the shell before starting Studio (`export ANTHROPIC_API_KEY=...`, `export OPENAI_API_KEY=...`) and
  start Studio from that shell;
- or (Anthropic) `ant auth login`, a stored profile the library reads by itself (needs the `ant` command-line
  tool, not installed on the owner's Mac).
An `export` in the Studio Terminal tab does not work: each command there runs in its own shell. Never write a key
into code, configs or the repository (this repository is public).

## Python use without the window

    from zulf_studio import StudioSession
    s = StudioSession({"motif": "methyl", "one_bond": {"C1": 136.0}})
    s.set_field(transverse_nt=37, z_nt=45)
    print(s.lines(min_relative=0.05))
