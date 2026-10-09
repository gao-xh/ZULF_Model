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

## Window

Left column (everything recomputes while a slider moves):
- **Structure**: a registered motif or a structure JSON (`zulf_hypothesis.structure_spec`), then Build. The
  natural-abundance 13C isotopologues are listed. N-H/O-H: fast (decoupled) or slow (kept).
- **Couplings (Hz)**: coarse and fine slider per coupling; "Add / set" adds a coupling by key `J(a,b)` (for
  example `J(C2,HC1)` makes the nitrile 13C isotopologue of acetonitrile visible).
- **Static field and line width**: B transverse and B z in nT (fine +-5 nT), decay rate in 1/s (log scale),
  the field magnitude and the Lorentzian FWHM (= rate / pi).
- **Experimental spectrum**: load a `series.json`; display phase (deg) and delay (ms) rotate the shown data only.
  **Auto phase**: "match the simulation" (phase and delay that best match the current simulation, delay within
  +-3 ms; needs the lines roughly in place) or "data only" (model-free, delay within +-0.5 ms); "also delay"
  off keeps the delay; Reset sets both to zero.

Centre: data, simulation and residual; line sticks coloured by isotopologue. Part: real, imaginary or
magnitude. **lock scale** freezes the display scale of the simulation (automatic: least squares on the
magnitudes). The status bar shows the simulation time, the lines in view, the scale and the rms residual.

The simulation is a quick look: complex Lorentzian lines with one decay rate, sign convention of the processed
spectra (numpy FFT). It is not the fit model.

Bottom tabs:
- **Lines**: every transition (isotopologue, frequency, amplitude, relative).
- **Fit**: runs `scripts/fit_joint_series.py` in the background from the current couplings and field (starts,
  workers, evaluations, trace frames, coupling precision (0.01 Hz default, D55), field on/off, rate families, rate bounds, extra options). Apply result
  copies couplings, field and decay rate into the sliders. **Fit progress**: drag through the trace frames
  (objective and couplings along the fit path); "Copy this frame's couplings" sets them. The field start
  matters: try several (D47).
- **Figure**: the publication figure of `scripts/paper_figure.py` (PNG at the chosen dpi, PDF, SVG and
  caption.txt), preview, Export to a folder, Open folder. If the sliders are exactly an applied fit, that fit is
  drawn; otherwise the current parameters are drawn and labelled "manual parameters (not a fit)".
- **Analysis**: blind analysis of an averaged FID (`scripts/analyze_sample.py`, W5): the FID of the loaded series
  (its `source_fid`) or one chosen with Browse; mode blind (hypotheses and search) or known structure (the current
  structure, with a labelling); workers; Run / Stop / Open report. The report (blind.md or structure.md) is shown
  next to the form when it is written. The ranking is a list of conditional candidates.
- **Jobs**: every background job of the session (fits, blind analyses, imports, figures) with state, starts
  finished / total, stage, elapsed time, best objective and output folder; Stop, Open folder, Apply fit result,
  Open report.
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
- **Open** (Cmd+O) by content: a session (`.zulfstudio`), a `series.json`, a fit run (folder or its fit.json), a
  folder with series.json, an averaged FID (`.npy`: imported) or an instrument scan folder (averaged, then
  imported). **Open recent** keeps the last 12; the dialogs remember their folders.
- **Import averaged FID** (Cmd+I): `make_series_entry.py` into `runs/studio/series/<name>` (sampling rate from
  scans.json or the .ini next to the FID; whole-grid fit ranges), loaded when done. **Import scan folder**:
  `average_scans.py` into `runs/studio/averages/<folder>` (averages kept for other work belong in
  `~/research/<project>/data/processed/<measurement>/`), then imported.
- **Save session** (Cmd+S) / **Save session as** (Shift+Cmd+S): structure, couplings, field, line width, view,
  display, and the paths of the loaded series and applied fit (JSON; data and fits stay where they are). Opening
  ignores unknown keys and reports a series or fit that no longer exists. The title shows the session name and a
  dot while there are unsaved changes; closing or opening another file asks first, and closing while jobs run
  asks too (they keep running).
- **Export**: parameters, couplings, lines and spectrum (`parameters.json`, `couplings.csv`, `lines.csv`,
  `spectrum.csv`, a session file and the applied fit's `fit.json` / `J_table.csv`; Cmd+E); the plot as PNG, PDF
  or SVG; the publication figure.

**View** menu: Cmd+1 ... switch the bottom tabs. **Run** menu: Start fit (Cmd+Return), Blind analysis (Cmd+B),
Stop the running job (Cmd+.), Jobs (Cmd+J). The window size and splitters are remembered.

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
