# ZULF Studio: user guide

ZULF Studio is a desktop program for zero-field NMR spectra: live simulation with sliders, fits, figures, and
an interface through which AI models can drive it. Design: D48 and D49 (docs/DECISIONS.md); code in
`zulf_studio/`; environment: docs/ENVIRONMENT.md.

## Start

    conda activate zulf
    python scripts/zulf_studio.py                                       # acetonitrile, zero field
    python scripts/zulf_studio.py --series runs/series/acn/series.json --fit runs/processed/acn_field_fam2
    python scripts/zulf_studio.py --structure '{"motif": "N-ethyl (Et3N)", "one_bond": {"C1": 131, "C2": 125}}'
    python scripts/zulf_studio.py --no-gui --series runs/series/acn/series.json     # API only, no window

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
- **Lines**: every transition (isotopologue, frequency, amplitude, relative); Export writes parameters.json,
  lines.csv, spectrum.csv and figure.png.
- **Fit**: runs `scripts/fit_joint_series.py` in the background from the current couplings and field (starts,
  workers, evaluations, trace frames, field on/off, rate families, rate bounds, extra options). Apply result
  copies couplings, field and decay rate into the sliders. **Fit progress**: drag through the trace frames
  (objective and couplings along the fit path); "Copy this frame's couplings" sets them. The field start
  matters: try several (D47).
- **Figure**: the publication figure of `scripts/paper_figure.py` (PNG at the chosen dpi, PDF, SVG and
  caption.txt), preview, Export to a folder, Open folder. If the sliders are exactly an applied fit, that fit is
  drawn; otherwise the current parameters are drawn and labelled "manual parameters (not a fit)".
- **Log**: everything the session did (sources session, api, ai, fit, figure, terminal); also written to
  `runs/studio/studio.log`.
- **Terminal**: shell commands in the repository directory, with the Studio's Python first on PATH.
- **Python**: an interactive console with `session` (StudioSession), `api` (StudioAPI) and `np`.
- **AI assistant**: a language model operates the session (below).
- **AI API**: the address of the JSON API and the tool list.

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

The AI assistant tab sends a request to a model, which then calls the Studio tools in a loop
(`zulf_studio/assistant.py`, executed in-process on the session). Every call and its result go to the transcript
and to the log (source "ai"). Large arrays are summarised before they reach the model. The model is told that a
fitted result is a conditional numerical result.

| Provider | Library | Credentials | Model |
|---|---|---|---|
| Anthropic (Claude) | `anthropic` | `ANTHROPIC_API_KEY`, `ANTHROPIC_AUTH_TOKEN`, or a profile from `ant auth login` | `claude-opus-5-5` (default); the server-side refusal fallback is on |
| OpenAI (e.g. Codex) | `openai` | `OPENAI_API_KEY` | required: a model your account offers (field "model" or `OPENAI_MODEL`) |

Usage: choose the provider and model, write the request (Ctrl+Enter sends), watch the tool calls; Stop ends the
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

Giving Studio the key (the AI assistant tab has the same guide under "Setup guide"):
- paste it into "API key (this session only)" and press "Use for this session": it stays in Studio's memory until
  Studio closes (not saved, not logged); "Forget" removes it; "Check" re-reads the credentials;
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
