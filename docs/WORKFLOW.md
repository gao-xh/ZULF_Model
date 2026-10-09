# Workflows

How ZULF data and models are used, from instrument scans to couplings and structures. The map below shows every
workflow; each section after it gives the steps, the command and the documents behind it. The fit of a known
structure (W3, or W4 in Studio) is described in full detail in the second part of this file, sections 1-10.

**Keep this file current** (AGENTS.md): a change that adds, removes or reorders a step a person or agent runs
(a script, a route, an option default, an output file) updates the matching diagram and command here in the same
commit. Diagrams are Mermaid (rendered by GitHub, plain text in git); node labels name the script or module.
An HTML page of the map and W1-W7 is generated from this file: `python scripts/workflow_page.py`
(writes runs/workflow/index.html; regenerate it after editing this file).

Last reviewed: 2026-10-09 (D56; checked against the code; key algorithms per workflow; sampling rate from the data; whole-grid fit ranges; guard rows D57; per-peak rate families D59; spin-system fits D60; Studio jobs, files and machine load; Studio modes D58; Studio live monitor; molecules D61; checked 2026-10-09).

## Map

```mermaid
flowchart TD
    scans["Instrument scans<br/>data/original (read-only)"]
    avg["W1 average and screen scans<br/>average_scans.py"]
    fid["Averaged FID<br/>data/raw or data/processed"]
    proc["W2 processing: drift, crop, window, FFT<br/>phase: edge delay + calibration<br/>baseline: modelled, not subtracted<br/>make_series_entry.py"]
    spec["Complex spectrum + series.json<br/>runs/series/NAME"]
    known{"Structure known?"}
    blind["W5 blind analysis<br/>analyze_sample.py"]
    jnet["Fitted couplings J, K<br/>fit.json, J_table.csv"]
    j2s["W6 J network to structure<br/>j_structure.py: A, B-J, B-K"]
    report["Figures (display-only baseline),<br/>analysis log, report<br/>paper_figure.py, docs/analysis"]
    nn["W7 learned candidate model<br/>zulf-model train / propose<br/>(built, not yet trained)"]

    subgraph known_fit["Known-structure fit: one engine, two ways to drive it"]
        fit["W3 command line<br/>fit_joint_series.py"]
        studio["W4 ZULF Studio<br/>run_studio.py: sliders, live simulation,<br/>fit and analysis jobs, imports, figures, AI API"]
        studio -- starts fit jobs --> fit
        fit -. loads results, trace .-> studio
    end

    scans --> avg --> fid --> proc --> spec --> known
    known -- yes --> known_fit
    known -- no --> blind
    blind -- ranked candidate structures --> known_fit
    studio -. starts blind jobs, imports FIDs .-> blind
    blind -. structure-free J network, planned .-> j2s
    fit --> jnet
    jnet --> j2s
    j2s -. candidate structures .-> known_fit
    jnet --> report
    studio --> report
    nn -. candidate hints .-> blind
```

Solid arrows: the path in use. Dashed arrows: optional or planned links. Couplings worth reporting come only from a
fit (W3 / W4): the blind analysis ranks hypotheses with coarse settings and hands its best candidates to that fit. W3 and W4 are the same fit (Studio runs
fit_joint_series.py and reads its fit.json back); Studio adds manual tuning, live simulation and figures.

| Workflow | Input | Output | Main command | Details |
|---|---|---|---|---|
| W1 scans to FID | instrument run folder | average FID, half averages, scan metrics | `scripts/average_scans.py` | below; skills/zulf-fid-processing |
| W2 FID to spectrum | average FID | processed, phased complex spectrum, series.json, fit ranges | `scripts/make_series_entry.py` | sections 1-2, 4, 10; skills/zulf-fid-processing, zulf-phasing |
| W3 known-structure fit | spectrum + structure | couplings with uncertainties, fit.json | `scripts/fit_joint_series.py` | sections 3-10 |
| W4 the same, interactive | FID or spectrum, optional fit | tuned or fitted parameters, fits, blind analyses, figures, sessions | `scripts/run_studio.py` | docs/STUDIO.md |
| W5 blind analysis | FID, no structure | ranked hypotheses and labellings | `scripts/analyze_sample.py` | skills/zulf-blind-analysis |
| W6 J network to structure | fitted couplings | ranked heavy-atom graphs (A, B-J, B-K) | `scripts/j_structure.py` | docs/J_TO_STRUCTURE.md |
| W7 learned candidates | generator configs | trained spectrum -> candidate model | `zulf-model train` | docs/PLAN.md Phases 1, 2, 4 |

## W1. Instrument scans to an averaged FID

```mermaid
flowchart LR
    run["data/original/MEAS<br/>n.dat, n.ini"] --> avg["average_scans.py<br/>decode, screen (robust z)"]
    avg --> out["average_fid.npy<br/>average_even / odd.npy<br/>scans.json, scans.png"]
```

    python scripts/average_scans.py ~/research/zulf/data/original/MEAS ~/research/zulf/data/processed/MEAS \
        [--exclude-z 3] [--keep FILE|0-5989,6005]

- The run folder is read-only; averages made by our code go to `data/processed/MEAS/`. Received averages (other
  tools, other groups) are in `data/raw/MEAS/`.
- The even / odd half averages are disjoint scans for held-out checks of a fit.

Key algorithms:
- Scan screening: per scan, the RMS deviation from the mean of all scans (each detrended) and the late-record
  noise; robust z = (value - median) / (1.4826 MAD); scans above `--exclude-z` are left out and the mean is
  recomputed once without them.
- Half averages by scan position (even / odd), so the two halves share drift but not noise.

## W2. Averaged FID to a processed spectrum: processing, phase, baseline

```mermaid
flowchart TD
    fid["average_fid.npy + ini"] --> diag["diagnostics<br/>switching edge, ringing end,<br/>signal extent"]
    subgraph time["time domain (zulf_processing, one rule per series)"]
        drift["drift removal<br/>SG filter on the full record"]
        crop["crop after the ringing<br/>record length"]
        win["mean removal,<br/>exponential window, zero fill"]
        drift --> crop --> win
    end
    diag --> drift
    win --> fft["FFT: complex spectrum"]
    subgraph phase["phase"]
        p1["first order: delay at<br/>this FID's switching edge"]
        p0["zero order: instrument calibration<br/>or data-only symmetry estimate"]
    end
    fft --> p1 --> p0
    p0 --> ranges["fit ranges<br/>signal bands, mains harmonics out"]
    ranges --> entry["runs/series/NAME<br/>series.json, phased.png"]
    subgraph baseline["baseline: never subtracted before a fit"]
        bfit["in the fit: model through the same<br/>processing + per-band background polynomial"]
        bdisp["for figures only: anchor spline,<br/>then AsLS under the line clusters"]
    end
    entry --> bfit
    entry -.-> bdisp
    bfit -. "residual phase and delay refined in the fit" .-> p0
```

Interactive: ZULF Studio, mode Process (scan selection, recipe with a live preview, save; docs/STUDIO.md).

    python scripts/make_series_entry.py --fid DATA/average_fid.npy --id NAME --out runs/series/NAME \
        [--crop 0.1] [--record 8] [--apodization 0.3] [--zero-fill 3] [--grid 20,380] [--sampling-rate HZ]
        [--exclude lo,hi] [--sg-window S] [--phase0-deg D] [--delay-ms M] [--crop-points N] [--record-points N]

Processing (time domain; details and settings in section 1 below):
- Diagnostics first: the switching edge (half height), the ringing end (near 50 ms on the NMRduino), the signal
  extent. Sampling rate from scans.json or the .ini next to the FID (`--sampling-rate` overrides; sequences run
  at 2000, 4000 or 8333 Hz, and a wrong rate scales every frequency); sample-count defaults are rescaled to it.
- Drift removal: a Savitzky-Golay low-pass on the full record, subtracted before the crop (the slow baseline of
  the FID; a time-domain baseline correction).
- Crop after the ringing (usually 0.1 s) for a record of about the signal extent (e.g. 8 s); mean removal;
  exponential window (0.3 1/s usual); zero fill 3; FFT.
- One processing rule for every spectrum of a series. `zulf_core.render.acquisition.process_record` is the only
  implementation of these steps; data and model both pass through it (AGENTS.md).

Phase (section 2 below; skills/zulf-phasing):
- First order: the delay at this FID's own switching edge plus the instrument offset.
- Zero order: the instrument calibration of the sequence (configs/confirmed_samples.json) or a data-only estimate
  by peak symmetry, checked against a model spectrum.
- The fit then refines a residual zero-order phase (shared-phase gain) and the delay; a fitted delay far from the
  edge signals another problem.
- Studio (W4): a model can come from a motif, structure JSON, a typed-in spin system (D60) or a molecule
  (From molecule: SMILES or mol file -> chain specification, D61; API tool structure_from_molecule).
- Studio (W4): auto phase against the current simulation ("model") or from the data alone ("data"), for display.

Baseline (frequency domain):
- Not subtracted from the data before a fit. The crop gives every line oscillating wings and a rolling baseline;
  the model is rendered through the same processing, so it carries them, and a per-band background polynomial
  is solved inside every fit evaluation (section 4). AsLS on the data before fitting biased 1J(C4,H4) of
  pyridine by -3.5 Hz.
- For figures only: `zulf_processing.anchor_spline_baseline` (spline through anchor points away from the lines)
  and then AsLS under the line clusters (`asls_baseline`), the same steps for data and model; the residual is
  unchanged (section 10, scripts/paper_figure.py).

Key algorithms:
- Switching edge: half height between the start plateau and the first extremum (2-20 ms), interpolated between
  samples (NMRduino: edge at 3.41-3.51 ms). Minus the edge time is the first-order phase delay.
- Drift removal: Savitzky-Golay smoothing (default window 201 samples, order 2) on the full record with mirrored
  edges, subtracted (`process_record`).
- Spectrum: finite-record Fourier sum on the requested grid (`evaluate_spectrum`, FFT when the grid allows it).
- Phase convention: multiply by exp(-i (phase0 + 2 pi f (delay + crop reference))) (`phase_correct`).
- Data-only phase (`estimate_phase`): at the strongest peaks the window-summed complex value has phase
  phase0 + 2 pi f (reference + delay) modulo pi (lines may be negative); delay on a grid, phase0 from the
  weighted circular mean of the doubled angles.
- AsLS baseline (Eilers and Boelens 2005, display only): minimise sum w (y - z)^2 + lam sum (second difference
  of z)^2 with asymmetric weights, lam = (smooth_hz / step)^4 so the setting is grid-independent.
- Fit ranges: the whole grid (default 20-380 Hz, `--grid`) minus +-0.4 Hz around the mains harmonics
  (n x 60.06 Hz) and instrument lines and minus `--exclude` intervals; a model line where the data are empty is
  then rejected by the fit. `--ranges bands` keeps only the bands above 5 noise levels (noise from 330-380 Hz),
  widened by 3 Hz.

Checks: look at `phased.png` before fitting; `zulf-model diagnose AVERAGE.npy 0.ini` and
`scripts/validate_processing.py` compare processing and phase choices.

## W3. Fit of a known structure

```mermaid
flowchart LR
    s["series.json<br/>+ structure motif"] --> m["spin-system model<br/>isotopologues, exchange<br/>through the same processing"]
    m --> f["fit_joint_series.py<br/>multi-start, rate families, field, gamma,<br/>precision; gains, residual phase, delay,<br/>background polynomial per band"]
    f --> d["band_diagnosis.py<br/>line_table.py, j_tuner.py"]
    d -- misfit --> f
    f --> u["reliability budget<br/>noise, near-equivalent solutions,<br/>processing variants, recovery"]
    u --> r["fit.json, J_table.csv (J and K)<br/>paper_figure.py, analysis log"]
```

    python scripts/fit_joint_series.py --series runs/series/NAME/series.json --structure '{"motif": ...}' \
        [--fit-field] [--field t,z --fit-gamma 13C] [--precision 0.01] --out runs/processed/NAME

- Steps, algorithms, settings and pitfalls: sections 3-10 below; a full command at the end of this file.
- Model: `--structure` takes a motif, a chain, or a typed-in spin system `{"spin_system": {"components":
  [{"name", "isotopes", "J" (Hz or variable names), "weight"}], "variables": {...}, "fixed_weights": false}}`
  (D60; variables are the fit parameters).
- Field: zero field by default; `--fit-field` fits the transverse and z components (D47); with a known field
  `--fit-gamma` names a heteronucleus (D53).
- Precision: fits stop when every coupling moves by less than `--precision` (default 0.01 Hz, D55).
- Decay-rate families: `--family-edges` edges in Hz, `auto` (model line clusters, D46), `peaks` (one family per
  data peak, D59) or `peaks+auto`.
- Output couplings carry K next to J (D51).

Model: spin Hamiltonian, signal and rendered spectrum (docs/CONVENTIONS.md; zulf_core.physics, section 4 below):

```math
H = \sum_{i<j} J_{ij}\, \mathbf{I}_i \cdot \mathbf{I}_j \;-\; \sum_{n} \gamma_n\, \mathbf{B} \cdot \mathbf{I}_n
\qquad (\text{in Hz}; \ \mathbf{B} = 0 \text{ by default})
```

```math
\rho(0) \propto \sum_i \gamma_i I_{z,i}, \qquad D = \sum_i \gamma_i I_{z,i}, \qquad
\rho(t) = e^{-2\pi i H t}\, \rho(0)\, e^{2\pi i H t}, \qquad s(t) = \operatorname{Tr}\!\big[\rho(t)\, D\big]
```

```math
s(t) = \sum_k \operatorname{Re}\!\big( A_k\, e^{2\pi i f_k t} \big), \qquad f_k = E_m - E_n > 0, \qquad
A_k \propto \langle n | \rho(0) | m \rangle \langle m | D | n \rangle
```

```math
M_s(\theta, c) = \mathcal{P}\Big[ \sum_{q} c_q \sum_{k \in q} A_k\, e^{2\pi i f_k (t + \delta) - R_k t} \Big]
  + \sum_{b} \mathbb{1}_b(f)\, \phi(f) \sum_{j=0}^{n_{\mathrm{bg}}} c_{bj}\, P_j\big(u_b(f)\big)
```

- $H$: spin Hamiltonian in frequency units (Hz), evolution $e^{-2\pi i H t}$; block diagonal in collective-spin and
  total-$M$ sectors, so it is diagonalised sector by sector. $J_{ij}$: scalar coupling between spins $i$ and $j$
  (Hz). $\mathbf{I}_i$: spin operator of spin (or collective spin) $i$. $\gamma_n$: gyromagnetic ratio
  $\gamma/2\pi$ of nucleus $n$ (Hz/uT; registry value, or fitted, D53). $\mathbf{B}$: static field (uT); only its
  transverse magnitude and $|B_z|$ are observable (D47).
- $\rho(0)$, $D$: sudden-drop protocol, preparation and detection both weighted by $\gamma$ along $z$.
  $s(t)$: detected signal. $E_m$, $E_n$: eigenvalues of $H$ (Hz); $f_k$: line frequency; $A_k$: complex line
  amplitude (per molecule: traces divided by the Hilbert-space dimension).
- $M_s$: model spectrum of spectrum $s$. $q$: component (isotopologue, molecule) with complex gain $c_q$;
  $\delta$: delay of the signal start; $R_k$: decay rate of line $k$ (its rate family, plus an exchange rate
  where modelled). $\mathcal{P}$: the processing operator of W2 (drift removal, crop, mean removal, window, zero
  fill, Fourier sum on the fit grid, phase reference), the same for data and model. $b$: fit band;
  $\mathbb{1}_b$: 1 inside band $b$; $c_{bj}$: complex background coefficients of order $j \le n_{\mathrm{bg}}$
  (1 by default in the known-structure fits); $P_j$: Legendre polynomial; $u_b(f) \in [-1, 1]$: frequency scaled
  across band $b$; $\phi(f)$: the phase correction of the data, so the background is smooth in the record's own
  frame (where baseline leftovers are smooth).
- Reduced coupling reported next to $J$ (D51): $K_{AB} = J_{AB} / (h\, g_A g_B)$, $g = \gamma/2\pi$ in Hz/T.

Loss function (what least_squares minimises; JointSeries.residual, MixtureForward.predict; section 5 below):

```math
L(\theta) = \sum_{s} \Big[ \lVert r_s \rVert^2 + \sum_{p} \rho_{sp}^2 + \sum_{b} \gamma_{sb}^2
  + a^2 \lVert r_s \rVert_{\mathcal{R}_s}^2 \Big]
  + \sum_{k} \Big( \frac{\bar{J}_k - \mu_k}{\sigma_k} \cdot \frac{\sqrt{w_{\mathrm{prior}}}}{\bar{n}} \Big)^2
```

```math
r_s = \frac{W_s \big( M_s(\theta, \hat{c}_s) - Y_s \big)}{n_s}, \qquad n_s = \lVert W_s Y_s \rVert, \qquad
\hat{c}_s = \arg\min_{c} \big\lVert W_s \big( M_s(\theta, c) - Y_s \big) \big\rVert
```

```math
W_{s,ii} = \frac{w_i}{\sigma_{\mathrm{noise}}}, \qquad
w_i = \Big[ 0.2 + 0.8 \exp\!\Big( -\frac{d_i^2}{2 \tau^2} \Big) \Big] \Big( \frac{h_i}{h_{\max}} \Big)^{-q}
```

```math
\rho_{sp} = \lambda \min\Big( \max_{|f - f_p| \le \epsilon} r_s(f),\; 0 \Big)
\quad \text{(hard, for ranking)}
```

```math
m = t \log \sum_{|f - f_p| \le \epsilon} e^{r_s(f)/t}, \qquad
\rho_{sp} = -\lambda\, t\, \operatorname{softplus}(-m/t), \qquad \operatorname{softplus}(x) = \log(1 + e^{x}),
\qquad t = \kappa\, h_p \quad \text{(smooth, while optimising)}
```

```math
\gamma_{sb} = \frac{\lambda_g\, t_g}{\sigma\, n_s}
  \operatorname{softplus}\!\Big( \frac{\max_{f \in b} \big[ (1 - \eta) E_s(f) - D_s(f) \big] - k_g \sigma}{t_g} \Big),
\qquad t_g = 0.5\, \sigma
\quad \text{(guard, D57)}
```

Every symbol:
- $L$: the objective; the value reported as "objective" is $L$ at the solution with the hard rows.
- $\theta$: the nonlinear parameters the optimiser moves: couplings $J_{ij}$, decay rates $R$, delay $\delta$,
  field $\mathbf{B}$ and gamma when fitted. $c$: the linear parameters (gains $c_q$ or a shared phase with real
  amplitudes, background $c_{bj}$, nuisance terms); $\hat{c}_s$: their exact least-squares solution at every
  evaluation (variable projection), so they never enter $\theta$.
- $s$: spectrum of the series (one for a single spectrum); $i$: frequency point of the fit ranges, at
  frequency $f_i$; $p$: data peak top; $k$: coupling that has a prior.
- $Y_s$: the processed, phased complex data spectrum on the fit ranges. $M_s$: the model spectrum above.
- $r_s$: the weighted, normalised residual vector, real and imaginary parts stacked (real part only for spectra
  processed elsewhere). $\lVert \cdot \rVert$: Euclidean norm over those entries. $n_s$: norm of the weighted data,
  so a perfect fit gives $L = 0$ and the empty model $L = 1$ per spectrum.
- $W_s$: diagonal weight matrix. $\sigma_{\mathrm{noise}}$: robust noise level (1.4826 x median absolute
  deviation of the narrow-feature excess). $w_i$: relative weight in $[0.2, 1]$. $d_i$: distance (Hz) from $f_i$
  to the nearest data peak core in its band (cores: narrow features above `--signal-threshold`, default 4
  $\sigma_{\mathrm{noise}}$). $\tau$: `--signal-taper` (default 2 Hz). $h_i / h_{\max}$: height of the nearest
  peak over the band maximum (at least 0.1); $q$: `--signal-height-power` (default 0, factor 1).
- $\rho_{sp}$: missing-peak row of data peak top $p$; zero while the model reaches the data top, negative when it
  stays below. $f_p$: position of the peak top; $h_p$: its height; $\epsilon$: `--peak-tolerance` (default
  0.15 Hz); $\lambda$: `--peak-penalty` (0 = off); for negative data tops (dips) the residual is mirrored and
  $\lambda$ = `--dip-penalty`. $m$: soft maximum of $r_s$ in the window; $t$: its smoothing scale; $\kappa$:
  `--peak-smooth` (0 = hard rows throughout).
- $\gamma_{sb}$: guard row of bin $b$ (D57): the model may not put a peak where the data have no line. Bins of
  `--guard-bin` (default 1 Hz) over the guard points: the whole processed grid (or `--guard-range`) minus
  `--guard-margin` (default 2 Hz) around the data line regions (4 $\sigma$ on the whole grid) and the fully
  weighted fitted points; with fit ranges of a few bands most of the guard lies outside them, with the whole
  grid it is the empty part. $E_s(f)$: the model line envelope (sum of line magnitudes with the solved gains,
  phases ignored); $D_s(f)$: the largest data magnitude within $\pm 0.5$ Hz (line tails and noise allowed);
  $\eta = 0.2$: the share of $E$ not counted (it exceeds the coherent model where lines partly cancel);
  $k_g$: `--guard-sigma` (default 3); $\sigma$: noise of the whole grid; $\lambda_g$: `--guard-penalty` (default
  5, 0 = off). Zero while the model stays under the data there; the data in the guard never enter otherwise.
- $a$: residual-peak strength and $\mathcal{R}_s$: the windows around assignable residual peaks
  (`--residual-peaks`, 0 = off); $\lVert r_s \rVert_{\mathcal{R}_s}$: the norm over those points only.
- Prior term: a Gaussian prior (maximum a posteriori) that penalises a coupling for leaving a fixed reference value,
  not for the size of an optimiser step. $\bar{J}_k$: mean of coupling $k$ over the series (the value itself for
  one spectrum); $\mu_k$: its reference, the structure's starting value (motif default, literature, or
  `--structure` / `--couplings`), fixed during the fit; $\sigma_k$: allowed deviation, `--prior-sigma-hh` for
  H-H and `--prior-sigma-ch` for C-H couplings (0 = no prior, the default); couplings with $|\mu_k| \ge 50$ Hz
  (1J) never get one. $w_{\mathrm{prior}}$: `--prior-weight` (default 1); $\bar{n}$: mean of $n_s$ over the series,
  so the term is on the scale of the normalised data term. A deviation of $\sigma_k$ costs as much as one data
  point off by one noise level. Check a prior fit against single-spectrum fits without it (section 5).
- Optional coarse-to-fine smoothing replaces $r_s$ by $S r_s$ (the same Gaussian kernel on data and model).
- The blind search (W5) scores hypotheses on a common yardstick instead: $\chi^2$ on the data cores and
  $\mathrm{BIC} = \chi^2 + n_{\mathrm{par}} \ln N$ ($n_{\mathrm{par}}$ free parameters, $N$ data points).

Key algorithms (sections 3-9 below give the detail):
- Spin physics: the Hamiltonian and protocol above, diagonalised per collective-spin and total-M sector; transition
  lists do not depend on line widths and are cached (TransitionCache) while rates and gains change. Intermediate exchange uses a Liouville model (D45).
- Rendering: each transition is a damped complex exponential synthesised on the acquisition times by a type-1
  NUFFT (Gaussian gridding, relative error about 1e-12) and pushed through `process_record`, so the model carries
  the crop wings, window and record length of the data.
- Variable projection: complex gains (or a shared phase with real amplitudes), per-band background polynomials
  and nuisance terms are solved by linear least squares inside every evaluation; the optimiser sees couplings,
  decay rates, delay, field and gamma only.
- Optimiser: scipy `least_squares` (trust-region reflective, bounds, analytic Jacobian, x_scale "jac"), many
  starts in parallel processes; the precision callback stops when couplings move by less than `--precision`
  and the cost no longer falls (D55).
- Objective: signal-weighted complex residual plus missing-peak rows (soft max / soft hinge while optimising,
  hard rows for ranking) and optional coupling priors and residual-peak rows.
- Decay-rate families: one rate per frequency family of each isotopologue's lines (`--family-edges`).
- Fine structure: line table with dJ derivatives, peak sources, residual peaks, local fits of the couplings that
  move the lines of one band, then a global refit (section 8).
- Reliability: linearised errors are 10-30 x too small; the budget adds the spread of near-equivalent solutions
  (3 % set), processing variants and a recovery test on synthetic data with the real residual (section 9).

## W4. The same fit, interactive: ZULF Studio

```mermaid
flowchart LR
    u["person"] --> app["Studio window<br/>sliders, fit tab, figures"]
    ai["AI agent / assistant"] --> api["JSON API 127.0.0.1:8766"]
    app --> sess["StudioSession"]
    api --> sess
    sess --> fit["fit_joint_series.py jobs"]
    sess --> fig["paper_figure.py figures"]
    sess --> bl["analyze_sample.py jobs (W5)"]
    sess --> imp["imports: average_scans.py,<br/>make_series_entry.py"]
    sess --> files[".zulfstudio session files,<br/>exports (JSON, CSV, PNG/PDF/SVG)"]
```

    python scripts/run_studio.py --series runs/series/NAME/series.json [--fit runs/processed/NAME]

- Studio is a front end of W3, not a later step: its fit jobs are fit_joint_series.py runs (same options, same
  fit.json), and a command-line fit can be loaded into it (`--fit`); while a fit runs the main view follows its
  best point and any point chosen in the Monitor tab. Use it to look, tune by hand,
  start and watch fits, phase and make figures; record results as for W3.
- Real-time simulation (structure, couplings, field, line width), auto phase, fits with a progress slider,
  figure export, log and terminal; AI keys in the macOS Keychain or environment only (docs/STUDIO.md).
- Jobs (fits, blind analyses, imports, figures) run as subprocesses with one BLAS thread per worker; the status
  bar shows the running job (stage, starts, elapsed, best objective, Stop) and the cores, load and analysis
  workers of the whole machine, and Studio asks before starting more workers than free cores.
- Files: open by content (session, series, fit run, FID, scan folder; drag and drop, recent list), import an
  averaged FID or a scan folder, save / open sessions (.zulfstudio), export parameters, couplings, lines,
  spectrum, the applied fit and the plot.

Key algorithms:
- Live simulation from the same spin physics as W3 (model cached per structure), drawn as complex Lorentzians
  a gamma / (gamma + i (f - f_k)) with one decay rate: a quick look, not the processed forward model of a fit.
  With an applied, unchanged fit the plot draws that fit's own forward model at fit.json `z_final` instead.
- Display scale by least squares on magnitudes (can be locked, so it does not jump while sliding).
- Auto phase "model": delay on a 0.005 ms grid, phase0 in closed form, so the data best match a positive multiple
  of the simulation; "data": the W2 data-only estimate.

## W5. Blind analysis of a new sample

```mermaid
flowchart LR
    fid["FID"] --> p["standard processing<br/>overview figure"]
    p --> g["group candidates per band<br/>X-Hn line patterns"]
    g --> h["fragment enumeration<br/>covers, trees, motifs;<br/>labelling with a motif"]
    h --> s["search: refine every hypothesis<br/>common chi2 yardstick, BIC,<br/>checks, extension moves"]
    s --> r["blind.md / labelings.json<br/>ranked table"]
    r -- best candidates --> w3["W3 / W4 refined fit<br/>reported couplings"]
```

    python scripts/analyze_sample.py FID.npy --id SAMPLE [--workers 4]                      # blind search
    python scripts/analyze_sample.py FID.npy --id SAMPLE --structure '{"motif": "ethyl", ...}' \
        [--labeling unknown]                                                                 # candidate motif

- Labelling (D54, with `--structure`): `natural` (default), `15N`, `2H-exchange`, or `unknown` (every applicable
  labelling, ranked on one scale).
- Results are conditional candidates: report the ranking, margins and flags, not one assignment.
- The output is a ranking of hypotheses (blind.md, blind.json), not a coupling table: the search uses coarse
  settings to compare many structures. Refit the best candidates with W3 / W4 to get couplings to report.
  Feeding a structure-free J network from this search into W6 is planned (PLAN Phase 7).

Key algorithms:
- Group candidates: an isolated X-Hn group has lines at fixed multiples of 1J (XH: J; XH2: 3J/2; XH3: J and 2J,
  0.8 : 1); every band peak taken as every pattern line fixes a J, and the other lines are checked as observed,
  unobservable or missing.
- Fragment enumeration: smallest covers of the bands by group candidates, then every tree over the groups
  (Pruefer sequences) plus isolated parts and symmetric methyl doubling; longer-range couplings from a prior by
  bond distance.
- Search: each hypothesis refined in amplitude variants (free, fixed, ratios) and exchange regimes; one common
  yardstick (chi2 on the data cores) and BIC = chi2 + k ln N; checks per result; extension moves accepted only if
  the BIC improves by at least 6; fits with an isotopologue switched off rank last.
- Global pattern search (zulf_core.solver.search): a phase-insensitive magnitude objective over the whole
  parameter box proposes starts, since local refinement converges only near the right line pattern.
- Labelling (D54): every labelling refitted, then all rescored on one common scale.

## W6. Fitted J network to structure candidates

```mermaid
flowchart LR
    j["observation JSON<br/>units, protons, couplings,<br/>isotopes"] --> a["route A<br/>rule ranges"]
    j --> bj["route B-J<br/>MLP on J"]
    j --> bk["route B-K<br/>MLP on K"]
    a --> c["three rankings<br/>agree / differ"]
    bj --> c
    bk --> c
```

    python scripts/j_structure_benchmark.py            # once: trains runs/models/j_edges_J.json and _K.json
    python scripts/j_structure.py OBS.json --explain 1 [--out ranking.json]

- All three routes run every time (D56). With 2H or 15N in the sample, B-K is the isotope-independent one.
- A network from a fit of the true structure is biased toward it; see docs/J_TO_STRUCTURE.md.

Key algorithms:
- Candidates: every connected graph over the carbon copies and up to two unseen atoms X, with degree caps, at most
  one ring and up to two double or triple bonds; copies of a unit kept equivalent and distinct units distinct by
  Weisfeiler-Lehman colours; isomorphic graphs merged.
- Score: sum over couplings of log p(J | bond count) + log prior (unseen atom -1, ring -2, double bond -1.5,
  triple -3, missing valence -3); bond count = C-C distance + 1 (C-H) or + 2 (H-H); 1J enters through the
  hybridization of its carbon.
- Route A: literature ranges convolved with the fit sigma. Route B: MLP p(bond count | J, kind, proton counts)
  turned into a likelihood by Bayes, p(class | J) / p(class); mode K feeds the equivalent 13C-1H coupling
  J g_13C g_1H / (g_A g_B), so it is unchanged under 1H -> 2H.

## W7. Learned candidate model (built, not yet trained)

```mermaid
flowchart LR
    g["generator<br/>graphs, J rules, isotopologues"] --> r["render pipeline<br/>noise, processing"]
    r --> sh["prerendered shards"]
    sh --> t["train<br/>CNN set / CNN+Transformer"]
    t --> pr["propose candidates"]
    pr --> ref["refine with zulf_core<br/>held-out ranking"]
```

    zulf-model prerender configs/run_cnn_set_phased_v1.json runs/shards/phased --count 200000 --workers 8
    zulf-model train configs/run_cnn_set_v1.json
    python scripts/smoke_pipeline.py                  # generate -> train -> propose -> refine, small

- Status: code and CPU verification runs only; the full training and the frozen test sets are open
  (docs/PLAN.md). Not part of any result yet.

Key algorithms:
- Encoder: a 1D CNN that downsamples the spectrum into tokens, each with Fourier features of its absolute
  frequency (frequency is physical information in ZULF spectra).
- Set model (DETR style): component slots cross-attend to the tokens; trained with Hungarian matching between
  slots and true components; groups predicted in canonical order, so magnetic equivalence is part of the output.
- Transformer: decodes the token grammar BOS (COMP weight groups SEP J tokens)+ EOS by beam search restricted to
  the grammar; soft targets over neighbouring J bins and a within-bin regression head.
- Evaluation always matches interpretations by permutation (`spinsystem.best_permutation`).

---

## W3 in detail: fitting a known structure (processed or raw ZULF spectra)

This is the working record of how a ZULF spectrum (or a concentration series) is turned into fitted couplings
with honest uncertainties, as practised on the Blake pyridine series and on isopropylamine (October 2026). It
names every step, the algorithm behind it, the code that does it, the settings that worked, and the pitfalls that
were found. Each analysis keeps its own log under docs/analysis/; this file is the method. Update it when a step
changes (AGENTS.md: keep the workflow current).

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

### 1. Diagnostics and processing (per data set)

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

### 2. Phase

Code: zulf_processing.phase (calibrated_phase), zulf_core.render.phasing; skill skills/zulf-phasing.

- First order: the delay at this FID's own switching edge (plus the instrument offset).
- Zero order: the instrument calibration of the sequence (configs/confirmed_samples.json, "phase_calibration")
  or a data-only estimate by peak symmetry with a model-spectrum bias check.
- In a complex fit the shared-phase gain refines a residual zero-order phase, and a model delay (per spectrum or
  shared) is fitted; with the data phased at the edge the fitted delay should stay near the edge (pyridine:
  -1.9 to -4.5 ms against -3.2 ms; a free per-spectrum delay that jumps signals another problem).

### 3. Structure to spin-system model

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

### 4. Forward model

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

### 5. Objective

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

### 6. Parameterization

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
- Coupling precision (`--precision`, D55; default 0.01 Hz): a fit stops when every coupling moves by less than this
  and the objective no longer improves; 0.001 for fine structure, 0 to run to the tolerances as before. Sets the
  decimals of J_table.csv.
- Static field (`--fit-field`, D47; default zero field): for measurements in a nonzero field (fringe field,
  residual shield field). Two parameters, `field_transverse_ut` and `field_z_ut` (>= 0, `--field-bounds`); the signal
  is stationary at zero field, so give nonzero starts (`--field-start transverse,z`) and try several, since the
  coupling multistart does not perturb spectrum parameters. Report the fitted field with the result (D19).
  `scripts/line_table.py` still lists zero-field lines.
- Unknown heteronucleus (D53): with a known field (`--field transverse,z`, held fixed) free its gamma
  (`--fit-gamma 13C --gamma-start ...`); fit.json `gamma_identification` names the nearest registered nucleus.
  Start from several candidate gammas (13C 10.71, 15N -4.32, 31P 17.24 Hz/uT); compare objectives.

### 7. Search

Code: scripts/fit_joint_series (main, _solve_start, _coordinate_scan).

- Multi-start: seeds (earlier solutions, literature, other groups' values), perturbed starts (`--spread`), starts
  drawn from the priors; `--workers` processes (one BLAS thread each).
- Nested starts for nested models: start the slow-exchange model from the fast solution with the exchangeable
  couplings at 0 (isopropylamine: generic slow starts ended 26 % above the fast fit, the nested start 13 % below).
- Branches: a series can sit in two near-equivalent coupling branches; start a run inside each branch
  (pyridine: branch A started from the single-spectrum fits gave a 22 % lower total than every earlier run).
- Optional coarse-to-fine smoothing schedule and coordinate grid scans for barriers along single couplings.
- Ranking always on the hard objective (hard missing-peak rows; residual-peak rows on one common window set).

### 8. Fine structure: from a misfit to the couplings that cause it

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
8. Watching a fit (`scripts/fit_monitor.py runs/processed` -> http://127.0.0.1:8770, or `--text` in a terminal; in
   Studio: the Monitor tab, opened by a double-click on a job, Live monitor in Jobs or Cmd+Shift+M).
   fit_joint_series records every run in OUT/monitor (on by default, `--monitor off`): the objective of every
   residual evaluation of every start with its stage and its parameter vector (start / smoothing / fit / coordinate scan / residual-peak
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

### 9. Reliability and uncertainty budget

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

### 10. Report

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

### Typical command (one spectrum, complex, rate families, fine structure)

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
