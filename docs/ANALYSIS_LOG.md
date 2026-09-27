# Analysis log

Chronological record of experimental analyses and the solver changes they
prompted. Experimental files are not in the repository; entries give the
upload name, processing, numbers and conclusions. Reusable lessons go to the
skills under `skills/`; this file keeps the history.

## 2026-09-26 Blind tests on NMRduino averaged FIDs (4000 Hz, 65516 points)

### 7dc9a043
- No J-coupled signal in 20-520 Hz with any window; instrument lines only.
- Later used as the no-signal reference for background checks.

### 7ad4aafb
- One band 148-192 Hz; pyridine 13C2/13C3/13C4 (literature-like starts)
  0.144 vs benzene 0.579 (complex, full record). Blind search from the model
  skeleton: 0.212 with unphysical couplings. Answer not revealed.

### e3d282da (confirmed: pure L-alanine in water)
- Lines 129.2/130.5, 144.5/146.0/148.5, 258.8 Hz, SNR 4-15.
- Window 0.05-1.05 s (zero fill 4): CH-CH3 fragment H7 0.358 (no-lines 0.711,
  floor 0.37); 1J 145.4 / 129.8, 2J -4.1 / -4.5, 3J(HH) 6.6 Hz; isopropyl
  rejected by abundance; methyl-acetate type 0.617.
- 15N slow-exchange improvement withdrawn: the 69-75 Hz feature is present in
  the no-signal dataset. Correct: NH3+ exchanges fast in water near pH 6
  (zwitterion; neutral -NH2 fraction about 2e-4).
- Phased real-only fits exposed a ripple of period 1/(crop + delay) from the
  record-start baseline; fixed by modelling the background in the record frame
  (0aeb959). With the fix, H7 remained the consistent model.
- Skills split by category and derivation paths recorded (c50b700, 6115e4a).

### Solver changes prompted
- 0aeb959 background in the record frame for phased spectra.
- d5543fc signal-focused weighting (band_weighting="signal"): noise-scaled
  chi-square, data-driven signal mask, signal_region_residual in results.
  On e3d282da the mask found 127-133.5, 142-151, 256.75-261.25 Hz and also
  72.25-76.5 Hz. Correction (user): that feature is not established as
  background. No 13C isotopologue of H7 has lines in 55-100 Hz; a 15NH3+
  isotopologue would put its J line near 73 Hz (the earlier H7 + 15NH3 fit
  gave 1J(N,H) -73.3 Hz and a 15N amplitude 0.42 vs 0.35 expected). The
  no-signal comparison matched only the mean level, not the shape, and
  intermediate NH3+ exchange would broaden rather than remove 15N lines.
  Status: unresolved. Decisive tests: water blank on the same setup; acidic
  sample (slower exchange sharpens 73 / 146 Hz); D2O (ND3+ moves the 15N line
  to about 11 / 22 Hz); shape-wise comparison with other datasets.

### Window and apodization (in progress)
- Peak SNR at a 4 s window is highest without apodization or with about
  0.3 1/s; 4 1/s (the H7 fitted rates) lowers it (130.5 Hz: 13.5 -> 7.3).
  The true line decay is slower than the H7 rates (2.8-4.8 1/s); those rates
  probably absorb unresolved fine structure (the 1.3 Hz doublet at 129-131 Hz).
- Running: H7 refined at windows 0.3/0.5/1/2/4 s, hard truncation and
  0.6 1/s exponential weighting; signal-weighted comparison of H7, H7 + 15NH3
  and self-consistent alanine.

### Signal-weighted re-evaluation of H7, H7 + 15NH3, self-consistent alanine
- Setup: complex, 0.05-1.05 s, zero fill 4, band_weighting="signal",
  background order 1, 3 starts.
- Mask 72.25-76.5, 127-133.5, 142-151, 256.75-261.25 Hz:
  H7 signal-region 0.378 / overall 0.452, Cb/Ca 0.39, no bounds;
  H7 + 15NH3 0.266 / 0.553, amplitudes 1 : 0.31 : 1.33, 1J(N,H) at bound;
  self-consistent 0.239 / 0.449, amplitudes 1 : 0.07 : 0.05.
  The 15N models gain on the unresolved region with abundance ratios far from
  1 : 1 : 0.35; not evidence. H7's Cb/Ca 0.39 points at the methyl-carbon lines
  (130 Hz doublet, 259 Hz) being modelled poorly.
- Solver loophole found: with a data-only mask and outside weight 0.2, a
  model can place spurious lines outside the mask cheaply (seen as strong
  predicted lines at 263-265 Hz that the data do not show). Planned fix: add
  the model-predicted line regions to the mask (iterate fit -> union mask ->
  refit until stable), so every predicted line is fully penalised.
- Process lessons (skills updated): never kill by command-line text; one queue
  launcher script; set module-level ranges after imports and print them.

### Solver: smooth signal weights and model-predicted cores
- User suggestions: focus the solver on signal; taper the weight around
  high-weight peaks instead of a hard mask edge.
- Implemented: weights 1 at peak cores with a Gaussian fall-off (2 Hz) to 0.2;
  a second pass adds model-predicted lines to the cores (closes the loophole
  seen at 263-265 Hz). Tests: smooth monotone fall-off; a 13CH3 candidate on a
  13C-H spectrum triggers the model pass for its unmatched 2J line.

### Window scan of H7 (signal weighting with the earlier hard mask)
- 1J robust: 145.3 +- 0.1 (C-alpha) and 130.0 +- 0.3 Hz (C-beta); for windows
  >= 1 s: 3J(H-alpha, H-beta) 7.0 Hz, 2J(C-alpha, H-beta) -4.4 Hz;
  2J(C-beta, H-alpha) unstable (-4.2 to -6.8 Hz).
- Cb/Ca 1.10 at 0.3 s (the expected 1 : 1) falling to about 0.45 at >= 1 s:
  H7 lacks couplings that shape the methyl-carbon lines at high resolution;
  the fit compensates by lowering the C-beta amplitude.
- Weak matched apodization (0.6 1/s): lower signal-region residual (0.25 vs
  0.29-0.36) and identical parameters for 1, 2 and 4 s windows; recommended.
- Next: re-evaluate H7, H7 + 15NH3 and self-consistent alanine with smooth
  weights and model-predicted cores, 1 s window, 0.6 1/s, full ranges
  (running).

### Re-evaluation with smooth weights and model-predicted cores (full ranges)
- Complex, 0.05-1.05 s, 0.6 1/s apodization, ranges 62-320 Hz (instrument
  lines only removed). Data cores: 74.25-74.5, 129.0-129.5, 130.25-131.5,
  144.0-144.75, 145.5-146.5, 148.25-149.0, 258.75-259.25 Hz.
- H7: signal-region 0.329, Cb/Ca 0.43, no bounds; 1J 145.40 / 130.28,
  2J -4.41 / -6.81, 3J(HH) 7.03 Hz (model pass +5 points).
- H7 + 15NH3: 0.316, 1 : 0.52 : 0.24, 1J(N,H) at bound (-84.8 Hz),
  3J(HN, H-alpha) 1.9 Hz (inconsistent with slow exchange) (+18 points).
- Self-consistent alanine: 0.288 but amplitudes 1 : 0.05 : 0.04 (C-beta and
  15N suppressed, one broad C-alpha component, rate 8.5 1/s) (+23 points).
  Its earlier reproduction of the 130 Hz doublet relied on spurious lines
  outside the old mask; with model-predicted cores that solution is gone.
- Verdict: H7 remains the physically consistent model; NH coupling not
  established; 72-76 Hz unresolved. Next solver step: amplitude ratios fixed
  by natural abundance (pure sample 1 : 1 : 0.35), so models cannot collapse
  onto one isotopologue.
- Check of the model pass (user question: are model-predicted peaks
  up-weighted?): yes, but only predicted narrow excess above 4 sigma joins
  the cores (H7 +5, H7 + 15NH3 +18, self-consistent +23 points, mostly edges
  of existing cores: 132, 146.75, 258.5-260 Hz; self-consistent also 71.75,
  132.25-132.5, 260.25-260.5 Hz). Not covered: weaker spurious lines (the
  264 Hz line of H7 + 15NH3 keeps weight about 0.26), predicted dips below
  the baseline (138.5 Hz in the self-consistent fit), and lines that appear
  only after the second pass (one pass, not iterated). The figure's yellow
  shading showed data cores only. Planned: model threshold relative to the
  noise but lower (model has no noise), absolute deviation from the
  background instead of positive excess, iterate until the cores are stable,
  plot the weights actually used.

### Solver: model lines straight from the transition lists
- User: the model already gives the line positions; use them instead of
  rendering the spectrum and peak-picking again.
- Implemented: `MixtureForward.model_line_heights` gives each line's height
  abs(g_c a_k) times the rendered peak of a unit line with the same rate
  (lines within half a grid step add up); lines above 2 sigma of either sign
  join the cores; passes repeat up to 3 times until no new cores. Results add
  `data_region_residual` (data cores only, same region for every model) and
  `signal_model_lines_hz`.
- On the saved e3d282da fits this finds lines the peak picking missed:
  H7 143.0 (7.4 sigma) and 266.0 Hz (8.5); H7 + 15NH3 also 148.0 (16.7, its
  2J(N,H) line) and 263.75 Hz (11.4); self-consistent alanine 69.75-76.25,
  132-132.5, 147-147.75, 255-267 Hz. Several are lines that partly cancel in
  the rendered spectrum.
- Re-evaluation with transition-list cores (9e46dbd; complex, 0.05-1.05 s,
  0.6 1/s, 62-320 Hz, 72-76 Hz included; two model passes each):
  H7: data-region 0.335, Cb/Ca 0.39, no bounds; 1J 145.41 / 130.28,
  2J -4.35 / -6.56, 3J(HH) 6.95 Hz; model lines 143, 250, 266 Hz joined.
  H7 + 15NH3: 0.321, 1 : 0.57 : 0.29, 1J(N,H) at bound (-84.7 Hz); its 15N
  lines at 78-91 and 162-177 Hz joined the cores where the data show nothing;
  the 74 Hz feature is not explained.
  Self-consistent alanine: 0.249, 1 : 0.07 : 0.04 relative to C-alpha, but in
  absolute terms C-alpha grew to 7.5x the H7 value with rate 7.7 1/s: a broad
  component that acts as background (it makes the 139 Hz dip; its broad lines
  stay below 2 sigma each, so no core there). C-beta : 15N = 1 : 0.68 with
  narrow lines (1.1 and 1.4 1/s) and 1J(N,H) -74.07 Hz, putting a 15N line on
  the 74 Hz feature (expected 15N : 13C 0.35).
- Conclusion: model-line cores work as intended (spurious lines penalised);
  the remaining failure mode is a component turning into broad background.
  Next: fix amplitude ratios to natural abundance and bound decay rates (or
  tie rates across isotopologues) so no component can become background.
- Why the 139 Hz dip of the self-consistent fit survives: the C-alpha
  component (rate 7.7 1/s, 7.5x the H7 amplitude) is about 50 overlapping
  broad lines; each is only about 0.2 sigma high, so none joins the cores,
  but together they carry about 9 sigma of smooth "signal" over 125-155 Hz,
  taking over the background (background term 0.0024 vs 0.006 in H7). Its
  lines at 139.04 / 139.61 Hz interfere destructively with the rest and cut
  a hole: misfit 8.9 sigma at 139.25 Hz, but weight 0.25 there, so it costs
  only 2.2 sigma. The per-line criterion misses broad components. Options:
  (1) cores from the incoherent line envelope sum_k abs(g a_k) profile_k(f)
  (catches broad overlapping lines and cancellation holes); (2) physical
  constraints: amplitudes at natural abundance, bounded or tied decay rates.

### Solver: fixed amplitude ratios; constrained self-consistent alanine
- Added `RefineSettings.amplitude_ratios` / `MixtureForward(amplitude_ratios=)`:
  components are merged into one column sum_c r_c col_c before the linear
  solve (one overall gain); the Jacobian merges the derivative columns the
  same way (tested against finite differences for shared-phase and complex
  gains; refinement keeps the ratio exact).
- User: add both constraints, run only the self-consistent alanine. Running:
  ratios 1 : 1 : 0.34, one decay rate tied over C-alpha, C-beta, 15N; three
  starts in parallel (candidate, earlier fit with the C-beta rate, random).
- Result (2a16dfe; ranges 62-320 Hz, 72-76 Hz included; three starts, about
  10 minutes in parallel): best data-region 0.290 (candidate start), 0.295
  (earlier fit), 0.376 (random); no bounds; shared rate 2.2 1/s. For
  comparison on the same data cores: H7 0.335, H7 + 15NH3 0.321, free
  self-consistent 0.249 (collapsed, unphysical).
  Robust between the two good starts: 1J(Ca,Ha) 145.0 / 145.9, 1J(Cb,Hb)
  129.4 / 130.3, 1J(N,H) -73.7 / -73.8 Hz (the 15N line at natural abundance
  lands on the 74 Hz feature). Not determined: 3J(Ha,Hb) 8.5 / 6.0,
  3J(Ha,HN) 10.4 / 9.0, 2J(Ca,Hb) -5.1 / -2.4, 2J(Cb,Ha) -6.7 / -7.4 Hz.
  The 139 Hz hole is gone; remaining misfits: 130.8 Hz and 148.5 Hz peaks too
  low, a broad shoulder at 131-134 Hz and a line at 251.5 Hz the data do not
  show; 259 Hz height now matches (H7 overshot it).
- Tension: 1J(N,H) lines need slow NH3+ exchange, while alanine in neutral
  water exchanges fast; the 74 Hz match is therefore suggestive, not proof
  (pH of the sample unknown; decisive tests as listed above).

### Solver: incoherent line envelope for the model cores
- User: start the change (option 1 above). `MixtureForward.model_envelope`
  replaces the per-line heights: E(f) = sum_k abs(g_c a_k) P_R(f - f_k) with
  the rendered unit-line magnitude profile, cut at 20 Hz. Tests: equals the
  rendered magnitude around an isolated line; equals the sum of component
  magnitudes for two cancelling lines while the coherent sum drops to about
  0.33 of it.
- On the earlier free self-consistent fit the envelope is 4.4-4.8 sigma over
  136-140 Hz, so the 139 Hz hole would now be fully weighted; 32 % of the
  fitted points fall in its model cores (broad C-alpha component).
- Running: constrained self-consistent alanine (1 : 1 : 0.34, shared rate)
  with envelope cores, three starts in parallel.
- Result with envelope cores (6595bad): data-region 0.295 / 0.330 / 0.354
  (candidate / earlier fit / random), overall 0.359 (0.410 with per-line
  cores); shared rate 2.6 1/s; the random start hit the phase-delay bound.
  The dips at 137 and 141 Hz are gone. Remaining: 130.8 and 148.5 Hz peaks
  too low, a line at 251.5 Hz the data do not show (now at full weight), a
  76 Hz dip deeper than the data.
  Robust: 1J(N,H) -73.6 / -73.8 Hz; 1J(Cb,Hb) 129.4 / 130.3 Hz. Less so:
  1J(Ca,Ha) 144.7 / 146.1 Hz. Not determined: 3J(Ha,Hb) 8.3 / 5.8,
  3J(Ha,HN) 10.7 / 8.4 Hz.
- Side effect: the envelope (2 sigma, profile tails to 20 Hz) covers most of
  62-88, 125-152 and 250-270 Hz, about 350 points, so there the weighting is
  close to uniform and the focus on peaks is reduced. Options: a higher
  model threshold, or a shorter profile cut.
- User: release the fixed amplitude ratios. Running: self-consistent alanine
  with free amplitudes, one shared decay rate, envelope cores; starts:
  candidate, best constrained envelope fit, random.
- Result (free amplitudes, shared rate, envelope cores): data-region 0.271
  (start: constrained fit), 0.292 (candidate), 0.330 (random, phase-delay
  bound). Amplitudes C-alpha : C-beta : 15N = 1 : 1.54 : 0.80,
  1 : 0.56 : 0.54, 1 : 1.34 : 1.01; shared rate 2.3-3.4 1/s.
  C-alpha : C-beta is not determined (0.56-1.54, around 1). 15N relative to
  the mean 13C site is 0.63-0.86 in every start, about twice the natural
  0.34: the 74 Hz feature is stronger than a 15NH3+ line at natural
  abundance, consistent with an extra (for example instrumental) part or
  with a 15N model that is not right.
  Robust: 1J(Ca,Ha) 144.5-145.4, 1J(Cb,Hb) 129.1-130.3, 1J(N,H) -73.4 to
  -74.0 Hz. 3J(Ha,Hb) 7.0 / 8.5 / 9.0 Hz (candidate start 7.04, close to the
  literature value about 7.2); other small couplings not determined.
  Best fit matches 142-147 Hz closely; still low at 130.8 and 148.5 Hz and
  shows the 251.5 Hz line; the candidate-start fit has 130.8 Hz right but
  misses 129.2 Hz and is low at 259 Hz.

## 2026-09-26 Blind test b683220d (NMRduino averaged FID, 65516 points, 4000 Hz assumed)
- Diagnostics: first-point anomaly, saturation plateau to 2.75 ms, ringing to
  44.75 ms; standard processing (crop 0.05 s, SG 201 / 2, mean removed).
- Lines (0.05-1.05 s, 0.6 1/s, SNR vs 300-500 Hz noise): 128.4 (7.5),
  130.0 (9.9), 131.8 (5.2), 144.3 (5.2), 145.3 (5.0), 146.4 (9.5),
  147.75 (16.7), 150.0 (8.3), 151.9 (5.2), 257-258.3 (3.4), 262.9 (2.3) Hz.
  Instrument lines as before (60 Hz harmonics, 294, 922.9 Hz).
- Same band layout as e3d282da (L-alanine) but shifted: the 130 Hz group
  about 0.8 Hz lower (128.4 / 130.0 vs 129.25 / 130.75), the 146 Hz group
  about 1.5-2 Hz higher (146.4 / 147.75 / 150.0 vs 144.4 / 146.0 / 148.6),
  259 Hz group about 0.8 Hz lower; lines about 2x stronger. No feature at
  72-76 Hz (e3d282da had 4.3-4.5 there) and nothing above noise at 60-100 Hz
  besides instrument lines.
- First reading: a CH-CH3 fragment (H7 type) with 1J(CH) about 147 Hz and
  1J(CH3) about 129 Hz, i.e. slightly different substituent or pH from the
  alanine sample. Not yet fitted.
- Running: H7 (CH-CH3, both single-13C isotopologues) on b683220d, complex,
  0.05-1.05 s, 0.6 1/s, 62-320 Hz, signal weighting with envelope cores;
  variants: amplitudes 1 : 1 with a shared rate, and free amplitudes/rates.
  Starts 1J 147 / 129, 2J -4.4 / -4.5, 3J(HH) 7.0 Hz.
- H7 result on b683220d (about 10 minutes, no bounds in either variant):
  free: data-region 0.280, amplitudes 1 : 1.045, rates 1.9 (CH) / 4.7 1/s
  (CH3 carbon); 1J 147.12 / 129.19, 2J(Ca,Hb) -4.38, 2J(Cb,Ha) -4.07,
  3J(HH) 6.52 Hz. Fixed 1 : 1, shared rate 2.9 1/s: 0.344, same 1J, 2J
  -4.37 / -3.95, 3J(HH) 6.38 Hz.
  Unlike alanine (Cb/Ca 0.39-0.43) the free amplitudes come out at the
  natural 1 : 1 by themselves. The 146-151 Hz group is matched closely; the
  methyl-carbon lines (128.4 / 130.0 doublet, 257-263 Hz) are again the weak
  part: the free fit broadens them (4.7 1/s) and the fixed fit overshoots at
  258 Hz. Same pattern as alanine: H7 lacks a small splitting on the
  methyl-carbon lines.
  Compared with alanine (H7, same settings): 1J(CH) 147.1 vs 145.4, 1J(CH3)
  129.2 vs 130.3, 2J(Cb,Ha) -4.1 vs -6.6, 3J(HH) 6.5 vs 7.0 Hz. No 72-76 Hz
  feature: a 15NH3+ in slow exchange (for example alanine at low pH) would
  show one; points to a fragment without slowly exchanging N-H.
- User: test the hypotheses. Running against H7 (same data, settings,
  weighting): H8 = H7 plus one slowly exchanging proton on X (CH(XH)-CH3,
  e.g. lactic-acid OH), ISO = isopropyl CH(CH3)2 (natural ratio 1 : 2);
  each free and with natural ratios plus one shared rate.
  Only two protonated-carbon bands are present, so compounds with a third
  protonated carbon (ethyl, aldehyde) are excluded before fitting.
- Hypothesis results (data-region residual, same data cores):
  H7 free 0.280 (1 : 1.04); H8 free 0.190 (1 : 0.87, rates 1.8 / 2.5 1/s,
  no bounds); H8 fixed 1 : 1 + shared rate 0.283; ISO free 0.352 with
  amplitudes 1 : 1.16 (natural 1 : 2), ISO fixed 0.366 -> isopropyl
  rejected (worse fit, abundance contradicts).
  H8 free reproduces the 128.4 / 130.0 doublet, the 144-152 Hz group and the
  split 256.7 / 258.3 Hz structure; the methyl-carbon lines no longer need
  extra broadening (H7 needed 4.7 1/s). Misfit on the data cores falls to
  0.46 of H7's with 4 added couplings.
  H8 free: 1J 147.07 / 129.18, 2J(Ca,Hb) -4.69, 2J(Cb,Ha) -4.40,
  3J(Ha,Hb) 7.06 Hz; extra spin: 2J(Ca,XH) 0.55, 3J(Cb,XH) 1.51,
  3J(Ha,XH) -1.21, 4J(Hb,XH) -1.76 Hz. H8 fixed found a different set
  (-6.94, 4.28, 8.25, 1.37 Hz): the extra spin's couplings are not
  determined, and the free-fit values (negative 3J(Ha,XH), a 1.8 Hz 4J) are
  not typical of a CH-OH proton.
  Conclusion so far: one more weakly coupled spin-1/2 on the CH-CH3 fragment
  is required by the data; its identity is not yet established.
- Final guess (user asked): lactic acid / lactate, CH3-CH(OH)-COOH.
  Reasons: exactly one CH and one CH3 protonated carbon (two bands,
  amplitudes 1 : 1, isopropyl rejected); 3J(HH) 7.06 Hz as for a CH-CH3
  unit; 1J(CH) 147.1 Hz, 1.7 Hz above alanine on the same instrument, as
  expected for O rather than N on the CH; no 15N line at 72-76 Hz; one more
  weakly coupled spin-1/2 needed (H8), compatible with the OH proton, though
  its fitted couplings are not determined. Runner-up: alanine at another pH
  (disfavoured by the larger 1J(CH) and the absent 74 Hz line); other
  CH(X)-CH3 acids (for example 2-halopropionic) less likely.
- Confirmed by the user: lactic acid. Lessons added to the blind-analysis
  skill and case notes (same-instrument overlay with confirmed samples,
  isotopologue-rate asymmetry as a missing-coupling marker, extra spin
  reported as required but unidentified).

### Programmatic hypotheses, stage 1 (D34)
- User: can the chemical-hypothesis experience be programmed; leave
  interfaces for future cases. Added `zulf_core.hypothesis` (fragments with
  symmetry, builder, checks, moves, knowledge base with L-alanine and lactic
  acid). Replay on b683220d (runs/blind/b683220d/replay_builder.py): builder
  H7 / H8 reproduce the hand-built predictions (max relative difference
  2e-13); checks flag H7 rate asymmetry (C-beta 2.5x) and the isopropyl
  abundance (1.16 vs 2). The H8 free / fixed pair is not "comparable" under
  the 10 % residual rule (0.190 vs 0.283), so no undetermined flag.

### Programmatic hypotheses, stage 2 and neural hints (D35)
- User: continue; add an interface for neural-model guesses as an insight
  provider between the stages. `propose_hypotheses` on the confirmed samples
  (runs/blind/model_hint_insight.py): lactic acid #1 bonded 13CH3-13CH
  (J 128.5 / 147.8, all bands explained), #2 isopropyl variant; alanine #1
  bonded 13CH3-13CH (J 129.2 / 146.0) with the 74 Hz band open, 15N
  alternatives further down. Iterations needed on the way: partner lines
  scored against their predicted height (a 15NH3 2J line had "explained" the
  whole 146 Hz band), a 15N pattern prior, a multiplet width for explained
  bands (resolved methyl 2J multiplets at high SNR), a noise floor for
  near noise-free input.
- Neural hints: transformer checkpoints give no one-bond couplings; the set
  model gives 7 group hints per sample, 1 agreeing (147.6 / 131.9 Hz), 6 at
  frequencies with no data lines. Hint weight lowered to 0.1 after one weak
  agreeing hint nearly lifted CH-CH above CH3-CH on alanine.

### Minor isotopologues and labelling schemes (D36)
- User: doubly labelled molecules are much rarer, weight them down; do not
  hard-code, natural abundance by default with an enrichment mode. Weight =
  exact label-set probability of a `Labeling`; minor sets follow a parent
  (amplitude map, tied rate); SNR gate 2 / peak SNR. On b683220d the
  13C-13C set is omitted (0.0108 < 0.0503); rankings unchanged.

### zulf_hypothesis split out (D37)
- User: do not mix the hypothesis system with the refiner and the models.
  Moved to the top-level package `zulf_hypothesis` (adapters for neural
  hints inside it, lazy); boundary tests added; real-data rankings unchanged.

### Stage 3 search on the confirmed samples (D38; runs/blind/search_real.py)
- b683220d (lactic acid), 5.4 min on 4 workers, no hand input: proposals ->
  CH3-CH refined (free: rate_asymmetry + misfit) -> move "coupled proton on
  X" triggered and accepted (Delta BIC -610 free, -1025 fixed) -> best
  CH3-CH + HX [free], amplitudes 0.87 : 1; knowledge: lactic acid 0.32 Hz RMS,
  L-alanine 1.32 Hz. Isopropyl and CH3-CH3 variants rank below with bounds,
  abundance and background findings. The automatic path reproduces the
  hand analysis (H7 -> H8, isopropyl rejected, lactic acid).
- Yardstick on this window: 26 decimated core points (N = 52), small; the
  BIC differences are large compared with the penalty terms anyway.
- e3d282da (L-alanine), 2.6 min: best clean CH3-CH [free] (amplitudes
  0.64 : 1); the extension + HX has the lowest BIC (Delta -187) but an
  abundance warning (0.55 : 1), so it is not chosen as best; the 74 Hz band
  stays open (15N alternatives were outside the top 3 proposals).
  Knowledge: L-alanine 0.88 vs lactic acid 0.90 Hz RMS - the CH-CH3
  couplings alone do not separate the two compounds on this data (as in the
  hand analysis). Isopropyl ranks below with bounds / background findings.
- Fix: the fixed variant is skipped when no ratio block has two components
  (isolated single-component parts gave duplicate rows).

### Recognition tests beyond CH-CH3 (runs/bench, synthetic via the builder)
- Stage 2 on synthetic structures: isolated methyl correct (#1 13CH3
  J 140.9); 13CH3-15NH3 (slow exchange) correct (#1); acetone-like
  CH3-C(O)-CH3 seen as CH3 only (equivalence and the carbonyl bridge are
  invisible); CH2-CH2 (symmetric) wrong (multiplet split into several
  groups); ethyl and propane-like give no proposal (9-11 bands, a coupled CH2
  splits into many bands). Cause: isolated X-Hn patterns in stage 2, trees
  only (no rings, no heavy-atom bridges), only methyl-leaf equivalence.
- Pyridine (couplings from memory for the reference fragment):
  synthetic truth: proposals miss the ring (#1 single 13CH J=158); with the
  reference fragment added, it ranks first (fixed variant, amplitudes
  1 : 1 : 0.5 as expected, no findings; all proposals > 4.7e5 BIC worse) -
  scoring works, proposal generation does not.
  Real 7ad4aafb: 6 bands (150.5-185.2 Hz); proposals are 3-CH chains; the
  reference pyridine wins by > 4.7e4 BIC over every proposal. Best clean:
  pyridine [fixed]; the free variant has a lower BIC but amplitudes
  0.86 : 0.79 : 1 against 1 : 1 : 0.5 (C4 too strong) - the from-memory
  couplings or the model are not right in detail.
- Next: motif library scan (coupled spectra of registered motifs, 1J from
  band positions, quick fits) so rings and coupled CH2 enter the proposals;
  a fixed synthetic benchmark of these structures.

### Motif scan and recognition benchmark (D39)
- Registered motifs incl. slow-exchange amines (user: add amines). Benchmark
  (SNR 60, natural abundance): motif rank 1 in 10 / 12 (CH3-15NH2 second,
  CH3-C(O)-CH3 third); group path rank 1 in 2 / 12.
- Real 7ad4aafb: the motif scan alone ranks the pyridine ring first (generic
  aromatic couplings, 1J 167 / 167 / 161 from band positions), 3.9e4 BIC ahead
  of the benzene ring.
- Cost note from the real-data run: the coupled-proton move proposed one
  extension per aromatic C-H of pyridine (5 x 2 variants of 7-spin
  isotopologues, tens of minutes). Ring atoms are now skipped as anchors
  (no free valence); CH-CH3 still gets its one extension.

### Full pipeline with motifs on the real samples (runs/bench/search_real_motifs.py)
- First run (plain BIC, clean-first best, free-amplitude motif scan):
  pyridine found by the motif scan and ranked first after refinement, but
  "best" was the benzene ring (both pyridine fits had warnings); the scan
  ranked amine motifs first on lactic acid and alanine (15N isotopologue with
  a free amplitude imitating the 147 Hz band); lactic acid unchanged.
- Fixed by D40 (quasi-likelihood scores, clean margin, fixed-ratio scan);
  rerun in progress.
- Rerun with D40: pyridine best = pyridine ring [fixed] (clean, 3.7 from the
  minimum; benzene +78); alanine best CH3-CH [free], extension rejected
  (+1.6, as expected). Lactic acid: extension rejected (+14) because it was
  grown from the motif copy of CH-CH3 and restarted from generic couplings
  (one start, worse local minimum), unlike the previous run. Fix: extensions
  start from the parent's refined couplings (warm start); a motif identical
  to a group proposal is refined once. Rerun in progress.
- Warm-start rerun (b5ab115): pyridine best = pyridine ring [fixed] again;
  alanine best CH3-CH [free], extension rejected (+9 / +11); lactic acid:
  extension accepted only in the fixed variant (-10 vs its parent), the free
  extension rejected (+6). Multi-start extensions (ecd032e: 4 starts, 1.5 Hz
  around the warm start) gave identical numbers.
- Diagnosis: the parent CH3-CH compensates the missing coupling with
  broadened methyl lines (the rate_asymmetry finding); a warm start inherits
  that compensation and stays in its basin. The first run reached the
  better basin (chi2 1153 -> 528) only because its extension started from
  the cold proposal values. Fix (69620d5): extensions start from both the
  warm (parent's refined couplings) and the cold (proposal) values, each
  with perturbations; refine keeps the best. Lactic-acid rerun in progress.
- Lesson (skill): a nested extension is not guaranteed to improve from the
  parent's optimum when the parent has compensated for the missing
  physics; start extensions from both the parent's optimum and the
  uncompensated proposal values.
- Warm + cold starts (69620d5) on lactic acid: unchanged in substance - the
  fixed extension is accepted (-7 vs its fixed parent), the free extension
  rejected (+14). The earliest success (chi2 1153 -> 528) used a base of 2
  starts in runs/blind/search_real.py: the better basin depends on the
  number and placement of starts (several comparable minima for the extra
  proton, as in the hand analysis). Trying 8 extension starts; next option:
  the solver's global pattern search for extension starts.
- Caveat: the knowledge match "lactic acid 0.002 Hz" is circular (the
  reference entry was fitted on this same sample); blind use must exclude
  entries from the sample under test.
- 8 extension starts (default now 8): lactic acid best = CH-CH3 + coupled
  proton [fixed], chi2 700.8 (the same value as the first successful run's
  fixed extension, so a reproducible basin), accepted against its fixed
  parent (chi2 1741, quasi-BIC -42); the free extension reaches chi2 889
  (first run 528) and is rejected by 1.1 on the quasi scale (four more
  parameters). Outcome now matches the hand analysis (H8 over H7). The free
  extension's better basin is still start-dependent; next option if needed:
  global pattern search for extension starts.
- User chose option 2: global pattern search (zulf_core.solver.search) for
  extension starts. `extension_global_search` runs differential evolution on
  the phase-insensitive pattern objective over the small couplings and rates
  (one-bond couplings held at the warm values); its distinct starts join the
  warm / cold / perturbed ones. Lactic-acid rerun in progress.
- The first global-search rerun changed nothing because the search never
  ran: its default tapers (0.3 s rise + 1.0 s end) exceed the 1 s window and
  the exception was swallowed. Fixed: tapers scale with the record (0.15 /
  0.3 of its length) and search failures are logged in the refine step.
  Standalone check on the lactic extension: 41 s, 4 distinct starts.
- Lesson: a helper that "only proposes starts" must still report failures;
  a silent fallback made two runs look like a result.
- Working global search (526a952), lactic acid, 16 min: "global search: 4
  starts" logged for both extension variants; the free extension reaches
  chi2 593 (889 before; best seen 528) and is accepted, the fixed one 700.8,
  accepted. Both variants now favour CH-CH3 + coupled proton, as the hand
  analysis did (H8 over H7).

### Blind sample fde3fbb2 (automatic run only; user: do not take part)
- Standard processing; the automatic run started with doubly labelled
  isotopologues admitted by the SNR gate (peak SNR 266) and was slow (ring
  motifs with about ten 7-8 spin components). User: default to natural
  abundance mode. Pipeline default is now single-label natural abundance
  (doubles opt-in); run restarted with it.
- Automatic result (single-label natural abundance, 9.0 min, no hand input):
  best = motif CH2-CH2 [free] (1J 132.35 Hz, 2J(C,H) -2.92, 3J(H,H) 5.61 Hz),
  only finding "misfit"; then single 13CH2 / 13CH (+103, bounds), pyridine
  ring (+136 free, +184 fixed, bounds), benzene ring (+1154). c_hat 810
  (reduced chi2 of the best fit): large unexplained structure; no extension
  move applies (no methine anchor); no knowledge match (no CH2-CH2 entry).
  Reported to the user as program output only (no interpretation by me).
- Neural checkpoints on fde3fbb2 (user: run the model): the three
  transformers return five identical beams each (verify_clean: 13C + 5 H,
  1J 218.9; verify_cpu: 13C + 5 H, no J >= 1 Hz; noiseless: 13C2 + 4 H,
  1J 247.0) - degenerate beams, no alternatives. verify_cpu (set model)
  gives five different 13C + H2 + H + H systems; its #1 has 13C-H2 124.1 Hz.
  Insight report: one hint agrees with the data (13CH2 J 124.1 vs the data
  candidate 13CH2 J 127.3 in the 186-192 Hz band); 13 of 14 group hints
  predict lines where the data show none. Proposals unchanged by the hints.
- Hand guess (user asked for it after the program run): succinic acid
  (CH2-CH2 between two carboxyl groups), reasoned from 1J(CH) 132 Hz and
  the absence of methyl or methine bands.
- Confirmed by the user: ethylenediamine (H2N-CH2-CH2-NH2). The program's
  core (symmetric CH2-CH2, 1J 132.35 Hz, 3J(H,H) 5.61 Hz) was right; the
  substituent in the hand guess was wrong.
- Lessons: (1) 1J(CH) of a CH2 next to N and next to a carboxyl lie within
  a few Hz of each other; 1J alone does not identify the substituent.
  (2) The evidence for N is 15N isotopologue lines or N-H couplings (slow
  exchange); with fast N-H exchange the amine protons decouple and the
  molecule looks like a bare CH2-CH2. The decisive measurements are a low-pH
  or D2O sample, or a same-instrument overlay with a known diamine / diacid.
  (3) A large c_hat (810 here) and unassigned weak bands (65.2 and 156.5 Hz,
  SNR 5 and 4) mean the model is incomplete; they must be reported as open,
  not filled with a guessed substituent.
- Added to the program: template CH2-CH2 (fragments), knowledge entry
  "ethylenediamine" (sample fde3fbb2, J(C1,H1) 132.35, J(C1,H2) -2.92,
  J(H1,H2) 5.61 Hz), slow-exchange motif H2N-CH2-CH2-NH2 (13C and 15N
  isotopologues, symmetric ends merged), and benchmark case
  "ethylenediamine (fast exchange)" (motif rank 1, found as CH2-CH2).
- Fit inspection (user asked to see the CH2-CH2 fit; runs/blind/fde3fbb2/
  fit_ch2ch2.py, same settings as the auto run): re-refinement stays at 1J
  132.35, 2J(C,H) -2.92, 3J(H,H) 5.61 Hz, one decay rate exp(log_rate0)
  = exp(1.565) = 4.8 1/s, reduced chi2 738 on the yardstick. The model reproduces the main
  band's pattern (lines 187.7, 192.0, 198.5/198.7, 200.2, 204.5 Hz) but
  residuals reach 75-130 sigma at 191-192, 197-199 and 204.5 Hz: the data
  lines are sharper and taller than the single-rate model (the dispersive
  edges at 198.5 and 204.5 Hz are underestimated). Weak features at 65-68 Hz
  and 156-158 Hz (about 5 sigma above the local level) have no model line.
  Outside the cores the model's tail and linear background deviate at
  230-290 Hz (low weight there). Open: unequal line widths (e.g. residual
  N-H coupling under intermediate exchange, or several rates) and the two
  weak bands; not tested yet.

### Blind sample 4322bdfc (automatic run only, same protocol as fde3fbb2)
- Standard processing (crop 0.05 s, SG 201 / 2, mean removed, fs 4000 Hz),
  compared with the confirmed samples. Strong band 121-133 Hz (peaks 123.3,
  125.0, 127.3-127.7 Hz, SNR up to 113), band 235-257 Hz (250 Hz strongest
  in the 0.3 s window), weak structure 190-210 Hz, and features at 15-60 Hz
  (outside the standard ranges, which start at 62 Hz). Automatic run started
  with the fde3fbb2 script unchanged (single-label natural abundance).
- Automatic result (18.7 min, no hand input). Bands: 127.8 (SNR 101, 123.8-
  130.2), 249.2 (79), 255.8 (35), 195.8 (13), 237.2 (10), 73.2 (6), 118.5
  (6). Motif scan first: CH-CH3 (1J 127.8 / 125.0). Search: c_hat 227;
  misfit triggered the coupled-proton move on the CH-CH3 motif, accepted in
  both variants (free -260, fixed -189). Best = "CH-CH3 + HX on X" [free]
  with warnings bounds, abundance (amplitudes 0.061 : 1, expected 1 : 1),
  background_component, rate_asymmetry, misfit; its fixed variant is +64.5
  (undetermined, misfit): 1J 127.30 / 124.98, 3J(H,H) 12.1 Hz. Then motif
  CH3-NH3+ (+238, 1J(N,H) -64.2), plain CH-CH3 (+254 fixed: 1J 127.33 /
  125.13, 2J -4.16 / -5.80, 3J(H,H) 6.52 Hz). No hypothesis without
  warnings. Nearest knowledge entries L-alanine / lactic acid at 8.2 / 8.8
  Hz rms (far). Reported as program output only.
- Confirmed by the user: triethylamine, N(CH2CH3)3. The program's best
  (CH-CH3 + a coupled proton on a neighbour) has the wrong carbon type on
  C1 (CH instead of CH2); the ethyl motif was ranked 6th of 13 in the quick
  scan (1J 130.5 / 123.8 from band positions) and not refined (top_motifs 3).
  The accepted "HX on X" extension (J(C1,HX) 5.9 Hz) points at the real
  extra couplings: in Et3N the labelled CH2 carbon couples through N to the
  four protons of the two other N-CH2 groups (3J(C,N,C,H)).
- A full Et3N fragment (16 spins, split into magnetically inequivalent
  groups by the builder) did not finish one transition calculation in 15
  min. Diagnostic in progress: ethyl motif and a reduced N-ethyl fragment
  (one ethyl + one 4H group on C1 through N) refined on the auto-run data,
  settings and yardstick (runs/blind/4322bdfc/tea_refine.py).
- Diagnostic result (same data, settings, yardstick; n = 102, c_hat 227):
  ethyl free chi2 37370 (k 11; CH2 carbon rate 15.9 vs 3.4 1/s), ethyl
  fixed 74812 (k 9), reduced N-ethyl free 37795 (k 12, stuck near the
  ethyl optimum), reduced N-ethyl fixed 39974 (k 10: 1J 130.98 / 125.01,
  2J -4.41 / -3.04, 3J(H,H) 7.10, 3J(C,N,C,H) 3.13 Hz). On the quasi-BIC
  scale: program best (CH-CH3 + HX free) 156, ethyl free 216, CH-CH3 + HX
  fixed 221, N-ethyl fixed 222, ethyl fixed 371, CH-CH3 motif fixed 410.
  Among fixed-abundance fits the correct reduced structure ties with the
  program's extension (+1); the through-N protons halve chi2 against the
  fixed ethyl. The program's first place rests on a near-empty CH component
  (amplitude 0.061 vs 1 expected), which the abundance check flagged.
- Why the program missed it: (1) the quick scan picks 1J from band
  positions with generic couplings; ethyl got 1J 130.5 / 123.8 and ranked
  6th, below the refinement cutoff (top_motifs 3); (2) no move changes the
  proton count of a carbon (CH -> CH2), so the CH-CH3 motif could not turn
  into ethyl; (3) a free variant that nearly switches off an isotopologue
  still ranks first (flag only); (4) the 13CH2 isotopologue of Et3N spreads
  over many weak lines (195.8 Hz band, SNR 13; possibly features below
  62 Hz, outside the ranges).
- Added: template and motif "N-ethyl (Et3N)" (reduced: other N-CH2 protons
  as one 4H group on a never-labelled pseudo-site), knowledge entry
  "triethylamine" (fixed-variant values above), benchmark case
  "triethylamine (reduced)": motif rank 2 behind an isolated methyl (weak
  13CH2 isotopologue), 11 of 14 cases at rank 1. On the real data the new
  motif ranks 3rd in the quick scan (1J picks 156.8 / 127.8, again from
  band positions).

### Blind sample e66a4b08 (automatic run only, D42 program)
- Standard processing, compared with the confirmed samples. Bands:
  120-135 Hz (strongest 125.0 and 131.5 Hz in the 1 s window, SNR 86 /
  109), 139-153 Hz weak, 188-210 Hz (198 / 203 Hz, SNR 25), 238-272 Hz
  (250 Hz SNR 80, 258-267 Hz SNR 30-40), features below 62 Hz (outside the
  standard ranges). Automatic run started with the standard script (single-
  label natural abundance), with the D42 program (motif screen,
  change_proton_count, demotion); the 5-sample regression rerun was
  reniced to lowest priority meanwhile.
- Automatic result (57.7 min, shared CPU with the regression rerun; no
  hand input). Bands: 125.2 (SNR 43), 131.8 (77), 202.8 (15), 249.8 (52),
  256.2 (12), 262.5, weak 187.2 / 209.8 / 237.2. Screen kept CH-CH3 (scan
  1st), ethyl (scan 8th) and (CH3)2CH-NH3+ (scan 5th). Best: group proposal
  13CH3(125) + 13CH(203) + 13CH3(131) (bonded 13CH3-13CH, 13CH3-13CH3)
  + HX, fixed variant, accepted (-28.4): 1J 124.86 / 204.1 / 131.31 Hz,
  J(H1,H2) 10.5, J(C2,H3) 11.1 Hz; findings undetermined, misfit; c_hat
  290. The three change_proton_count proposals were rejected (+72 to +244).
  Free variants with a switched-off CH (0.04-0.18) are demoted. Motif
  CH-CH3 +129 / +145, ethyl +262 / +270. Nearest knowledge entries lactic
  acid / alanine at 25 Hz rms (no match). Reported as program output only.
- Regression rerun with D42 (runs/bench/search_real_v2.log): triethylamine
  best motif ethyl [free] (1J 130.98 / 125.0, 3J(HH) 7.12 Hz; before: CH-CH3
  + HX); ethylenediamine CH2-CH2 unchanged; pyridine ring [fixed]
  unchanged; lactic acid CH-CH3 + HX [fixed] (accepted -49; unchanged
  structure); alanine: the coupled-proton extension was rejected (-2.3,
  under accept_delta 6) but was reported as best because the ranking
  ignored the rejection. Fixed: rejected extensions rank after the other
  fits and cannot be best (D42 addendum); alanine's best is then its parent
  CH3-CH. The proton-count move was rejected on every sample except
  triethylamine's fixed variant, where ethyl -> CH-CH3 was accepted (-18.9);
  the overall best there remains ethyl [free].
- Blind check of D42 on 4322bdfc without the Et3N motif (added after the
  reveal; runs/blind/4322bdfc/auto_blind_v2.log, 15.7 min): the screen kept
  isopropyl (scan 4th), ethyl (scan 6th) and CH-CH3 (scan 1st); best =
  motif ethyl [free] (1J 130.98 / 125.00, 2J -4.39 / -3.06, 3J(H,H)
  7.12 Hz; amplitudes 1 : 0.94), next ethyl -> CH / CH3 variants (+25,
  +30, rejected) and isopropyl (+60, abundance). The program now finds the
  ethyl unit on its own; the tertiary-amine environment (through-N protons)
  still needs the Et3N motif or a remote-proton move.
- Confirmed by the user: N-ethylmethylamine (CH3-CH2-NH-CH3). The program's
  best had the right carbons except the CH2, read as a CH at 1J 204 Hz (its
  1.5 J band at 202.8 Hz); change_proton_count kept 1J at 204 when turning
  it into CH2 and was rejected (to fix: rescale 1J so the main X-Hn line
  stays put).
- Known-structure fits (complex route, fit phase): fast N-H exchange,
  fixed ratios and one rate: chi2 37342 (k 12), 115-121 s on one core,
  1J 124.89 / 134.93 / 131.71 Hz, 3J(H,H) 7.22, 3J(C,N,C,H) 3.88 / 3.86 Hz,
  phase0 -179 deg, delay -3.51 ms. Fixed ratios with free rates: 35769
  (k 14, 76 s). Slow exchange (NH kept) free: 22742 (k 27, 850 s),
  1J(15N,H) -65.3 Hz, but the CH2 carbon amplitude 2x natural; on the
  quasi-BIC scale the variants tie (241-244). Warm refits with the remote
  couplings freed (variants_fig.png): fixed ratios + free rates 23353
  (k 18), slow NH free 12408 (k 32; amplitudes far from natural). The fit
  stays well above noise (reduced chi2 90-270): the model still misses
  something (intermediate N-H exchange broadening is a candidate).
- Phase-first route: a model-free phase (hand-tuned minimum derivative
  entropy; delay -3.84 ms) was 35 deg off at the fitted optimum (residual
  zero-order phase -145.6 deg), and the real-only fit was worse (reduced
  chi2 745, 252 s). The entropy estimator is biased on this spectrum even
  noise-free (delay -3.30 vs -3.51 ms; flat profile) and failed the
  synthetic library test (correlation 0.56), so it was not added to the
  library (kept as runs/bench/phase_entropy_experimental.py). AsLS (standard
  and two-sided) did not separate baseline from lines; the broad humps are
  reproduced by the fitted model (processing response), so they must not
  be removed before fitting.
- Convention check (synthetic): the solver's shared gain phase and
  phase_delay give the absorption spectrum through phase_correct(values,
  f, phase0, delay, acquisition); cos FID -> absorptive, verified.

## Autonomous optimisation session (6 h; user: "optimise on your own, keep records")
- Scope (docs/PLAN.md Phase 3b A-H, R): 1J rescaling in change_proton_count,
  automatic variants (ratios, exchange, round-0 global starts), fit_structure
  entry, Voigt line shape, remote couplings, intermediate exchange, automatic
  phasing library function, automatic report, regression script.
  Budget: blind <= 30 min, known structure <= 10 min on 4 cores.
- 08:44 start. B: change_proton_count rescales 1J (main X-Hn line kept).
  A1-A4, H: ratios variant, thorough round-0 starts, fit_structure (exchange
  regimes x variants), report (table, J matrices, fit-phased figure). New
  moves: free_remote_couplings, gaussian_line_shape (model-level move);
  extensions refined in the parent's best variant only (budget).
- C (phasing): minimum-entropy and per-line estimators fail on dense J
  multiplets (noise-free rendering of the e66a4b08 fit: single-line phase
  errors 10-40 deg; real spectrum: delay 4.7 ms off). Phase-first route now
  uses the best complex fit's phase (fit_structure route="both");
  estimate_phase_lines kept for resolved spectra with its limits (D43).
- Found on the way (profiling one e66a4b08 refinement): 85 % of the time in
  NUFFT rendering over the full 65516-sample record although only samples
  200-4200 (+ SG half window) matter. Acquisition.local_record(): exact to
  4e-12, refinement 37.7 s -> 10.7 s (same evaluations).
- Motif scan: 1J of the best combination refined by coordinate search
  (+-4 Hz, 0.5 Hz): pyridine motif rank 2 -> 1; other confirmed samples
  unchanged; scan 5-27 s. The scan alone still ranks CH3-NH3+ first on
  alanine and lactic acid (the screen and refinement fix that downstream).
  Added motif / template / knowledge entry / benchmark case for
  N-ethylmethylamine (CH3CH2-N-CH3, fast exchange).
- R: configs/confirmed_samples.json + scripts/regression_confirmed.py.
- fit_structure on e66a4b08: v3 (8 starts, 240 s global search, all moves)
  29 min; found a warm-start bug (held 4J/5J couplings silently freed in every
  extension; free_remote_couplings never triggered): fixed, test added. v4
  (4 starts, 60 s search, structure-preserving moves only) 12.5 min: best
  fast exchange + remote J [fixed] chi2 23868 (k 16) vs fast fixed 35946
  (k 12); Gaussian width alone rejected (36266): the v3 "Gaussian" gain came
  from the freed remote couplings. Slow exchange fixed 37976 (v3 with 8
  starts: 25640): that model is start-sensitive; noted.
- Low-frequency lines: the fitted e66a4b08 model puts strong lines of the
  13CH2 isotopologue at 5-30 Hz (strongest 20.9 Hz, 0.54 of its maximum) and
  weaker ones at 30-31 Hz; the data show lines at 15-23 and 32 Hz. The
  standard ranges start at 62 Hz, so the worst-fitted isotopologue is fitted
  without its most informative lines. Test planned: ranges from 12 Hz.
- Residual after fast + remote J: mostly the 13CH2 isotopologue (194-216 Hz,
  line positions off by 1-2 Hz, e.g. 198.8 in the data vs 200.5 Hz), the
  131.8 Hz peak height, and a 235/237 Hz doublet the model shows as one line.
- Staged exchange (v5, standard ranges, 2 workers, 17 min): fast + remote J
  [fixed] best again (23868); slow exchange warm-started from the fast fit
  but with generic N-H couplings ended at 36308 (not better than fast) and
  slow + remote J at 54223 (worse than its parent: not a nested start).
  Changed: the slow model starts with its N-H couplings near zero (then it
  reproduces the fast fit at the start) plus a generic alternate start.
  Also added scripts/analyze_sample.py (one command for a new sample).
- Low-frequency ranges (v5low, ranges from 12 Hz): 32 min on 2 workers; the
  model follows 12-40 Hz roughly but not a 48-51 Hz feature (instrument or
  unmodelled), and the 13CH2 region is not improved. Not made a default.
- Nested starts (v6): slow exchange from the fast fit with N-H couplings near
  zero: 26921 (generic start: 36308); slow + remote J 23889 vs fast + remote J
  23868 with 9 more parameters (rejected). So slow N-H exchange does not
  explain the remaining misfit either. Gaussian width started at 0.3 Hz ended
  worse than its parent (not nested): now starts at 0.02 Hz.
- Added protonated forms (amine -> ammonium: one more H on each N with room)
  as warm-started slow-exchange extensions; e66a4b08 is likely protonated in
  water (pKa about 10.8). Run v7 pending.
