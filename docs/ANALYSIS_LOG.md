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
- Protonated form (v7, 39 min on 2 workers): slow exchange protonated
  (NH2+) [fixed] accepted over fast (22649, k 21); protonated + remote J
  reached chi2 17456 (k 25), 27 % below fast + remote J (23868, k 16), but
  on the quasi-BIC scale only 2.8 better than its parent (accept needs 6):
  a statistical tie, reported next to the best; chemically the ammonium
  form is expected in water. fit_structure budget trimmed afterwards
  (variants fixed + ratios, 3 starts, 45 s global search).
- Blind regression (2 workers): alanine and lactic acid skeletons found at
  rank 1 (7 min each). Running.
- 11:14 container restart killed both regressions (blind: 4 of 6 done, all
  skeletons correct at rank 1, 1.5-12 min each on 2 workers; known: alanine
  225 s, lactic acid 147 s). Resumed for the remaining samples.
- Found in the known-structure tables: (1) identical extensions refitted in
  the second round (fixed: skipped); (2) two-step improvements rejected step
  by step (lactic acid slow exchange + remote J was 7.4 better than the kept
  fit, both steps < 6): acceptance now compares with the nearest kept
  ancestor. Alanine: slow exchange (NH3+ kept) chi2 506 vs 852 fast, k 17
  vs 10, tie on the quasi-BIC scale (+0.3).
- Regression round 1 (code before the fixes below; 2 workers each):
  blind: alanine, lactic acid, pyridine (+ Gaussian width), ethylenediamine
  (+ Gaussian width) skeletons correct at rank 1; triethylamine wrong: the new
  CH3CH2-N-CH3 motif + remote J won (chi2 10074, k 18), ethyl screened out
  (6th of 8 in a fixed-variant screen), Et3N motif with a 1J of 157.8 Hz from
  the scan. Known: alanine 225 s, lactic acid 147 s, pyridine 1001 s (all
  fast exchange [ratios] best); ethylenediamine crashed (one component: no
  variant applied with variants fixed + ratios).
- Fixes: screen in the ratios variant (Et3N: ethyl 1st, was 6th); models with
  one component fitted as free; extensions compared with the nearest kept
  ancestor; identical extensions not refitted. Final regression (both modes,
  4 workers) started 11:47.
- Final regression (13:29, both modes, 4 workers; runs/regression/final):

  | sample | compound | known best | red. chi2 | time s | 1J | blind best | skeleton ok | time s |
  |---|---|---|---|---|---|---|---|---|
  | e3d282da | L-alanine | L-alanine (fast exchange) [ratios] | 17.76 | 122 | [129.87, 145.39] | 13CH3(129) + 13CH(146) (bonded: 13CH3-13CH) [ratios] | True (1) | 288 |
  | b683220d | lactic acid | lactic acid (slow exchange) + remote J [ratios] | 18.81 | 142 | [129.06, 146.79] | 13CH3(128) + 13CH(148) (bonded: 13CH3-13CH) + HX on X [ratios] | True (1) | 269 |
  | 7ad4aafb | pyridine | pyridine ring (fast exchange) + Gaussian width [fixed] | 70.59 | 1012 | [163.08, 163.24, 178.68] | motif pyridine ring [A2=173.5, A3=160.2, A4=161.5] + Gaussian width [fixed] | True (1) | 542 |
  | fde3fbb2 | ethylenediamine | H2N-CH2-CH2-NH2 (fast exchange) + Gaussian width [free] | 680.02 | 643 | [132.29] | motif CH2-CH2 [C1=132.5] + Gaussian width [free] | True (1) | 54 |
  | 4322bdfc | triethylamine | N-ethyl (Et3N) (fast exchange) + remote J [ratios] | 114.94 | 585 | [125.02, 130.92] | motif ethyl [C1=127.0, C2=123.2] [ratios] | True (1) | 325 |
  | e66a4b08 | N-ethylmethylamine | N-ethylmethylamine (fast exchange) + remote J [fixed] | 195.64 | 1351 | [124.9, 131.75, 135.63] | 13CH3(125) + 13CH(203) + 13CH3(131) (bonded: 13CH3-13CH, 13CH3-13CH3) + HX on X [ratios] | False (2) | 866 |

  Blind: 5 of 6 skeletons correct at rank 1 (triethylamine now ethyl, was
  wrong in round 1); N-ethylmethylamine correct at rank 2: the CH -> CH2
  move with rescaled 1J was accepted (-83) but an added coupled proton on the
  CH form scores better (20879 vs 31014, k 23 vs 17). Known structure:
  lactic acid now slow exchange (OH kept) + remote J, as in the hand
  analysis (H8); triethylamine reduced chi2 115 (was about 430 with fixed
  ratios); fitted delays -3.4 to -3.9 ms on all samples except
  ethylenediamine (+2.2 ms: one narrow band, delay and phase trade off) ->
  optional instrument delay prior added (processing.phase_delay_bounds_s).
  Budget: known-structure fits 2-4 min for two-carbon molecules, 10-22 min
  for pyridine and the 3-isotopologue amines (11-spin slow/protonated forms);
  blind 1-14 min.
- Two extension rounds in the blind search, e66a4b08 (27.7 min, concurrent
  with the full test suite): best = CH -> CH2 (1J rescaled) + remote J,
  correct skeleton (chosen over the +HX fit by the clean-margin rule; +HX
  scores 6.6 lower but carries warnings). Full test suite: 179 tests OK.
  Blind regression with --blind-rounds 2 on the other five samples started
  13:59; the default stays one round until it passes.
- Blind regression with two extension rounds (runs/regression/blind_rounds2):
  alanine 308 s, lactic acid 521 s, pyridine 575 s, ethylenediamine 53 s,
  triethylamine 323 s, all right skeleton at rank 1; with e66a4b08 (rank 1,
  1660 s alongside the test suite) 6 of 6. Two rounds made the blind default.
- 14:30 end of the 6 h session. Open items (PLAN): known-structure budget for
  the 11-spin exchange forms (10-22 min), the physics double diagonalisation
  (compute_transitions and transition_derivatives on the same Hamiltonian),
  instrument delay prior values, intermediate exchange, low-frequency range
  (a 48-51 Hz feature not modelled).

## Phase per dataset (user: phase correction is data processing and must be recalibrated for every dataset)
- Fitted phases of the six confirmed samples cluster: phase0 about -174 deg,
  delay -3.62 ms (sd 0.17 ms). A leave-one-out mean calibration is within
  12 deg of each sample's own complex-fit phase at its line frequencies; the
  per-spectrum model-free estimators (minimum entropy, per-line fits) in a
  narrow window around it were worse (15-89 deg): not usable per dataset.
- The raw FID (samples before the 0.05 s crop) shows the field switch-off:
  a plateau to 2.75 ms, a falling edge (half height 3.41-3.51 ms in every
  dataset), then ringing to about 25 ms. The pulse sequence in the uploaded
  NMRduino settings is standard_zf_4000Hz_no_dead.seq. The edge time equals
  the fitted delay within 0.05 ms for alanine, triethylamine and
  N-ethylmethylamine; lactic acid and pyridine fits differ by 0.34 / 0.44 ms.
  So the delay can be measured per dataset from the raw FID (no model).
- Zero-order phase from the J-spectrum alone, with the delay fixed at the
  edge, is still off by 12-80 deg on several samples (lines of both signs
  overlap): a per-dataset phase0 needs another handle. Test running: lactic
  acid and pyridine refitted with the delay fixed at the edge time.
- Edge vs fitted delay: lactic acid refitted with the delay fixed at the edge
  (-3.43 ms): chi2 1184.6 vs 1151.9 free (-3.78 ms), 1J unchanged, phase0
  shifted by 18 deg (compensation) -> the edge value is consistent with the
  data. Pyridine: the free delay went to +1.0 ms (chi2 5881; the regression
  run found -3.9 ms), edge-fixed 7326 with 1J(A2) 178.7 vs 176.4 Hz: the free
  delay absorbs model error in a narrow band. Added
  zulf_core.render.switching_edge_delay (per-dataset delay from the raw FID,
  no model; test with a synthetic edge). Per-dataset zero-order phase remains
  open (the J-spectrum alone does not fix it). Also fixed: a delay prior that
  excludes zero crashed the parameterisation (start now inside the bounds).

## zulf_processing package (user: global search then fine-tune; SG and crop per spectrum; processing as its own part)
- New package zulf_processing (raw, diagnostics, plan, phase, dataset; D44).
  Plan per dataset: crop start = max(default, end of this dataset's ringing)
  (e66a4b08: 203 instead of 200), stop moved to keep 1 s; SG / apodization
  still defaults (next step). Phase: grid over delay (+-0.5 ms around the
  switching edge, 10 us) x phase0 (0-180 deg, 1 deg), 4 best local minima,
  Nelder-Mead fine-tune, sign by the strongest point.
- Speed: process_dataset took 546 s per dataset, all in the 1-3 exponential
  curve fits of zulf_core.diagnose_fid over 64k points (unused by processing).
  Skipped in diagnose_raw -> 0.1 s. The first validation runs had timed out
  for this reason. Test added (diagnostics < 20 s).
- Validation (scripts/validate_processing.py, runs/processing_validation):
  phase error at 125 / 200 / 250 Hz vs each sample's complex fit (deg):
  calibration leave-one-out (median delay, ethylenediamine excluded from the
  phase mean): alanine 8/8/7, lactic -14/-7/-2, pyridine -12/-2/5,
  ethylenediamine -23/3/80 (its reference fit sits at +2.2 ms, ambiguous),
  triethylamine 2/2/2, N-ethylmethylamine -1/-1/-1.
  entropy:edge_fixed: -80/-79/-78, 33/42/48, 16/28/36, -9/17/-87, -2/-2/-1,
  13/12/12. lines:edge_fixed: -53/-51/-50, 26/35/41, -18/-6/2, 84/-70/7,
  42/43/43, -22/-22/-23. edge_prior delays run 0.25-0.5 ms past the edge
  (to the window bound): the criteria trade delay against phase0.
  The earlier "leave-one-out mean" row (52-90 deg errors) was an artefact:
  the ethylenediamine +2.2 ms delay pulled the mean delay by 1 ms.
- Decision: default phase = instrument calibration + this dataset's edge
  (calibrated_phase; config phase_calibration: phase0 176.3 deg, delay =
  edge - 0.033 ms). calibration:edge errors: 7/7/8, -10/-2/4, -12/-1/6,
  -25/0/76, 1/1/0, -3/-5/-5 (not independent: fitted on all samples).
  analyze_sample uses it (--phase entropy|lines to compare). Complex fits
  are unaffected (they fit their own phase); the processed phase is for the
  phased route and display.

## Window length per dataset (user: a too-short window loses spectral information)
- The records are 16.4 s; the fixed processing used 1 s (samples 200-4200).
  Signal extent per dataset (zulf_processing.signal_extent: power inside the
  signal bands per 0.25 s block vs the same bands in the last 30 % of the
  record): energy inside 1 s: alanine 88 %, lactic acid 88 %, pyridine
  99.6 %, ethylenediamine 98 %, triethylamine 89 %, N-ethylmethylamine 86 %;
  signal end (power stays below 3x floor for 1 s): 1.25, 2.75, 1.75, 3.5,
  4.75, 4.5 s. The missing 11-14 % are the narrow lines, and frequency
  information grows with t^2, so the tail matters for small couplings.
- For a forward-model fit a longer window does not lose information (the
  model applies the same processing; a noise-only tail adds no bias); costs
  are compute (record length) and exposure to late drift. Zero filling adds
  no information and makes the fitted points correlated, so with a longer
  window zero_fill is reduced to keep 4 points per Hz (same number of fitted
  points).
- Implemented: diagnostics.signal_extent (test against the closed form for
  decaying cosines in white noise: end 5.75 / 1.50 s vs 5.76 / 1.50 s), plan
  window_mode "signal_extent" (window = max(default, end), cap max_window_s
  8 s, zero_fill = round(4 / window_s)); default stays "fixed" until fits
  confirm. regression_confirmed.py --processing dataset --window-mode.
  Resulting windows: 2.0, 3.0, 1.75, 3.5, 5.25, 4.5 s.
- Running: known-structure fits of lactic acid, triethylamine and
  N-ethylmethylamine with the signal-extent windows
  (runs/regression/window_extent), to compare with runs/regression/final.
- Results (signal-extent window vs the 1 s window of runs/regression/final):
  lactic acid (3.0 s, zero_fill 1): same best model, 1J 129.15 / 147.04
  (was 129.06 / 146.79), delay -3.528 ms (was -3.746; edge -3.407),
  phase0 -178.9 deg (was -159.7; calibration 176.3), phased-route reduced
  chi2 34.6 (was 85.5), 803 s (was 142 s).
  Triethylamine (5.25 s): best model now slow exchange, protonated
  (ammonium) + remote J (was fast exchange), 1J 124.98 / 131.02 (was
  125.02 / 130.92), delay -3.459 ms (edge -3.472), phased-route chi2 578
  (was 1138), 2195 s (was 585 s).
  Reading: 1J stable within 0.25 Hz; delay and phase0 much better fixed
  (closer to the switching edge and to the calibration, no delay-phase
  trade-off); exchange regime changes with the late narrow lines. Complex
  reduced chi2 not comparable (zero_fill 4 -> 1: uncorrelated points).
  Cost grows with the record length (render + SG over the whole record):
  3.8-5.6x.
  N-ethylmethylamine (4.5 s): same best model (fast exchange + remote J,
  fixed variant), 1J 124.92 / 131.75 / 136.98 (was 124.90 / 131.75 /
  135.63: the third 1J moved 1.35 Hz), delay -3.363 ms (was -3.495; edge
  -3.510: now further off), phase0 168.1 deg (was 177.4), phased-route chi2
  589 (was 881) but its residual phase is 29.7 deg (not self-consistent),
  5804 s (was 1351 s). Mixed: 25 parameters over a 4.5 s record; the delay
  moved away from the edge. Default window stays "fixed" for now. Next:
  renderer cost for long records, then a soft (flat + Gaussian edge) crop
  window and weighted SG (user proposal), then re-run all six.

## Blake pyridine series (processed spectra, real part; x_pyridine 0.02-1.00, pyridinium in 6 M HCl)
- Data: frequency / amplitude arrays, 0-4000 Hz, 88458 points (spacing
  4000/88457 Hz: a plain FFT of 88457 samples at 4 kHz, 22.1 s). Given as
  phased, baseline-corrected real parts (user). Fitted range 140-200 Hz.
  Tool: scripts/fit_processed_spectrum.py (known structure, coupling
  overrides; --record 4000,88457 renders the finite record so every line
  shape has analytic derivatives: with record=None the Gaussian-width
  extensions fell back to finite differences and a fit took 65 min).
- Start: pyridine motif with literature couplings (J23 4.9, J34 7.7,
  J24 1.8, J35 1.4, J25 0.9, J26 -0.1; 1J 178 / 162.5 / 161.5, long-range
  CH from tables); predicted lines already close to the data.
- First fit (x = 0.50, exchange auto): best "slow exchange, protonated +
  remote J", k 53, reduced chi2 21.5; neutral symmetric pyridine k 26,
  chi2 3.8x higher. Cause: `protonated()` dropped every symmetry generator,
  so C2/C6 and C3/C5 became separate isotopomers with free couplings
  (double the parameters), which fitted line-shape detail; the N-H
  couplings themselves were below 1.1 Hz. Not physical for neutral pyridine
  in water. Fix: protonated() keeps generators that map the new N-H onto an
  existing label (pyridinium keeps the mirror; test added). This also
  affected the protonated forms in the confirmed regression (triethylamine,
  N-ethylmethylamine, ethylenediamine): re-run the regression.
- Re-run: pyridine samples with exchange "fast" (given protonation state
  only, symmetric), pyridinium as the protonated symmetric fragment
  (exchange auto: N-H decoupled vs coupled). C3/C4 1J of pyridinium from a
  scan: 171 Hz; C2 ambiguous (179 or 190 Hz), both starts to be fitted.
- Speed: with record=None the Lorentzian route has analytic derivatives
  (Jacobian 0.66 s); with --record 4000,88457 every direction needs an FFT
  of 88457 = 53 x 1669 points (5.3 s). Gaussian-width extension dropped for
  this series (it did not win on x = 0.50). A single known-structure fit:
  about 10 min.
- Baseline (symmetric neutral pyridine, 4 starts, spread 0.5 Hz): reduced
  chi2 4.7 (x 1.00) to 105 (x 0.50); 1J consistent (C2 178-179, C3 161-163,
  C4 162-164 Hz) but H-H and long-range C-H couplings jump between samples
  (J(H2,H4) 5.7 vs 3.1 Hz; 3J(C2,H6) 13.9 vs 8.7 Hz): local minima.
  Pyridinium (protonated, symmetric; best: N-H decoupled): 1J 197 / 172 /
  175 Hz, other couplings implausible (J(H3,H4) 10.3 Hz, bounds hit).
- Checks on x = 1.00: (1) delay held at 0: chi2 1.75x worse, so the delay
  stays a processing parameter; (2) the simulated absorption spectrum at
  accepted couplings has no negative lines, and correlates better with the
  data than the magnitude (0.56 vs 0.45): the data are real parts, as
  given; (3) user: do not hold couplings, they are the result. Residual on
  the data cores: couplings at the starting (accepted H-H, unverified
  long-range C-H) values 0.47; H-H held, C-H free 0.25; all free from the
  accepted values 0.214 (J(H3,H4) 3.7 Hz); staged (H-H held first, then all
  free) 0.210 with J(H3,H4) 7.24 Hz (accepted 7.7): the direct all-free fit
  ended in a worse minimum. Multi-start (24 starts, spread 2 Hz) recorded
  only 3 starts (evaluation cap 20000): relative scores 1.00 / 1.31 / 2.08
  with very different couplings.
- Added: zulf_hypothesis.uncertainty (linearised standard errors and
  correlations by variable projection; Monte Carlo test), every start's
  solution kept in the refinement summary (start_solutions, both
  signal-mask passes), scripts/fit_staged.py. Running: staged fits of all
  seven spectra (12 stage-2 starts, spread 1 Hz), runs/processed/*_staged.
- Staged fits of all seven (12 stage-2 starts, spread 1 Hz; 36-72 min each,
  4 in parallel): relative residual on the data cores 0.171 (x 1.00),
  0.186, 0.206, 0.254, 0.236, 0.260 (x 0.02), pyridinium 0.195; lower than
  every earlier fit. But in every spectrum the best minimum was reached by
  one start of twelve (next ones at 1.03-1.41x the score), the fitted delay
  is 3.7-6.4 ms, and the couplings other than 1J scatter by several Hz
  between neighbouring concentrations (J(H3,H4) 6.5-10.2 Hz, 3J(C2,H6)
  6.8-16.3 Hz) while their linearised errors are 0.02-0.4 Hz. 1J also
  scatters (C2 176.7-179.9, C4 157.8-166.3 Hz). Conclusion: with one
  spectrum at a time these data do not fix the small couplings (many
  near-equivalent minima); the linearised errors describe one minimum, not
  the ambiguity. Pyridinium: 1J 187.5 / 169.0 / 171.2 Hz.
  Next: joint fit of the six pyridine spectra (couplings shared or
  constrained across concentration), or the raw FIDs (complex data).
  Overview figure: runs/processed/overview_staged.png.
- User: fit the small couplings around the literature values with a weight.
  Added Gaussian priors to refine (RefineSettings.priors, prior_weight:
  residual rows sqrt(w) (x - mean) / sigma / norm, so w = 1 prices a
  one-sigma deviation like one data point one noise sigma off; test against
  the linearised closed form). Running: staged fits with priors centred on
  the starting couplings, sigma 0.5 Hz (H-H, accepted values) and 1.5 Hz
  (long-range C-H, starting values not verified against the literature:
  Shiner and Wyllie 1973 not reachable from here), weight 10 (the fits'
  reduced chi2 is 5-10, so a noise-sigma prior would be too weak);
  runs/processed/*_prior.
- Single-spectrum fits with priors (sigma 0.5 / 1.5 Hz, weight 10): residual
  x 1.00 0.155, x 0.50 0.227, x 0.02 0.227 (all lower than without priors),
  pyridinium 0.207. User: couplings should change monotonically with the
  concentration.
- Joint fit of the six pyridine spectra (scripts/fit_joint_series.py): every
  coupling J(x) = a + b (x - mean x), per spectrum its rates, delay, gains,
  phase; priors on a as above; 4 starts, 11 min. Scores 0.271 / 0.279 /
  0.282 / 0.656. Residuals 0.225, 0.231, 0.219, 0.177, 0.176, 0.166 (x 0.02
  ... 1.00): as good as or better than the single-spectrum fits for five of
  six spectra with two parameters per coupling, so a monotonic (linear)
  concentration dependence is consistent with the data.
  Trends (J at x 0.02 -> 1.00, Hz): 1J(C2,H2) 178.3 -> 177.7, 1J(C3,H3)
  164.5 -> 162.8, 1J(C4,H4) 161.9 -> 158.2; J(H2,H3) 3.7 -> 4.9 (neat
  literature 4.9), J(H3,H4) 7.9 -> 7.8 (7.7), J(H2,H4) 2.6 -> 3.5 (1.8),
  J(H2,H5) 2.4 -> 1.8 (0.9), J(H3,H5) -0.8 -> 3.2 (1.4); large slopes for
  J(C2,H3) (-9.7 Hz per mole fraction), J(C3,H5) (-5.4), J(H3,H5) (+4.1):
  suspicious, likely trading off against each other. Output:
  runs/processed/joint_linear.
- User: not linear, only monotonic. Joint monotone fit (every coupling free at
  each concentration, one-direction steps; direction from the linear fit's
  slope; priors on the series average; start from the linear fit; 6 starts,
  15 min): scores 0.218 / 0.224 / 0.225 / 0.226 / 0.228 / 0.228 (starts
  within 4 %). Residuals x 0.02 ... 1.00: 0.181, 0.192, 0.177, 0.162, 0.165,
  0.131, the lowest of every fit so far. Several couplings change in steps
  (flat segments where the step is 0). Plausible: 1J all decrease with
  pyridine fraction (C2 178.6 -> 177.8, C3 164.2 -> 162.8, C4 161.1 ->
  157.8), J(H3,H4) 8.7 -> 8.1, J(H2,H3) 3.3 -> 5.3. Still large changes:
  J(C2,H3) 7.0 -> -2.4 Hz, J(C3,H5) 10.4 -> 5.7, J(H3,H5) -1.3 -> 3.2,
  J(C4,H3) and J(C4,H2) at x 1.00 with errors > 1 Hz. Table:
  runs/processed/joint_monotone/J_table.csv.
- User: never linear. fit_joint_series rewritten: J_k(x_i) = v_k + A_k c_k,i
  with c rising 0 -> 1 through softmax step shares (direction = sign of A,
  shape free); no linear stage anywhere. From the previous monotone result
  (A): 6 starts, 5 at 0.2181, residuals 0.178, 0.193, 0.178, 0.162, 0.165,
  0.131 (same solution as before). From the 7ad4aafb couplings (user,
  symmetrised; B): 0.316-0.519, worse minima.
- Weights: the signal weighting gave weak peaks (148-152, 175-178 Hz) the
  outside weight 0.2 at the 4 sigma core threshold. Figures
  runs/processed/weights_threshold*.png, weights_taper.png. User choice:
  2.5 sigma threshold and a slower fall-off between peaks: taper 4 Hz
  (lowest weight between peaks 0.60-0.69, noise regions 0.27-0.31).
  Running: runs/processed/joint_mono_thr25_taper4 (from A).
- 2.5 sigma / 4 Hz weighting (from A): scores 0.3265, then 0.3335-0.3377 (best
  reached once). Criterion (user): the residual is not the measure; judge the
  couplings by reproducibility, robustness to processing choices, errors and
  correlations, prediction, plausibility. Weighting comparison (max |J diff|
  over x between the two weightings; formal errors 0.01-0.16 Hz):
  < 0.3 Hz: 1J(C2,H2) 0.18, 1J(C3,H3) 0.11, J(H3,H4) 0.17;
  0.5-1 Hz: 1J(C4,H4), J(H2,H3), J(H2,H4), J(H2,H5), J(H2,H6), J(C2,H4),
  J(C2,H5), J(C3,H2), J(C3,H4), J(C3,H6), J(C4,H3);
  > 1 Hz (not determined): J(C2,H3) 2.5, J(C4,H2) 2.5, J(H3,H5) 1.3,
  J(C2,H6) 1.25, J(C3,H5) 1.2. Directions agree for all 19 couplings. The
  weighting sensitivity is 5-20x the linearised errors: those understate
  the uncertainty. Figure runs/processed/J_weighting_comparison.png.
- Robustness under the 2.5 sigma / 4 Hz weighting: 16 starts (spread 1 Hz):
  best 0.318, four starts within 3 %; leave-one-out (5-spectrum joint fits,
  4 starts each; bracketed prediction of the left-out spectrum). Prediction
  residuals x 0.33 / 0.50 / 0.66 / 0.75: 0.242 / 0.240 / 0.279 / 0.251 vs
  0.240 / 0.224 / 0.247 / 0.239 inside the full fit: no gross over-fitting
  of the concentration trend. Per coupling (max over x): spread among the
  near-best starts, jackknife error, weighting difference, linearised error
  (runs/processed/robustness_summary.json). Reliable on every measure
  (< 0.3 Hz): 1J(C2,H2) 178.5 -> 177.6 Hz and 1J(C3,H3) 164.1 -> 162.8 Hz.
  Moderate: J(H2,H5). All other small couplings (and 1J(C4,H4)) differ by
  1-4.6 Hz between near-equal starts: the processed real spectra do not fix
  them; the linearised errors (0.01-0.16 Hz) are not meaningful here.
- Joint fits with the NMRduino pyridine FID 7ad4aafb (50 % in water; user
  unsure whether by volume or by mole: both tried; complex data, confirmed-
  sample processing, 140-179 and 181-200 Hz; 8 starts each):
  by volume (x 0.183, own node between 0.02 and 0.33): best 0.363,
  residual 7ad4aafb 0.206, Blake 0.18-0.26 (Blake x 0.50: 0.239);
  by mole (x 0.50, sharing the Blake x 0.50 node): best 0.433, residual
  7ad4aafb 0.305 and Blake x 0.50 0.289 (was 0.224 without 7ad4aafb): the
  two spectra do not fit one set of couplings at x 0.50, so the volume
  reading is the consistent one (or the instruments differ in line shape).
  The small couplings moved by several Hz against the Blake-only fit (e.g.
  J(H3,H4) at x 0.02 13.2 vs 9.5 Hz, J(H2,H4) 6.1 vs 2.8 Hz) with 8 starts:
  adding the FID did not pin them down; 1J(C2,H2) and 1J(C3,H3) stayed
  within 0.1-0.8 Hz.
- 48-start joint fit (24 perturbed around the 16-start best, 24 drawn from
  the priors; 4 worker processes, 28 min): best 0.3084 (16 starts: 0.3181),
  8 starts within 3 %, 1 within 1 %. Spread of the near-best (3 %)
  solutions, max over x: 1J(C3,H3) 0.14 Hz, 1J(C2,H2) 0.61 Hz, all other
  couplings 1.2-6.0 Hz (1J(C4,H4) 3.4 Hz). More starts lower the best score
  slightly but show more near-equivalent minima: the small couplings are not
  unique from these spectra. Figure runs/processed/joint_new_ms48/
  near_best_couplings.png. Next: synthetic recovery test with this best
  solution as the truth (information limit vs model mismatch).
- Synthetic recovery test (scripts/synthetic_recovery_series.py): the 48-start
  best solution as truth, rendered on the measured grids with the measured
  noise (robust sigma 19-136), fitted like the data (from the accepted
  values, 16 perturbed + 16 prior-drawn starts). Best found 0.497 while the
  truth scores 0.353 (0.318 without priors): the local fits never reached
  the truth basin. Errors of the best found: 1J(C2,H2) 0.37, 1J(C3,H3)
  0.52, 1J(C4,H4) 3.6, small couplings 0.7-6.7 Hz. So at least part of the
  non-uniqueness on the real data is a search failure, not only missing
  information. Added a coarse-to-fine continuation (Gaussian smoothing of
  data and model with the same kernel, 1.5 -> 0.8 -> 0.4 Hz -> none; test
  against finite differences); running it on the synthetic series.
- Smoothing continuation on the synthetic series (16 starts): best 0.452
  (plain 0.497, truth 0.353) but median coupling error 4.3 Hz (plain 1.8):
  lower scores without getting closer to the truth. Basin test from the
  truth: spread 0.3 Hz: 8/8 starts reach 0.3434 (below the truth's 0.353,
  noise fitted) with max coupling error 0.67 Hz over all couplings and
  concentrations; spread 1 Hz: 2/8. So the spectra carry the information
  (errors < 0.7 Hz in the right basin), but the basin is about 1 Hz wide in
  19 coupling levels: random starts almost never land in it. Added a
  coordinate grid scan (each small coupling level on a +-4 Hz grid, others
  held, then a local fit; test: a 3 Hz displaced coupling is recovered);
  running on the synthetic series from the accepted values.
- Coordinate scan on the synthetic series (8 starts, 2 cycles): best 0.637,
  worse than plain multi-start (0.497): greedy moves along one coupling at a
  time drift into poorer regions. Not used. Next: differential evolution per
  spectrum (scripts/global_single_spectrum.py; couplings within the starting
  value +- 5 Hz, rates and delay held, gains solved inside; then a local fit),
  tested on synthetic x 1.00 and 0.02.
- Differential evolution (+-5 Hz, popsize 10, 60 generations, 11590
  evaluations, 7.5 min on 2 workers) on synthetic x 1.00 / 0.02: scores
  0.107 / 0.120 vs the truth's 0.052 / 0.074; median error 1.9 Hz. The
  synthetic truth (the 48-start best on the real data) has couplings more
  than 5 Hz from the starting values (x 1.00: J(C2,H3) -2.5 vs 3.0, J(C2,H6)
  17.2 vs 11.2), outside the search box; repeated with +-8 Hz and 100
  generations.
- DE with +-8 Hz, 100 generations (19190 evaluations, 12 min): x 1.00 score
  0.104 (truth 0.052), median error 0.93 Hz but J(C2,H3) and J(C2,H6) still
  6-7 Hz off; x 0.02 0.126 (truth 0.074), median 1.7 Hz. Not a fix.
  Caveat of these tests: their truth is the real-data best solution, which
  itself has implausible values (J(C2,H6) 17 Hz vs about 11 Hz accepted).
  Fairer test running: a plausible truth (small couplings = accepted values
  + N(0, 1 Hz) offsets with small monotone trends, 1J with the fitted trends),
  fitted from the accepted values with the standard 24-start joint fit
  (runs/processed/synthetic_plausible*).
- Decisive identifiability test (plausible truth): small couplings = accepted
  values + N(0, 1 Hz) offsets with small trends, 1J with the fitted trends,
  measured noise; standard 24-start joint fit from the accepted values.
  Best found 0.651 (8 starts within 0.651-0.659) with coupling errors up to
  6.4 Hz (median 1.2-2.5 Hz; 1J(C4,H4) trend wrong); the truth scores 0.666
  with the priors and 0.6473 data-only vs 0.6460 for the best found. So
  couplings several Hz apart fit such spectra equally well at this noise
  level: an information limit, not a search failure. (The earlier synthetic
  truth, the real-data best solution, was a noise-fitted special case.)
  Conclusion for the Blake series: only 1J(C2,H2) and 1J(C3,H3) (and their
  decrease with the pyridine fraction) are determined; the small couplings
  (H-H and long-range C-H) cannot be determined from these processed
  spectra, even with a correct global search. Needed for them: more
  information per spectrum (raw FIDs with known processing, longer records,
  higher SNR, other isotopologues / fields) or reliable literature values.
- CORRECTION of the entry above: with --hold-small-first (stage 1 of every
  start: small couplings held at the accepted values, 1J trends, rates and
  delays fitted; stage 2 all free) the same synthetic series gives best
  0.6418 (below the truth's 0.666 and the plain fit's 0.651) with coupling
  errors median 0.72 Hz, max 1.46 Hz, 1J errors 0.17-0.60 Hz (1 of 8 starts;
  the others 0.66). So the near-truth minimum is the best one and reachable;
  the plain fits failed because constant 1J starts (the 1J change with the
  concentration, up to 4.5 Hz for C4) pushed the small couplings into
  compensating minima. The margin between minima is small (1-3 % of the
  score), hence many starts are still needed. Not an information limit.
  Running: the real series with --hold-small-first from the accepted values.
- Real Blake series with --hold-small-first from the accepted values (24
  starts, 21 min): best 0.3358, next 0.364-0.379 (1 start within 3 %); the
  48-start fit (started from earlier solutions) had reached 0.3084. The small
  couplings of the staged best differ from the 48-start best by several Hz
  (e.g. J(H3,H4) 7.6 -> 4.5 vs 8.6 -> 7.9 Hz). On the synthetic series the
  staged fit reached the truth basin; on the real data it does not reach the
  best known score: the real spectra contain something the model does not
  (line shape of the unknown processing, 14N effects on the C2 lines,
  baseline), and the small couplings absorb it. 1J(C2,H2) and 1J(C3,H3)
  agree between all fits (178.5-179.5 -> 177.3-177.9 Hz, 164.1-164.5 ->
  162.5-162.8 Hz). End of the 5 h session; open items in PLAN Phase 3d.
- Literature read: Wilzewski, Afach, Blanchard, Budker, J. Magn. Reson. 284
  (2017) 66-72 (density-matrix fit of ZULF spectra, 13C/15N-labelled neat
  liquids, sub-mHz precision). Points relevant to the Blake series:
  (1) they fit the magnitude spectrum because finite pulses leave small
  phase errors; (2) one exponential decay for all coherences; (3) methyl
  formate showed non-Lorentzian lines that a fitted residual field
  (transverse 3.1(3) nT, longitudinal 2.8(4) nT) explained, reduced chi^2
  1.96 -> 0.98; (4) the magnetometer's frequency-dependent amplitude and
  phase response biased the residuals of distant multiplets in opposite
  directions (calibration measured 4-400 Hz); (5) 2 Hz cut around every
  mains harmonic; (6) uncertainties chi^2/(n-p) C^-1 (as ours) with
  reparametrisation (mean, difference) of strongly correlated pairs.
  Checked: the Blake spectra have no narrow line at exactly 150 or 180 Hz
  (maxima within +-1 Hz are molecular lines, at 149.6-151.0 and 180.4-181.0
  Hz). Our protocol accepts a fixed static field (Protocol.field_ut) but the
  solver does not fit it; a fitted residual field and a smooth frequency-
  dependent gain are candidate model extensions for the real-data mismatch
  of the staged fit (PLAN Phase 3d). Their precision comes from labelled
  neat samples with high SNR; natural-abundance mixtures are not comparable.
