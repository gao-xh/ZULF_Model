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
