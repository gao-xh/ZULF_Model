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
