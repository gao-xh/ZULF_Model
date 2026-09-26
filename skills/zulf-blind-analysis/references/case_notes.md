# Case notes (blind tests, 2026-09-26)

Experimental files are not in the repository; paths refer to the uploads of
that session. Numbers are relative residuals on the stated scale.

## 7dc9a043 (no molecular signal)
- No band above noise in 20-520 Hz with any window; only instrument lines and
  a slow 1.6-2.2 Hz decaying oscillation also present in other datasets.
- Correct outcome: report "no J-coupled signal"; do not fit.
- Useful afterwards as the no-signal reference for background checks.

## 7ad4aafb (aromatic, pyridine-consistent)
- One band 148-192 Hz, rich fine structure, no 2J band: only CH carbons.
- Pyridine 13C2 + 13C3 + 13C4 with literature-like starts: 0.144 (complex,
  full record); benzene-13C1: 0.579. Without symmetry ties the refined
  couplings broke the ring symmetry: tie symmetry partners.
- Models: skeleton 13C + 5 H, 1J 163.5 Hz, small couplings unresolved.
- Blind global search from that skeleton (3 components, wide bounds): 0.212
  with unphysical couplings (1J 159/167/168, long-range +-15 Hz): search did
  not find the basin.

## e3d282da (weak; CH-CH3 unit) - CONFIRMED: L-alanine
- Lines 129.2/130.5, 144.5/146.0/148.5, 258.8 Hz (SNR about 4-15).
- Full record fits were insensitive (no-lines 0.856, floor 0.72); the
  0.05-1.05 s window separated hypotheses (no-lines 0.711, floor 0.37).
- CH-CH3 fragment, both 13C isotopologues, shared 3J(HH): 0.358, all
  couplings physical (1J 145.4 / 129.8, 2J -4.1 / -4.5, 3J 6.6 Hz).
  Two isolated methyls (methyl acetate type) 0.617; isopropyl rejected by
  abundance; 15N (slow NH exchange) not supported after the background check.
- Candidates: lactic acid or alanine with fast NH exchange; the data cannot
  separate them. Decisive data: different pH (slow exchange reveals 15N-H
  near 73 Hz for alanine) or higher SNR for the carboxyl isotopologue.
- Models: no usable structure at this SNR (best refined model candidate 0.563
  on the window).

Confirmation (given by the user after the analysis): the sample was pure
L-alanine dissolved in water (single compound, no mixture). Our final pair was lactic acid or alanine with fast NH exchange,
so the answer was retained. Checks against the answer:
- 13C-alpha (CH, bonded to N) 1J 145.4 Hz and 13C-beta (CH3) 1J 129.8 Hz, one
  carbon each: consistent (amplitudes 1 : 0.77, expected 1 : 1).
- 3J(H-alpha, H-beta) fitted 6.6 Hz (complex) and 7.7 Hz (phased): the
  literature value from memory (about 7.2-7.3 Hz) lies between them.
- NH coupling: fast NH3+ exchange near pH 6 was assumed. The 15N
  improvement (1J(N,H) -73.3 Hz, amplitude 0.42 vs 0.35 expected) was
  withdrawn on a mean-level comparison with the no-signal dataset, which is
  not decisive: no 13C isotopologue has lines at 55-100 Hz, and intermediate
  exchange would broaden rather than remove the 15N lines. Unresolved;
  decide with a water blank, an acidic sample or D2O.
- The carboxyl 13C and 15N isotopologues only have low-frequency lines,
  outside the analysed range; they would be the route to separating alanine
  from lactic acid (15N exists only in alanine), together with a low-pH
  measurement that slows NH3+ exchange.
- Chirality (L vs D) is not observable in a J-spectrum.
- Pure sample: isotopologue abundances are fixed at 13C-alpha : 13C-beta :
  13C' : 15N = 1 : 1 : 1 : about 0.34. The free fit gave 1 : 0.77; fitting
  with the ratio fixed is the first test for abundance-constrained amplitudes.
- Ask for the sample state early (solvent, pH, pure or mixture): it decides
  the exchange model and the abundance constraints.

## b683220d (confirmed: lactic acid)

- Lines (0.05-1.05 s): 128.4 / 130.0 / 131.8, 144.3-151.9 (strongest
  147.75, SNR 17), 257-263 Hz; nothing at 60-100 Hz besides instrument lines.
- H7 (CH-CH3) free: 0.280, amplitudes 1 : 1.04, 1J 147.12 / 129.19,
  2J -4.38 / -4.07, 3J(HH) 6.52 Hz; methyl carbon broadened (4.7 1/s).
- H8 (H7 plus one weakly coupled proton on X): 0.190, 1 : 0.87, 1J 147.07 /
  129.18, 2J -4.69 / -4.40, 3J(HH) 7.06 Hz; extra-spin couplings not
  determined. Isopropyl 0.352 with 1 : 1.16 (needs 1 : 2): rejected.
- Final guess lactic acid / lactate, runner-up alanine at another pH.
  Confirmed by the user: lactic acid. Solvent and pH were not given; the
  extra spin (OH in slow exchange, or another small effect) stays open.
