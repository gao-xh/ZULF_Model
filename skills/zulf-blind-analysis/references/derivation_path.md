# Derivation paths (observation -> inference -> test -> outcome)

These are the actual reasoning chains of the blind tests, including the steps
that turned out wrong. Follow the same structure for a new sample: every
inference names the observation it rests on and the test that could refute it.

## Sample e3d282da (weak signal; final: CH-CH3 unit; confirmed L-alanine)

1. Observation: FID of 65516 points, no ini; plateau to 2.75 ms, ringing to
   45 ms, late noise RMS 0.92.
   Inference: same NMRduino as earlier data, 4000 Hz; crop at 0.05 s.
2. Observation (crop 0.05 s, SG 201, full record): weak lines at 130.5 Hz
   (8x noise), 145-148 Hz cluster (10x), 260 Hz (3x, only in shorter windows).
   Instrument lines 60/120/180/294/323/443/923 Hz identified by comparison
   with other datasets.
3. Inference from patterns: J and 2J (130.5 / 260) suggest a 13CH3 with
   1J about 130 Hz; 145-148 Hz is a second group.
   Test: 2J partner of the 146 Hz group near 292-294 Hz.
   Outcome: absent (and the 294 instrument line is weak here) -> the 146 Hz
   group is not an isolated methyl.
4. Models (4 checkpoints, 2 phasings, 82 distinct candidates): no couplings
   >= 1 Hz from the transformers; unphysical ones from the set model.
   Refinement of all candidates: 0.834-0.856 on the full record (no-lines
   0.856). Outcome: models uninformative at this SNR.
5. Hypotheses H1-H4 (full record): best H3/H2 about 0.76; fitted lines too
   broad; H2's methyl predicted a 2J line stronger than J, opposite to the
   data. Inference: the 130 Hz group may not be an isolated methyl either;
   both groups show about 1.4 Hz fine structure -> small couplings to other
   protons.
6. User question: why not the 0.05-1.05 s window?
   Test: recompute scale and refit. Outcome: range doubled (no-lines 0.711,
   floor 0.37); ranking sharpened: H6 0.354, H7 (CH-CH3 fragment, both 13C
   isotopologues, shared 3J(HH)) 0.358, H5 0.381, H2 0.483, H1 0.617.
7. Check of H6 vs H7: H6 predicts nothing at 259 Hz, which the window shows
   clearly; H7 reproduces 259 Hz (methyl isotopologue) and the 144-149 Hz
   multiplet with physical couplings (1J 145.4 / 129.8, 2J -4.1 / -4.5,
   3J 6.6 Hz). Inference: prefer H7.
8. User note: natural abundance -> several isotopologues of one molecule.
   Inference: amplitudes are abundances; one CH and one CH3 carbon expected
   1:1 (fitted 0.77); an isopropyl group needs methyl:CH = 2:1.
   Test: isopropyl fit. Outcome: 0.373 and ratio 0.55 -> rejected.
   Absence of other bands -> no other protonated carbon.
9. User note: 15N has low gamma, so its lines are strong.
   Check: simulated 15N-H vs 13C-H amplitude 2.2x; net with abundance 0.75x.
   Test: H7 + alanine 15NH3+ (slow exchange) on 62-320 Hz. Outcome: 0.402 ->
   0.372, 15N amplitude 0.42 (expected 0.35), 1J(N,H) -73.3 Hz.
10. Consistency test: self-consistent alanine (NH protons in every
    isotopologue). Outcome: no better (0.402), abundances wrong (Ca 0.45,
    N 0.91 of Cb). Background check: the 69-75 Hz feature the 15N component
    used is equally present in the no-signal dataset. Outcome: 15N evidence
    withdrawn.
11. Phased cross-check. First with the H7-derived phase (user: phase without
    a hypothesis) -> model-free phase by peak coherence: 0.98 rad, -1.45 ms;
    delay weakly determined (narrow band). The phased spectrum showed a
    ripple; tests (no-signal set, crop 0.10 s, cancelling delay) proved it is
    the record-start baseline rotated by the first-order correction.
    Code fix: background modelled in the record frame (0aeb959).
    Outcome with the fix: H7 0.484 with no bounds; 15N and slow-exchange
    alanine gain only by unphysical parameters -> conclusion unchanged.
12. Conclusion: one CH-CH3 unit (1J about 145 / 129-130 Hz, 3J(HH) 6.6-7.7 Hz),
    no other protonated carbon, no evidence for coupled NH. Lactic acid or
    alanine with fast NH exchange; not separable with these data.
13. Answer revealed: pure L-alanine in water. The conclusion held; the decisive steps were
    6 (window), 8 (abundance test) and 10 (background check). The wrong turns
    were step 3's first methyl reading of the 130 Hz group and step 9's 15N
    improvement; both were corrected by explicit tests, not by residuals.

## Sample 7ad4aafb (strong signal; pyridine-consistent)

1. Observation: one band 148-192 Hz, about ten line groups, no band at 2x.
   Inference: only CH carbons (no methyl, no CH2); spread of 40 Hz ->
   several coupled protons; possibly several 13C isotopologues with different
   1J (aromatic).
2. Hypothesis benzene-13C1 (literature-like J): spans 140-175 Hz, wrong.
   Hypothesis pyridine 13C2/13C3/13C4 (1J about 178/162/162): span and main
   features match. Refinement: 0.144 vs benzene 0.579.
3. Check: symmetry partners (H2-H3 vs H5-H6) diverged in the refinement ->
   symmetry ties are needed; values from memory are not measurements.
4. Model check: transformer skeleton 13C + 5 H, 1J 163.5 Hz (consistent with a
   pyridine isotopologue topology; small couplings unresolved; spin count
   limited to 6 by training).
5. Blind test: global search from the skeleton, 1-3 components, wide bounds.
   Outcome: 1 component 0.338, 3 components 0.212 with unphysical couplings.
   Inference: the search does not find the basin without chemical starts.

## Sample 7dc9a043 (no molecular signal)

1. Observation: nothing above noise in 20-520 Hz for windows 0.3 s to full;
   the isopropylamine reference processed identically shows lines at 100-400x.
2. Tests: exponential baseline instead of SG for the low-frequency region ->
   only a decaying 1.6-2.2 Hz oscillation, also present (weaker) in the
   reference: instrumental or field settling.
3. Conclusion: no J-coupled signal detectable; possible reasons listed
   (no heteronuclear coupling, too weak, different acquisition).
