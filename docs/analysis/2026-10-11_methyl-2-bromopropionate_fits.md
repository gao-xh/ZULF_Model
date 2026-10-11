# Methyl 2-bromopropionate: known-structure fits, passes 1-7, couplings scan and a mixture test (2026-10-10 / 11)

- Data (drive only): `original/2026-04-12-methyl-2-bromopropionate/Methyl 2-bromopropionate/`, 8296 scans, 4000 Hz,
  zero field (14N_amines set of the collaborator). z5 average:
  ~/research/zulf/data/processed/2026-04-12-methyl-2-bromopropionate/z5/; whole-grid series runs/series/batch_methyl-2-bromopropionate
  (crop 0.1 s, record 8 s, window 0.3 1/s, grid 20-380 Hz).
- Question (Xuehan): fit the known structure CH3-CHBr-COO-CH3; then improve the fit (triplet at 155-158 Hz, fine
  structure of the methyl bands).
- Structure: chain from SMILES CC(Br)C(=O)OC (C1 methyl, C2 CHBr, C3 carbonyl, C4 OCH3); natural-abundance 13C at C1,
  C2, C4; bromine as a site without nuclear spin.
- Code: batch passes with scripts/run_batch.py (configs/batch_2026-10-10*.json); single runs listed below.

## Passes

| pass | change | stage-2 objective | 1J C1 / C2 / C4 (Hz) | result |
|---|---|---|---|---|
| 1 (batch) | 1J starts from the sp3 rule (125 Hz everywhere) | - (not run, superseded) | - | - |
| 2 (`_v2`) | neighbour rule (D61 amendment): 125 / 156 / 141 | 0.194 | 123.46 / 153.33 / 140.09 | C1 and C4 bands missed; C2 at the 155-158 Hz lines |
| 4 (`_v4`) | 1J starts from the observed bands (D64): 130.9 / 147.0 / 147.65 | 0.209 | 130.54 / 147.07 / 147.60 | CH3 J and 2J bands and OCH3 2J band fitted; the 155-158 Hz triplet missed (C2 under the OCH3 J line at 147.6 Hz) |
| 5 (`_v5`, stage 1 only) | C2 started at 153.5 Hz, ratios fixed | stage 1 0.411 (pass 4: 0.334) | 130.50 / 155.39 / 147.65 | worse with natural-abundance amplitude ratios; stopped |
| multi-start probe | candidates per site, joint scoring, short refits | - | prefers C2 147.0 (0.354) over 155.85 (0.396) | with fixed ratios the data prefer C2 near 147 Hz |
| 6 (`_v6`) | C2 from the triplet (155.4 Hz), free isotopologue amplitudes (`--variant free`) | **0.0387** | 130.67 / 155.78 / 147.66 | triplet fitted; amplitude of 13C@C2 0.025 vs 0.036 / 0.034 (about 72 %): the carbon bound to Br loses signal (scalar relaxation by the fast-relaxing quadrupolar 79Br / 81Br) |
| 6, stage 3 (`_v6_s3`) | rates <= 10 1/s, 8 starts, component search start and end | 0.0384 | 130.68 / 155.81 / 147.66 | no gain: not an optimisation problem |
| 7 (`_v7_longrange`) | long-range couplings across the ester added (free, start 0.5 Hz, no prior) | **0.0282** | 130.68 / 155.86 / 147.67 | third triplet line (158.2 Hz) fitted; 4J(C4,H2) -0.25, 4J(C2,H4) +0.45, 5J(C1,H4) -0.27, 5J(H2,H4) +0.26 Hz |

Pass 7 is the best fit (overlay figure runs/figures/overlay_metabolites.png; scripts/overlay_panels.py RUNS).
Other couplings of pass 7: 2J(C1,H2) -2.56, 3J(H1,H2) 6.48, 2J(C2,H1) -4.46 Hz. Conditional numerical results of
the given structure, not a determined assignment.

## Checks of what remains

Residual of pass 7 (rms / noise per window: 120-135 Hz 5.4, 135-152 3.9, 152-165 9.9, 240-256 1.4, 256-270 3.9,
285-305 5.9). With the model lines drawn (runs/figures/batch/methyl-2-bromopropionate_lines.png) two experimental peaks
have no model line: 128.0 Hz (C1 methyl J band; nearest model line 128.5 Hz) and 263.0 Hz (C1 methyl 2J band;
nearest 263.5 Hz); a shoulder at 294.5 Hz (OCH3 2J band) is also missing.

| hypothesis | test | outcome |
|---|---|---|
| wrong combination of the C1 methyl couplings | grid scan 2J(C1,H2) -8..8 Hz x 3J(H1,H2) -2..12 Hz, step 0.25 Hz, all else held (runs/figures/batch/mbp_c1_scan.png; runs/figures/batch/coupling_scan.py) | current values the only deep minimum (0.029; next 0.110), also in the C1 bands alone: excluded |
| missing long-range couplings | pass 7 | helps the C2 triplet, not 128.0 / 263.0 Hz; the couplings that reach C1 (5J, 6J) are far too small for a 0.5-1 Hz shift |
| bromine as a coherent spin | not done: 79Br / 81Br relax in microseconds (quadrupolar) and are self-decoupled; no Br splittings are expected, and the solver has no quadrupolar relaxation | its effect is the weaker C2 isotopologue (pass 6) |
| slow O-H / N-H exchange | not applicable: no exchangeable protons | - |
| second species: 2-bromopropionic acid (hydrolysis) | ester + acid as one model (`--structure` list; runs/processed/mbp_mix_acid_s1, _s2; runs/processed/mbp_mix_acid.sh), free amplitudes | stage 1 0.102; the acid's part is broad and does not reach 128.0 / 263.0 Hz; acid amplitudes as large as the ester's and 2J(C1,H2) -12 Hz (implausible); no methanol band near 282 Hz as hydrolysis would give: very unlikely (stage 2 pending) |

## Conclusion

Conditional numerical result for CH3-CHBr-COO-CH3 at zero field: 1J(C1,H) 130.7, 1J(C2,H) 155.9, 1J(C4,H) 147.7 Hz,
2J(C1,H2) -2.6, 2J(C2,H1) -4.5, 3J(H1,H2) 6.5 Hz, long-range couplings across the ester below 0.5 Hz, the
Br-bound carbon's isotopologue at about 72 % of the natural-abundance amplitude. Two peaks (128.0, 263.0 Hz) and a
shoulder (294.5 Hz) remain unexplained by this structure.

## Lessons (recorded in DECISIONS D64 and skills)

- 1J starts off by 5-8 Hz trap the fit; the band scan (--one-bond-from-data) fixes it, but overlapping bands of
  different sites can mislead it (C2 locked on the OCH3 line); the multi-start probe scores combinations jointly.
- A carbon bound to a quadrupolar nucleus (Br, Cl, 14N) can be weaker than natural abundance predicts: free the
  amplitudes when its band is clearly visible but not fitted with fixed ratios.
- The chain rule stops at three bonds; across esters / ethers long-range couplings of a few tenths of a hertz matter.

## Open points

- Origin of 128.0 / 263.0 / 294.5 Hz: sample purity and preparation (collaborator); a high-field 1H spectrum would show
  impurities directly.
- Stage 2 of the mixture test (running when this was written).

Commit: (this commit)
