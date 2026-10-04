# Isopropylamine: complex known-structure fit with the series fitter (2026-10-03)

- Data: 10000-scan average `/home/user/zulf_data/isopropylamine/average.npy`
  with `0.ini` (NMRduino, 4000 Hz, 65516 points, 16.4 s; sequence
  standard_zf_4000Hz_no_dead.seq; not committed). Isopropylamine is the
  development molecule, never a blind test.
- Question (user): fit isopropylamine with the current fitter (complex data,
  model rendered through the same processing, missing-peak penalty, analytic
  Jacobian, multi-start), and see whether the NH2 group leaves a trace (fast
  against slow N-H exchange).
- Code: scripts/fit_joint_series.py (new: one spectrum with `--shape free`,
  per-entry fit `ranges`, `--exchange slow|fast`; same options in
  reliability_series.py); scratchpad ipa/explore.py, scan.py, apod_view.py,
  make_series.py.

## Processing (checklist of skills/zulf-fid-processing)

1. Sampling rate 4000 Hz (ini); 65516 points.
2. Diagnostics: switching edge 3.497 ms (half height), ADC plateau to
   2.75 ms, ringing to 54.75 ms, first point anomalous.
3. Crop-start scan (record 8 s, no window, zero fill 4, magnitude of the
   strongest line at 133.3 Hz): FWHM 1.09 Hz at 30 ms (ringing, peak/noise 4),
   0.53 Hz at 50 / 75 / 100 ms, 0.56 at 150 ms, 0.62-0.66 at 200-500 ms.
   Chosen: 0.1 s (plateau, 45 ms after the ringing).
4. Window (crop 0.1 s): noise/max 0.0084 / 0.0050 / 0.0043 / 0.0040 at
   0 / 0.3 / 0.6 / 1 1/s; the magnitude FWHM of the 133.3 Hz line grows from
   0.53 to 1.0 Hz at 0.3 1/s because its shoulder merges (apod_view.png).
   Chosen: 0.3 1/s (as for the pyridine series), record 0.1-8.1 s, zero
   fill 3 (0.0417 Hz grid).
5. Phase: instrument calibration (configs/confirmed_samples.json, same
   sequence): delay = edge - 0.033 ms = 3.464 ms, phase0 176.3 deg. The low
   band is absorptive; the 245-257 Hz band alternates in sign (phased.png).
   The complex fit refines a residual zero-order phase in the gain.
6. Bands: 115-160 Hz and 235-262 Hz above noise; mains harmonics at 60.06
   and 240.06 Hz. Fit ranges 85-119.6, 120.4-165, 232-239.6, 240.5-265 Hz
   (120 and 240 Hz left out; 85-115 Hz kept for a possible 15NH2 line near
   1.5 |1J(N,H)|).

## Fits

    python scripts/fit_joint_series.py --series series_ipa.json --real-only false --shape free \
      --exchange fast|slow --range 85,265 \
      --structure '{"motif": "(CH3)2CH-NH2", "one_bond": {"C1": 133, "C2": 125, "N1": -65}}' \
      --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 3 \
      --starts 64 --spread 2.0 --seed 71 --workers 2 --out runs/processed/ipa_fast|ipa_slow

Fast exchange: the NH2 protons dropped (13C at C1 and at the two equivalent
methyl carbons; 6 free couplings). Slow exchange: NH2 kept (also the 15N
isotopologue; 12 free couplings). No priors. Peak tops for the missing-peak
rows: 16. Started 2026-10-03 04:37 UTC.

## Results

| fit | starts | best total | relative residual | next best |
|---|---|---|---|---|
| fast exchange | 64 (618 s) | 0.1456 (2 starts) | 0.376 | 0.3178 |
| slow exchange, generic starts | 64 (1791 s) | 0.1828 | 0.426 | 0.2158 |
| slow exchange, nested start (the fast solution, N-H couplings 0) | 1 (66 s) | **0.1271** | 0.353 | - |

The slow model contains the fast one (N-H couplings 0 decouple the NH2 protons
from the 13C isotopologues), so the generic slow starts missed that basin: their
N-H couplings started at the motif values (3J(H1,N-H) 5.5, 2J/3J(C,N-H)
-4.5 / 4.5 Hz).

Couplings (Hz):

| coupling | fast | slow (nested) |
|---|---|---|
| 1J(C1,H1) (CH) | 133.45 | 133.43 |
| 1J(C2,H2) (CH3) | 124.44 | 124.42 |
| 3J(H1,H2) | 6.11 | 6.18 |
| 2J(C1,H2) | -4.21 | -4.28 |
| 2J(C2,H1) | -1.69 | -1.76 |
| 3J(C2,H3) (to the other methyl) | 5.18 | 5.20 |
| 2J(C1,HN) | - | 0.83 |
| 3J(H1,HN) | - | -0.71 |
| 3J(C2,HN) | - | 0.56 |
| 1J(15N,HN) | - | -64.89 |
| 2J(15N,H1) / 3J(15N,H2) | - | -6.39 / -2.18 |

Decay rates: fast 1.86 (13C at C1) and 3.08 1/s (methyl 13C); slow 1.53, 1.95
and 7.3 1/s (15N). Model delay -3.62 / -3.57 ms (edge -3.46 ms).

All couplings of the carbon skeleton agree within 0.1 Hz between the two
regimes. The slow fit is 13 % lower with 6 more couplings, but uses the freedom
for couplings of the NH2 protons below 1 Hz (slow exchange would give
3J(H1,N-H) of about 5-7 Hz) and for a broad 15N line (7.3 1/s) near 98-100 Hz
that the data do not show clearly (85-112 Hz panel). Remaining misfits are the
same in both: the lines at 121.6 / 122.3 and 133.3 Hz are sharper and taller
in the data, 131.8 Hz is too high in the model, the 251 Hz line is narrower in
the data; one decay rate per isotopologue cannot give lines of different
width.

Figures: runs/processed/ipa_fast/spectrum_fit.png,
runs/processed/ipa_slow_nested/spectrum_fit.png (real part with a display
AsLS baseline subtracted from data and model, imaginary part, 85-112 Hz);
run logs runs/processed/ipa_{fast,slow,slow_nested}/RUN_LOG.md.

## Comparison with the 9.22 subgroup slides (user)

User: the earlier result (9_22_subgroup.pptx, slides 4-20) is clearly better.
That analysis: reference J 1J(CH) 132.79, 1J(CH3) 124.27, 3J(H,H) 7.30,
2J(C1,H2) -7.96, 2J(C2,H1) -7.28, 3J(C2,H3) 2.03 Hz; 13CH + 13CH3 + 15NH2
hypothesis; J fixed, then 114 individual in-band decay times (one per
transition, within-cluster shrinkage beta 0.1), common cluster amplitude and
phase, frozen prediction of validation groups; real part, no AsLS. Its
prediction follows the data closely in both bands; several 15N-hypothesis
transitions sit in 113-129 Hz, and many transitions end at short T2* (about
0.05-0.1 s, some at the bound).

Same objective as here (fast exchange, one decay rate per isotopologue,
scratchpad ipa/fixedJ.py, couplings held, rates and delay refitted):

| couplings | score (data + hard peak rows) | decay rates (1/s) |
|---|---|---|
| this fit | 0.1457 | 1.85, 3.06 |
| 9.22 reference J | 0.7519 | 1.10, 20.0 (methyl at the bound) |
| 9.22 J as start, all free (runs/processed/ipa_fast_from922) | 0.4811 (local minimum) | |

With one width per isotopologue the 9.22 coupling set does not reproduce the
spectrum (the methyl lines are switched off at the rate bound); the closer
match on the slides comes from the per-transition decay times and cluster
amplitudes/phases (and the 15N lines), which reshape the lines at fixed J. A
free width per line can remove lines (as in test A of
2026-10-01_blake-pyridine_peak-penalty.md). The two results are therefore not
comparable by their figures; a fair comparison needs the same line-shape
freedom on both coupling sets.

## Several decay rates per isotopologue (2026-10-04)

User: release the one-rate-per-isotopologue limit. The solver already splits
transitions into rate families by frequency (`ParameterPolicy.family_edges_hz`,
analytic derivatives, D31); fit_joint_series / reliability_series now take
`--family-edges` and `--rate-bounds` (test: three families per isotopologue,
Jacobian against finite differences). Fitted transitions of the fast model
(scratchpad ipa/transitions.py): 13C@C1 at 131.8-138.7 Hz (plus low
frequencies < 6 Hz), methyl 13C at 117-130 Hz and 235-256 Hz (plus < 17 Hz).

- Level 1, edges 50, 190 Hz: per isotopologue one rate below 50 Hz, one for the
  J band, one for the 2J band (3 x 2 rates).
- Level 2, edges 50, 118.5, 121, 123.5, 126, 128.8, 131, 132.7, 134.2, 136,
  137, 138, 190, 241, 245.5, 249.8, 253.5 Hz: one rate per line cluster (18
  families per isotopologue; empty families have no effect).

Rate bounds 0.2-15 1/s (T2* >= 67 ms, so a line cannot be switched off by an
extreme width). Seeds: this fit's J and the 9.22 reference J (same line-shape
freedom for both), 14 perturbed starts (spread 1 Hz), seed 73, outputs
runs/processed/ipa_fast_fam1 and ipa_fast_fam2. Started 2026-10-04.

| rates | best total | relative residual | from the 9.22 J | next best |
|---|---|---|---|---|
| one per isotopologue (above) | 0.1456 | 0.376 | 0.48 (free) | 0.3178 |
| level 1 (3 per isotopologue), 16 starts, 680 s | 0.1069 | 0.323 | 0.4648 | 0.2908 |
| level 2 (one per line cluster), 16 starts, 2731 s | **0.0500** | **0.223** | 0.5413 | 0.1402 |

Couplings move by at most 0.16 Hz from the one-rate fit (level 2: 1J(CH)
133.50, 1J(CH3) 124.41, 3J(H,H) 5.95, 2J(C1,H2) -4.08, 2J(C2,H1) -1.47,
3J(C2,H3) 5.20 Hz). Started from the 9.22 reference J with the same freedom,
the fit ends 9-11 x higher, so the better figure on the slides was the line
shapes, not the couplings.

Level-2 rates of the families that hold lines (1/s): 13C at C1 131.8-132.0 Hz
3.12, 133.3-133.45 Hz 1.33, 134.8-134.9 Hz 2.05, 136.6 Hz 2.35, 137.2 Hz
1.67, 138.7 Hz 2.77; methyl 13C J band 117.1 Hz 2.82, 119.6-119.8 Hz 1.78,
122.0-122.7 Hz 3.10, 124.8 Hz 1.82, 127.2-127.7 Hz 2.78, 129.9 Hz 3.08; 2J
band 235-240 Hz 7.11, 242-245 Hz 4.77, 246-248.5 Hz 4.81, 251-253 Hz 5.45,
254-256 Hz 5.28. No rate of a family with lines in the fit ranges is at a
bound (the methyl family below 50 Hz is at 14.7, outside the fit ranges and
undetermined; empty families do nothing). The 2J-band lines decay about twice
as fast as the J-band lines of the same isotopologue (1.8-3.1 against
4.8-7.1 1/s); within the methine 1.3-3.1 1/s. A transition-dependent
relaxation is expected physically (rates depend on the eigenstates of each
coherence; residual-field broadening scales with each line's g factor), so
the frequency families are an empirical stand-in for a relaxation model.
Remaining misfits: the 121.6 / 122.3 Hz pair (sharp dip between them in the
data), a feature at 239.4 Hz and the 117 Hz shoulder. Figures
runs/processed/ipa_fast_fam{1,2}/spectrum_fit.png.

## The 121.69 / 122.28 Hz pair (2026-10-04)

User: check from the simulation whether the pair exists. Data (crop 0.1 s, real
part): maxima 121.69 (0.56) and 122.28 Hz (0.99 of the band maximum), dip
121.97 Hz (0.20 at 0.3 1/s; 0.07 with the full record and no window: two
separate narrow lines). Fitted model (level 2): methyl-13C transitions at
121.85 (0.58), 122.00 (0.47) and 122.68 Hz (0.26); with narrow lines (all
rates 0.3 1/s) they form a strong peak near 122.1 and a weak one near 122.45 Hz
with a dip between (runs/processed/ipa_fast_fam2/doublet_122.png; scratchpad
ipa/doublet.py). The model has a doublet there, 0.3 Hz too high and with the
intensities reversed; line widths cannot fix positions.

df/dJ of the three lines: 2J(C2,H1) -0.26 / -0.11 / +0.07, 3J(C2,H3)
-0.26 / -0.30 / -0.43, 3J(H1,H2) -0.21 / -0.11 / +0.20, and the methyl-methyl
4J(H,H) -0.27 / -0.48 / -0.84 Hz per Hz. That coupling is not in the motif and
was fixed at 0. Running: level-2 fit with 4J(H2,H3) free (starts -0.8, 0,
+0.8 Hz from the level-2 best), runs/processed/ipa_fast_fam2_4j.

Result (1125 s): all three starts end at 4J(H2,H3) = +0.19 Hz, total 0.0475
(-5 % against 0.0500), residual 0.217; the other couplings move by at most
0.13 Hz (3J(C2,H3) 5.07). The methyl-methyl coupling is small and does not
make the pair.

Residual-peak rows (user: penalise peaks in the residual, both signs; a
spread-out residual hardly matters, a localized one marks a structural
mismatch). New in fit_joint_series: `--residual-peaks STRENGTH` with
`--residual-peak-sigma/-halfwidth/-assign/-rounds/-candidates`
(JointSeries.find_residual_peaks: prominence of |model - data| over a robust
local noise level, narrower than 1.5 Hz, assignable when a model transition
lies within 0.6 Hz; outer refit loop with the windows weighted; candidates
ranked on one common window set; tests for a missing line, a negative feature,
a far line reported only, a broad misfit ignored, Jacobian). Running on the
4J-free level-2 solution: strength 3, 1 candidate, 3 rounds,
runs/processed/ipa_fast_fam2_rp.

Result of the residual-peak run (2684 s): peaks per round 0: 122.13 (9.7 sigma,
model above data), 133.25, 116.88, 242.58 Hz; 1: 121.67 (model below, not
assignable: the lines had moved away), 133.67, 116.83; 2: 122.13, 133.25,
116.88; final: 121.67 (not assignable), 116.83, 124.79. Chosen: the refitted
solution, plain objective 0.0689 (+45 % against 0.0475), weighted 0.2425
(-9 %); couplings within 0.15 Hz (1J(CH3) 124.55). The pair is not resolved:
emptying the dip at 122.1 Hz leaves the data line at 121.69 Hz unmatched. Two
weaknesses of the stage: assignability flips when lines move (keep windows,
judged on the starting transitions) and the ranking accepted a large plain
loss (cap it).

Source attribution (JointSeries.peak_sources; scratchpad ipa/sources.py) of the
four peaks of the 4J-free solution: 122.12 Hz: methyl-13C lines 122.13 /
122.15 / 122.36 Hz, split-type couplings 4J(H2,H3) 0.25, 3J(H1,H2) 0.20,
2J(C2,H1) 0.17, 3J(C2,H3) 0.11 Hz/Hz (1J(CH3) only shifts); 133.25 Hz:
methine lines 133.32 / 133.50 Hz, 2J(C1,H2) 0.31, 3J(H1,H2) 0.16; 116.88 and
242.58 Hz: single methyl lines (shift only). Line table of the three 122 Hz
lines (scripts/line_table.py --second-order): contributions e.g. 122.153 Hz =
124.33 (1J(CH3)) - 1.77 (3J(H1,H2)) - 0.91 (3J(C2,H3)) + 0.54 (2J(C2,H1)) Hz;
the two lower lines are 0.026 Hz apart (near-degenerate), second-order and
cross terms up to 0.1 1/Hz, and a 1 Hz Taylor step is off by 0.15 Hz even at
second order.

Targeted local fit of the 122 Hz window (120.9-123.4 Hz; free 4J(H2,H3),
3J(H1,H2), 2J(C2,H1), 3J(C2,H3) and the two methyl rate families; 27 starts;
window-only forward with line band +-2 Hz and Jacobian of the free columns:
257 s, where the full-spectrum version had not finished after 50 min;
scratchpad ipa/local122.py): window cost 0.0054 against 0.0260 (-79 %) with
4J(H2,H3) 1.43, 3J(H1,H2) 8.18, 2J(C2,H1) -4.41, 3J(C2,H3) 6.18 Hz; held at
these values the whole spectrum scores 0.92-1.09 (baseline 0.0475): the other
lines reject them. Global refits from the local candidates: running.

Speed (isopropylamine, one thread): full model 0.08 s, full Jacobian 5.4 s
(one NUFFT pass per rate family: 40 passes with 18 families x 2
isotopologues); window forward 0.006 s, Jacobian 0.86 s with the line band
and 0.21 s with only the free columns.

## Conclusion

The carbon-skeleton couplings of isopropylamine are the same in both exchange
regimes: 1J(CH) 133.4, 1J(CH3) 124.4, 3J(H,H) 6.1-6.2 Hz. No resolved NH2
structure: no clear 15N line and N-H couplings under 1 Hz, consistent with
fast (or intermediate) N-H exchange in the neat amine. The 13 % gain of the
slow model is not evidence of slow exchange; it can absorb line-shape errors.

## Open points

- Line widths: done empirically with frequency families (above); next a
  physical relaxation model with few parameters (residual-field g-factor
  broadening, a random-field / dipolar rate, N-H exchange).
- Slow-exchange variant with rate families (does the 13 % gain survive?).
- Exchange rate fitted for the NH2 group (Liouville model, D45).
- Uncertainty budget as for the pyridine series.

Commit: see git log (this file).
