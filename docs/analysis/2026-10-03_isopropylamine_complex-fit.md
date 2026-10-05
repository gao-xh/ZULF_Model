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
lines reject them. Global refits (all parameters free) from the four best
distinct local candidates:

| global total | 122 Hz window cost | 3J(H1,H2) | 2J(C2,H1) | 3J(C2,H3) | 4J(H2,H3) |
|---|---|---|---|---|---|
| 0.0475 (baseline) | 0.0260 | 5.94 | -1.57 | 5.07 | 0.19 |
| 0.1639 | 0.0228 | 9.30 | -6.13 | 5.24 | 2.24 |
| 0.4198 | 0.0248 | 9.40 | -5.23 | 12.49 | 2.13 |
| 0.4255 | 0.0260 | 9.41 | -5.25 | 12.49 | 2.15 |
| 0.4392 | 0.0420 | 8.41 | -4.97 | 5.88 | 1.93 |

Every refit ends 3.4-9 x higher and gives back most of the window gain; the
window optimum lies in another basin (3J(H1,H2) about 9, 2J(C2,H1) about -6
Hz, near the 9.22 reference values 7.30 / -7.28) that the rest of the
spectrum rejects. Within this model (fast exchange, 13C isotopologues, rate
families) the 121.69 / 122.28 Hz pair cannot be produced by the couplings
without breaking the other lines: a model error. Candidates: NH2 protons
(slow or intermediate exchange) coupled to the methyl-13C cluster, the 15N
isotopologue, a line-shape effect beyond one rate per family, an impurity.

Speed (isopropylamine, one thread): full model 0.08 s, full Jacobian 5.4 s
(one NUFFT pass per rate family: 40 passes with 18 families x 2
isotopologues); window forward 0.006 s, Jacobian 0.86 s with the line band
and 0.21 s with only the free columns.

## NH2 protons and the 15N isotopologue (2026-10-04)

User: add the NH2 protons (slow or intermediate N-H exchange: they couple to
the methyl-13C cluster and can split it) and the 15N isotopologue. Slow model:
NH2 protons kept, 15N@N1 added at natural abundance (gain ratio 0.34 to
13C@C1, fixed by the "ratios" gain model), new couplings J(C2,HN1) (3J
C2-C1-N-H), J(HC1,HN1) (3J H-C-N-H), J(HC2,HN1) (4J), J(C1,HN1) (2J; acts
only on the 13C@C1 lines near 133 Hz, since C1 is 12C in the methyl-13C
isotopologue), J(N1,HN1) start -65, J(N1,HC1) -1, J(N1,HC2) +1 Hz; the rest
from the fast best (ipa_fast_fam2_4j, 0.0475); level-2 rate families.

Line table first (scripts/line_table.py): with the NH couplings at 0 the 122 Hz
lines are unchanged; with J(C2,HN1) 3, J(HC1,HN1) 6, J(C1,HN1) -3,
J(HC2,HN1) 0.5 Hz the methyl lines split into many lines over 120.7-123.2 Hz.
The 15N lines lie at 94-99 Hz (1J(N,H) x 3/2), none near 122 Hz.

Nested baseline (NH couplings 0): total 0.0743 against 0.0475. The 122 Hz
window is unchanged (0.0264 / 0.0260); the whole loss lies at 95-100 Hz
(cost per 5 Hz bin 0.018 against 0.000; scratchpad ipa/slow_where.py), where
the 15N isotopologue puts sharp lines at natural abundance and the data have
none. So the 1J(15N,H) of about 65 Hz is averaged (or broadened away) by N-H
exchange, and couplings of a few Hz to the NH protons would be averaged even
more.

Test 1 (with 15N; scratchpad ipa/local122_slow.py): local fit of 120.9-123.4
Hz with the three NH couplings and the split-type couplings 4J(H2,H3),
3J(H1,H2), 2J(C2,H1), 3J(C2,H3) and the two methyl rate families free, 24
starts (527 s). Window cost 0.0021 against 0.0264 (-92 %) with J(C2,HN1)
-1.33, J(HC1,HN1) 4.41, J(HC2,HN1) -0.56, 3J(H1,H2) 5.07, 3J(C2,H3) 5.43 Hz:
the pair and the dip are reproduced (ipa/local122_slow_best.png). In the full
model (one set of gains for the whole spectrum) the same values give 0.98: the
pair keeps its shape but not its height, the methine lines at 126-140 Hz break
(J(HC1,HN1) and 3J(H1,H2) also act on 13C@C1). Global refits, all free:

| global total | window | J(C2,HN1) | J(HC1,HN1) | J(HC2,HN1) | 3J(H1,H2) | 2J(C2,H1) |
|---|---|---|---|---|---|---|
| 0.0475 (fast best) | 0.0260 | - | - | - | 5.94 | -1.57 |
| 0.0743 (slow, NH = 0) | 0.0264 | 0 | 0 | 0 | 5.94 | -1.57 |
| 0.2857 | 0.0203 | 4.00 | -0.39 | 0.93 | 5.53 | -5.65 |
| 0.3182 | 0.0546 | -3.55 | 4.19 | 0.26 | 6.04 | -4.43 |
| 0.3968 | 0.0319 | -3.36 | 1.84 | -0.29 | 7.56 | -1.24 |
| 0.6645 | 0.0219 | -1.06 | 4.36 | -0.59 | 4.94 | -0.81 |

Test 2 (15N left out with min_ratio 0.3, only the three NH couplings and the
two methyl rate families free in the window, the C-H / H-H couplings held at
the fast best; scratchpad ipa/nh_only.py). The nested start equals the fast
best exactly (0.0475 / 0.0260). Local fit, 24 starts (301 s): best window cost
0.0198 (window forward) with J(C2,HN1) -3.38, J(HC1,HN1) 9.50, J(HC2,HN1)
-1.12 Hz: the NH couplings alone do not open the pair. Global refits:

| start | start total | global total | window | J(C2,HN1) | J(HC1,HN1) | J(HC2,HN1) | J(C1,HN1) |
|---|---|---|---|---|---|---|---|
| NH = 0 | 0.0475 | 0.0475 | 0.0260 | 0.00 | 0.00 | 0.00 | 0.00 |
| local 0.0198 | 0.8500 | 0.1573 | 0.0359 | -3.29 | 9.46 | -1.26 | -0.36 |
| local 0.0274 | 0.9020 | 0.1690 | 0.0194 | -0.82 | 4.62 | -0.95 | -0.68 |
| local 0.0542 | 0.6194 | 0.1994 | 0.0352 | -3.72 | -2.40 | 0.58 | 0.04 |

Every refit with nonzero NH couplings ends 3.3-4.2 x above the fast best.
From NH = 0 the optimizer does not move: a coupling between an otherwise
uncoupled proton and the rest splits every line symmetrically, so the spectrum
is stationary in it at 0 (zero gradient) and a fit started there stays there;
nonzero starts are needed (they were given here and also lose).

Line width (user: narrow the lines of this peak; scratchpad ipa/narrow122.py,
ipa/narrow122.png). Fast best, only the 121-123.5 Hz rate family of 13C@C2
changed, gains / phase / background re-solved:

| R (1/s) | FWHM (Hz) | window | total |
|---|---|---|---|
| 3.16 (fitted) | 1.01 | 0.0260 | 0.0475 |
| 1.58 | 0.50 | 0.0806 | 0.1218 |
| 0.79 | 0.25 | 0.1535 | 0.2912 |
| 0.25 | 0.08 | 0.2017 | 0.5556 |

Narrow lines give one sharp peak at 122.1 Hz (lines 122.127 / 122.153 /
122.364 Hz), never the pair: the model has no line within 121.5-121.9 Hz, and
the fitted 1 Hz width is the compromise that covers both data peaks. The pair
needs a line near 121.7 Hz with the others at 122.2-122.3 Hz, i.e. a position
change of about 0.6 Hz, not a width change.

Two lessons for the local stage. (1) The window forward solves gains, phase
and background again on the window; its cost is not the full-model window cost
(test 2: local 0.0198, the same values in the full model 0.122 in the window),
so local gains overstate what survives globally. Hold the gains and phase at
the global values in the local fit (or fit only one scale). (2) Couplings
starting at 0 to an uncoupled group need nonzero starts (stationary point).

Result: neither the NH2 protons nor the 15N isotopologue explain the 121.69 /
122.28 Hz pair. 15N: lines at 94-99 Hz, absent from the data (evidence for
exchange averaging of the N-H couplings). NH2: they can reproduce the window
only together with changes of the C-H / H-H couplings that the rest of the
spectrum rejects (best global 0.157 against 0.0475), and they are physically
unlikely when the 65 Hz 15N-H coupling is averaged.

## Conclusion

The carbon-skeleton couplings of isopropylamine are the same in both exchange
regimes: 1J(CH) 133.4, 1J(CH3) 124.4, 3J(H,H) 6.1-6.2 Hz. No resolved NH2
structure: no clear 15N line and N-H couplings under 1 Hz, consistent with
fast (or intermediate) N-H exchange in the neat amine. The 13 % gain of the
slow model is not evidence of slow exchange; it can absorb line-shape errors.

The 121.69 / 122.28 Hz pair stays unexplained: not the couplings of the fast
model, not 4J(H2,H3), not the NH2 protons, not the 15N isotopologue (absent
lines at 94-99 Hz), not the line width (narrow lines give one peak at 122.1
Hz). Left: line-shape physics beyond one rate per family, a second species or
an impurity line.

## Open points

- Line widths: done empirically with frequency families (above); next a
  physical relaxation model with few parameters (residual-field g-factor
  broadening, a random-field / dipolar rate, N-H exchange).
- Slow-exchange variant with rate families: done (2026-10-04); with the 15N
  isotopologue it loses (0.0743), without it the NH couplings stay at 0.
- Exchange rate fitted for the NH2 group (Liouville model, D45).
- Uncertainty budget as for the pyridine series.
- Local fit with the gains and phase held at the global values (lesson above).
- 122 Hz pair: impurity / second species check (other amine spectra, a
  spectrum of the same sample at another time).

Commit: see git log (this file).

## Paper-style figure (2026-10-05)

User: show the fit like the methyl-formate panel of a ZULF paper, thinner lines,
wider range, nothing cut, rolling baseline fixed for display (colleague: "the
stuff Blake does to fix the rolling baseline"), detail panels overlaid.
Display processing: crop 0.1 s, the whole 16.2 s record, window 0.1 1/s, zero
fill 4; the fitted model (ipa_fast_fam2_4j) rendered through the same
processing, gains solved on the fit bands and kept for 5-300 Hz. Display
baseline, the same steps for data and model: (1) cubic spline through anchor
points more than 2.5 Hz from every model line and away from narrow data peaks
and power-line harmonics, knots every 6 Hz (zulf_processing.
anchor_spline_baseline, new, with line_mask; test on a 10 Hz ripple under lines
of both signs); (2) only under the line clusters (edges tapered over 1 Hz) an
AsLS 1.5 Hz baseline, which lifts the negative crop-wing valleys between close
lines; elsewhere unchanged. Power-line harmonics (60.06 Hz x n) are kept and
marked. The ripple of period about 1/(crop + delay) = 10 Hz and the broad
negative regions are processing effects of the 0.1 s crop with first-order
phase correction (skills/zulf-phasing), not signal; the fit is unaffected
(complex data, model through the same processing).
Figure: runs/processed/ipa_fast_fam2_4j/paper_style4_w0.1_k6_p2.5_overlay.png
(scratchpad ipa/paper_fig4.py). The 121.69 / 122.28 Hz pair is clearly resolved
in the data and is a single line in the model.

Final display settings (user: keep the flat noise, lift the line clusters):
anchors more than 1 Hz from every model line above 20 Hz, away from the
power-line harmonics and from narrow data peaks within 5 Hz of the lines; knots
every 1.5 Hz; protected runs longer than 2 knot spacings bridged by a straight
line (a cubic piece over a cluster overshot); the noise-free model without
outlier rejection; then AsLS 1.5 Hz within 2.5 Hz of the model lines only.
Figure runs/processed/ipa_fast_fam2_4j/paper_style4_w0.1_k1.5_p1_l2.5_overlay.png.

Script version (repo): scripts/paper_figure.py with the fit's options, --fit
runs/processed/ipa_fast_fam2_4j/fit.json --fid <isopropylamine average FID>
--segments "104,158;222,276" --gains 1,2.5 --colors '{"13C@C1": "#2f8f5b",
"13C@C2 (x2)": "#c0469e"}' --title Isopropylamine --insets with
[{"component": "13C@C2 (x2)", "smiles": "CC(N)C", "atom": 0, "segment": 0,
"rect": [103.6, 0.26, 15.0, 1.0]}, {"component": "13C@C1", "smiles": "CC(N)C",
"atom": 1, "segment": 0, "rect": [142.3, 0.26, 15.5, 1.0]}, {"component":
"13C@C2 (x2)", "smiles": "CC(N)C", "atom": 0, "segment": 1, "rect": [221.8,
0.26, 15.5, 1.0]}]. Structure insets drawn with RDKit (user: the hand-drawn
ones were ugly). Figure runs/processed/ipa_fast_fam2_4j/paper_figure.png.


## Slow exchange with 1J(15N,H) = 83 Hz (2026-10-05, user)

The 119.91 Hz line (17 sigma, full record) is not a power-line line: only
isopropylamine shows it (ethylenediamine, triethylamine, N-ethylmethylamine on
the same instrument do not), the 120.12 Hz harmonic is invisible in all four,
and the model has a methyl-13C line at 119.93 Hz although 119.6-120.4 Hz was
left out of the fit. The missing branch of the pair is the 121.60 Hz line.

Run (scratchpad ipa/run_ipa_slow83.sh, runs/processed/ipa_slow83): --exchange
slow (NH2 kept, 15N isotopologue included), structure one_bond N1 83
(J(N1,HN1) starts at -83 Hz; 15NH2 lines near 124.5 Hz), 4J(H2,H3) and
J(HC2,HN1) free, seeds from ipa_fast_fam2_4j with (J(C1,HN1), J(HC1,HN1),
J(HC2,HN1), J(C2,HN1)) = (-4.5, 5.5, 0.3, 4.5), (0, 4.4, -0.6, -1.3),
(-3, 3, 0.5, 2), (2, 7, -0.5, -3); level-2 families; 8 starts, max_nfev 150.

Coupling diagram (colleague: draw the structure and join groups with lines
whose weight follows J; user: mark 1J and the size of the uncertainty):
scripts/coupling_diagram.py on ipa_fast_fam2_4j, variants ipa_fast_fam1,
ipa_fast_fam2, ipa_fast_fam2_rp, ipa_slow_nested, noise scale sqrt(6.6) = 2.56
(integrated residual correlation length 6.6 points on the 0.042 Hz grid):

| coupling | J (Hz) | sigma (Hz) | noise part | variant part | class |
|---|---|---|---|---|---|
| 1J(C1,H1) | +133.50 | 0.06 | 0.008 | 0.061 | reliable |
| 1J(C2,H2) | +124.43 | 0.07 | 0.009 | 0.068 | reliable |
| 2J(C1,H2) | -4.07 | 0.12 | 0.019 | 0.123 | reliable |
| 2J(C2,H1) | -1.57 | 0.15 | 0.034 | 0.145 | trend |
| 3J(C2,H3) | +5.07 | 0.10 | 0.028 | 0.091 | reliable |
| 3J(H1,H2) | +5.94 | 0.15 | 0.017 | 0.148 | reliable |
| 4J(H2,H3) | +0.19 | 0.03 | 0.033 | 0.004 (2 fits) | trend |

Figure runs/processed/ipa_fast_fam2_4j/coupling_diagram.png. The variant spread
dominates; the 122 Hz misfit is not in these numbers (all variants share it).

1J marked on the structure insets of the paper figure (user): 1J(C2,H2) 124.43 +- 0.07 Hz,
1J(C1,H1) 133.50 +- 0.06 Hz (sigma from coupling_diagram.json).
