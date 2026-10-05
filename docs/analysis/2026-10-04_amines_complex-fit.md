# Ethylenediamine, triethylamine, N-ethylmethylamine: complex known-structure fits (2026-10-04)

- Data: confirmed samples fde3fbb2 (ethylenediamine, H2N-CH2-CH2-NH2), 4322bdfc
  (triethylamine, N(CH2CH3)3) and e66a4b08 (N-ethylmethylamine,
  CH3-CH2-NH-CH3); NMRduino average FIDs, 4000 Hz, 65516 points (uploads, read
  from ZULF_DATA_DIR, not committed). Identities given by the user after the
  blind analyses (docs/ANALYSIS_LOG.md); these are known-structure fits.
- Question (user): fit the other amines with the current fitter, as for
  isopropylamine (docs/analysis/2026-10-03_isopropylamine_complex-fit.md):
  fit quality, couplings, NH / 15N traces, and whether localized misfits like
  the isopropylamine 121.69 / 122.28 Hz pair also appear.
- Code: aa07a76 at the start; run logs runs/processed/{eda,tea,nema}_fast*/RUN_LOG.md;
  scratchpad amines/make_series.py, auto_edges.py, run_stage1.sh, run_stage2.sh.

## Processing (same as isopropylamine)

Crop 0.1-8.1 s (start sample 400), exponential 0.3 1/s, zero fill 3 (0.0417 Hz
grid), phase from the instrument calibration (configs/confirmed_samples.json:
phase0 176.3 deg at the switching edge - 0.033 ms); edges 3.470 / 3.439 / 3.477
ms. The complex fit refines a residual zero-order phase in the gain. Fit ranges
(mains harmonics n x 60.06 Hz +-0.4 Hz out; the regions where 15N lines would
lie kept even without signal; the broad low-frequency features below 55 Hz,
which no model line explains, left out):

| sample | ranges (Hz) |
|---|---|
| ethylenediamine | 85-119.6, 120.4-179.6, 180.5-237.5 |
| triethylamine | 85-119.6, 120.4-142.1, 180.5-215.1, 231.6-239.6, 240.5-263.2 |
| N-ethylmethylamine | 55-59.6, 60.5-77.3, 102.2-119.6, 120.4-153.5, 177.9-179.6, 180.5-215.9, 231.6-239.6, 240.5-276 |

Model lines (line table, motif defaults): ethylenediamine 13C@C1 (x2) at
190-206 Hz, slow exchange 15N lines at 98-110 Hz; triethylamine 13C@C1 (CH2)
185-196 Hz, 13C@C2 (CH3) 123-127 and 236-258 Hz; N-ethylmethylamine 13C@C1
(CH3 of ethyl) 123-127 / 237-258, 13C@C2 (CH2) 192-206, 13C@C3 (N-CH3)
128-134 / 259-266 Hz, slow exchange 15N lines at 62-67 Hz.

## Stage 1: fast exchange, one rate per isotopologue

fit_joint_series --real-only false --shape free --exchange fast, signal
threshold 2.5 / taper 4, peak penalty 5 (smooth 0.03, min sigma 3), 32 starts
(spread 2 Hz) from the earlier regression couplings (runs/regression/final;
2J(C2,H1) of N-ethylmethylamine started at -4.5 instead of -11.9 Hz):

| sample | objective | relative residual | time |
|---|---|---|---|
| ethylenediamine (eda_fast) | 0.1015 | 0.302 | 266 s |
| triethylamine (tea_fast) | 0.3115 | 0.544 | 710 s |
| N-ethylmethylamine (nema_fast) | 0.3129 | 0.548 | 700 s |

(isopropylamine at this stage: 0.146.) Couplings (Hz):

- ethylenediamine: 1J(C,H) 132.40, 2J(C,H) -2.93, 3J(H,H) 5.63.
- triethylamine: 1J(CH2) 130.71, 1J(CH3) 124.97, 2J(C1,H2) -4.61, 2J(C2,H1)
  -3.18, 3J(H1,H2) 7.22, 3J(C1,HX) (through N) -2.21.
- N-ethylmethylamine: 1J(CH3) 124.91, 1J(CH2) 138.33, 1J(N-CH3) 131.76,
  2J(C1,H2) -3.49, 2J(C2,H1) -5.91, 3J(H1,H2) 7.25, 3J(C2,H3) 11.45 (through
  N, implausibly large), 3J(C3,H2) 3.73.

Misfits (runs/processed/*_fast/spectra.png): one rate per isotopologue cannot
carry both the narrow and the broad lines (triethylamine: 125 Hz lines too
weak, 250 Hz too strong; N-ethylmethylamine: 195-205 Hz CH2 lines too weak).

## Stage 2: rate families

Running (run_stage2.sh): edges from the stage-1 lines (auto_edges.py: lines of
relative amplitude >= 0.05 clustered at gaps > 0.5 Hz, at most 18 families),
remote couplings at 0 in the motif freed (triethylamine J(H1,HX), J(H2,HX),
J(C2,HX); N-ethylmethylamine J(C1,HC3), J(HC1,HC3), J(HC2,HC3), J(C3,HC1);
seeds 0 and +-0.8 Hz), rate bounds 0.2-15 1/s, 4 starts.

Ethylenediamine stage 2 (eda_fast_fam, 4 families 195.3 / 199.5 / 202.4 Hz,
123 s): objective 0.0719, relative residual 0.253; 1J 132.06, 2J -2.86, 3J(H,H)
6.35 Hz. The data show sharp pairs (190.6 / 192.1, 197.8 / 198.8, 200.2 Hz,
widths 0.3-0.5 Hz); the model has five lines (190.6, 198.09, 198.55, 199.5,
204.6 Hz) and covers the pairs with broad lines (rates 7.9 / 5.2 / 11.0 / 2.5
1/s, FWHM 2.5 / 1.6 / 3.5 / 0.8 Hz): broadening is too cheap where positions
are missing. The fitted delay is off: +2.3 ms (stage 1) and the -10 ms bound
(stage 2), against -3.6 ms for isopropylamine, triethylamine and
N-ethylmethylamine (near 200 Hz a delay is ambiguous by about 1 / f inside one
band). New option `--phase-delay-bounds lo,hi` (ms) in fit_joint_series
(instrument prior -4.6,-2.6); rerun eda_fast_d / eda_fast_fam_d running.
Structural candidates next: two vicinal couplings J, J' of the AA'BB' CH2-CH2
unit (line table: J 3 / J' 9.7 Hz moves the 199.5 Hz line to 200.14 Hz, the
data line is at 200.2 Hz; no line at 192.1 Hz yet) and slow NH2 exchange.

Ethylenediamine with the delay prior (`--phase-delay-bounds=-4.6,-2.6`; note
the `=`: argparse reads "-4.6,..." as an option otherwise): stage 1
(eda_fast_d) 0.1020, delay at the -2.6 ms bound; stage 2 (eda_fast_fam_d)
0.0736 (against 0.0719 free), delay at the -4.6 ms bound, rates 7.2 / 5.6 /
9.7 / 2.7 1/s, 1J 132.24, 2J -2.97, 3J(H,H) 5.86 Hz. The delay is not the cause
of the broad lines; it drifts because the model cannot place the pairs.

## Stage 2 results (rate families)

| sample | families | objective | relative residual | delay (ms) | time |
|---|---|---|---|---|---|
| ethylenediamine (eda_fast_fam) | 4 | 0.0719 | 0.253 | -10 (bound) | 123 s |
| triethylamine (tea_fast_fam) | 18 | 0.0490 | 0.215 | -3.54 | 1003 s |
| N-ethylmethylamine (nema_fast_fam) | 19 | 0.0752 | 0.276 | -3.55 | 2099 s |

(isopropylamine best 0.0475 / 0.217.) Couplings (Hz):

- triethylamine: 1J(CH2) 130.75, 1J(CH3) 125.00, 2J(C1,H2) -4.66, 2J(C2,H1)
  -3.06, 3J(H1,H2) 7.20, 3J(C1,HX) -2.15, J(H1,HX) -0.20, J(H2,HX) 0.43,
  J(C2,HX) -0.50.
- N-ethylmethylamine: 1J(CH3) 124.93, 1J(CH2) 138.37, 1J(N-CH3) 131.74,
  2J(C1,H2) -3.37, 2J(C2,H1) -6.20, 3J(H1,H2) 7.20, 3J(C2,H3) 12.08 (through
  N; implausible, the 190-200 Hz CH2 lines are poorly fitted), 3J(C3,H2) 4.24,
  remote couplings within 0.3 Hz of 0.

Figures: runs/processed/{eda_fast_fam,tea_fast_fam,nema_fast_fam}/spectrum_fit.png
(scratchpad amines/plot_fit.py).

A common misfit near 121 Hz. Real-part peaks in 120.4-123 Hz: isopropylamine
121.67 / 122.29 Hz, triethylamine 121.08 Hz (31 sigma), N-ethylmethylamine
121.33 Hz (55 sigma), ethylenediamine (no CH3) nothing above 3 sigma: not an
instrument line, and present in triethylamine, which has no N-H. It belongs to
the methyl-13C isotopologue of a CH3-CHn unit: the models put lines there
(triethylamine 13C@C2 121.01-121.73 Hz, amplitudes 0.08-0.34;
N-ethylmethylamine 13C@C1 121.22-121.87 Hz) but spread over 0.5-0.7 Hz, and
the fit covers the sharp data line with a broad family. So the isopropylamine
pair is one case of a systematic misfit of this cluster; NH2 protons are not
needed for it (triethylamine has none).

## Tools added during this analysis (2026-10-05)

User: an interface to adjust J by hand (coarse and fine), start an automatic
refinement, and see the residual and the score live while dragging. Added
scripts/j_tuner.py + j_tuner_ui.html (docs/WORKFLOW.md section 8, item 7),
fit_joint_series.make_parser / build_problem (the tuner builds the fitter's own
problem), scripts/make_series_entry.py (FID -> series entry with this recipe;
reproduces the triethylamine series exactly). Checked on isopropylamine
(ipa_fast_fam2_4j): the page shows objective 0.047526 as the fit; 3J(H1,H2)
+0.25 Hz by the fine slider gives 0.1432 (+201 %); a refinement of the small
couplings and the rate families in view returns to 0.047526 in 11 evaluations
(20 s). Update time about 0.3 s per change (residual peaks included).

Fit monitor (user: show the output of a running fit without changing the
algorithm): fit_joint_series now records every run in OUT/monitor;
scripts/fit_monitor.py shows it live. Checked: the ethylenediamine stage-2 run
(4 starts, 2 workers, max_nfev 15) gives identical scores, couplings and
spectrum parameters with --monitor on and off; demo on triethylamine
(runs/processed/tea_monitor_demo, 6 starts, 3 workers, max_nfev 60).

## Slow N-H exchange with 1J(15N,H) = 83 Hz (2026-10-05, user)

User: ethylenediamine and triethylamine are not good; fit the NH2 with
1J = 83 Hz; N-ethylmethylamine again as well. Triethylamine has no N-H (no
NH test there). Line table (ethylenediamine, slow, 1J(N,H) -83): with the NH
couplings at 0 the 13C@C1 lines are those of the fast model and the 15N@N1
lines sit at 123.3-125.7 Hz (1.5 x 83); with 2J(C,HN) -3, 3J(C,HN) 1.5,
3J(H,HN) 3 Hz the methylene-13C lines split into about 90 lines over
187.7-207.7 Hz, the 15N lines spread over 120.8-127.3 Hz.

Runs (run_slow83.sh): --exchange slow, 15N isotopologue included, NH coupling
seeds nonzero (0 is a stationary point): ethylenediamine (eda_slow83; seeds
(2J(C1,HN1), 3J(C1,HN2), 3J(HC1,HN1)) = (-4.5, 4.5, 5.5), (-3, 1.5, 3),
(3, -1.5, -3), (-1.5, 1, 7); 13 family edges incl. 110 / 130 / 150 Hz for
the 15N lines; 8 starts) and N-ethylmethylamine (nema_slow83; series with the
77-102 Hz range added for the 15NH line near 83 Hz; 4 seeds, 3J(C2,H3)
reset to 4.5 in three; 20 families; 8 starts, max_nfev 120). Timing: residual
0.05 / 0.08 s, Jacobian 1.7 / 7.7 s.

Slow-exchange results (1J(15N,H) -83 Hz): ethylenediamine (eda_slow83, 1043 s)
best 0.0435 (fast 0.0719), relative residual 0.195; 1J(C,H) 131.39, 2J(C,H)
-3.71, 3J(H,H) 4.53, 3J(H,HN) 6.33, 2J(C,HN) -0.09, 3J(C,HN) -0.97,
1J(N,H) -83.1 Hz (next starts 0.0515, 0.0725). N-ethylmethylamine (nema_slow83,
4245 s) best 0.1653 against 0.0752 fast: rejected. 1J(N,H) of aliphatic amines
is about 64-67 Hz in the literature (85 Hz is typical of sp2 N), so a gain at
83 Hz may come from 15N lines hidden in a line cluster; to check against a run
with 65 Hz.

Sharp structure rows (user: stop fitting several splittings with one broad
line, a broad line refined is fine; penalise narrow negative residual peaks):
fit_joint_series --dip-penalty (rows at sharp data valleys) and
--peak-max-width (rows only at sharp extrema); tests in test_joint_series.
Test run on ethylenediamine (fast, families): eda_fast_fam_sharp with
--peak-penalty 20 --dip-penalty 10 --peak-prominence 0.02 --peak-min-sigma 5
--peak-max-width 1.0.


## The 15N-H line test and a weighting gap (2026-10-05)

**The 15N-H line test.** In slow N-H exchange, the 15N isotopologue must show the 15N-H line, with its intensity
fixed by natural abundance (ratio 0.338 to a 13C component). The fast-exchange model has no such line, so the slow
model does not nest the fast one, even with every NH coupling at 0.

- **N-ethylmethylamine (one NH).** The line falls at |1J(N,H)|.
  - Slow model at the nested start (fast-family couplings and rates, NH couplings 0): the line is 0.112 at 83 Hz
    and 0.087 at 65 Hz. The strongest 13C lines are 0.042.
  - Data: 0.0012 / 0.0020. The largest peak between 62 and 102 Hz is 6-9 sigma, against 224 sigma for the main
    lines.
  - Slow exchange is therefore excluded at either value, which explains nema_slow83 (0.165). The planned nested
    slow runs (nema_slow83n / nema_slow65n) were cancelled.
- **Ethylenediamine (NH2).** The 15N lines fall near 1.5 x J: 124.5 Hz for 83, 97.5 Hz for 65.
  - Data: at most 3 sigma at 97.5, 124.5, 83 and 65 Hz.
  - eda_slow83 nevertheless puts its 15N lines at 118-128 Hz, up to 0.0085, against a data maximum of 0.0019 there
    (residual 3.8 x the data norm in that window). It also puts ripples at 210-238 Hz that the data do not show.
    Figure: runs/processed/amine_overview/ethylenediamine.png.
  - Its objective (0.0435 against 0.0719) improves anyway. The gain comes from the 13C lines at 185-210 Hz, where
    the free NH couplings absorb misfit; the residual there drops from 0.25 to 0.19.

**The weighting gap.** fit_joint_series weights by data signal only. Far from data peaks, the relative weight is
signal_outside_weight 0.2, so a squared residual there counts 0.04. A model line where the data show none is
therefore cheap. The single-spectrum solver closes this gap with signal_model_passes: the model's own lines join
the cores, so a line placed where the data show none is fully penalised. The joint fitter never had that pass.

**Consequence.** The slow-exchange gains of eda_slow83 and ipa_slow83 are not evidence for slow exchange or for
1J(N,H) = 83 Hz until they are rescored with model-line weighting. The eda_slow65 control (same recipe, 65 Hz) is
running for comparison.

Overview figures of the current fits: runs/processed/amine_overview/{ethylenediamine,triethylamine,
N-ethylmethylamine}.png (scratchpad amines/plot_runs.py).
