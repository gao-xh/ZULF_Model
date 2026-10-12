# Amines: refits with background off, 1J from bands, noise filter (2026-10-11)

- Question (Xuehan): improve the amines overlay with the improvements found on the metabolites.
- Command: `python scripts/run_batch.py configs/batch_2026-10-10_v13.json --parallel 3 --workers 2 --tag _v13
  --table docs/analysis/2026-10-10_batch_results_v13.md` (all entries: `background_order -1`,
  `one_bond_from_data`, 1J starts from the previous best pass; window 1.5 1/s and record 4 s for diethylamine,
  n,n-dimethylethylenediamine, sec-butylamine). Applied to every entry as a probe; the overlay switched only
  where better.
- sec-butylamine stage 2 (6.4 h, 24 evaluations) and dipropylamine stage 1 (7 h) were stopped by PID;
  sec-butylamine fit.json is the best monitor point (objective 0.0781, `stopped_early`); dipropylamine has no v13.
- Figures: runs/figures/batch/amines_v13_compare.png (previous vs v13; runs/figures/batch/compare_versions.py),
  runs/figures/overlay_amines.png (version printed in every panel).

## Results

| entry | previous | v13 objective | verdict |
|---|---|---|---|
| tert-butylamine | v1 s2 (0.121) | 0.169 | no gain, keep v1 |
| isopropylamine | v2 s2 (0.063) | 0.133 | no gain, keep v2 |
| triethylamine | v1 s2 (0.090) | 0.090 | same, keep v1 |
| ethylenediamine | v1 s2 (0.062) | 0.078 | small lines 190-193 Hz no better, keep v1 |
| diethylamine | v2 s2 (0.559, almost no model signal) | 0.855 (other data vector) | **v13**: 125, 251, 255 Hz lines fitted; model lines near 200 Hz not seen in the data |
| n,n-dimethylethylenediamine | v4 s2 (0.124) | 0.074 | **v13**: all bands follow the data |
| sec-butylamine | v2 s2 (0.070) | 0.078 (stopped early) | **v13**: bands follow the data with less noise |
| dipropylamine | v2 s2 (0.066) | - | keep v2 |

Objectives of passes with the noise filter are not comparable with the 0.3 1/s passes.

## Open points

- The 2J (XA3) band near 245-255 Hz is too weak in the model for isopropylamine, triethylamine and (less)
  tert-butylamine, in every pass. It belongs to the same isotopologue as the fitted J band, so free amplitudes
  cannot help; candidates: a faster 2J-band decay limited by the rate families, field sensitivity of the 2J
  lines, or the early switching transient (compare crop 0.3 s for these slow-decaying spectra).
- diethylamine: spurious model lines near 200 Hz.

## Slow N-H exchange and the 15N isotopologue (v14, v15; owner's suggestion)

With fast exchange the N-H protons are decoupled and the 15N isotopologue has only small couplings to
carbon-bound protons, so its lines fall below the fitted range; with slow exchange the N-H protons stay in the
spin system and both the 13C and the 15N isotopologues gain lines (1J(N,H)). Triethylamine has no N-H: neither
applies. Command: `python scripts/run_batch.py configs/batch_2026-10-10_v14.json --parallel 2 --workers 2 --tag
_v14` (exchange slow, static; start from the chosen pass of each entry). v15 (one fitted exchange rate k,
`nh_exchange` 50 1/s) did not finish one evaluation of stage 1 in 46 min (exact Liouville model with the
component search) and was stopped.

| entry | previous (same data) | v14 objective | 1J(N,H) (Hz) | verdict |
|---|---|---|---|---|
| tert-butylamine | v1 s2 0.121 | 0.110 | -75.76 +- 0.05 | **v14**, 2J band slightly better |
| isopropylamine | v2 s2 0.063 | 0.116 (stage 2 stopped) | -75.80 (stage 1) | worse, keep v2 |
| ethylenediamine | v1 s2 0.062 | **0.033** | -65.00 +- 0.31 | **v14**: the small lines at 190-193 Hz now fitted |
| diethylamine | v13 s2 0.855 | **0.535** | -75.89 +- 0.15 | **v14**: 251 / 255 Hz lines better, spurious 200 Hz lines weaker; 2J(C2,HN1) -11.6 Hz implausibly large |

Figure runs/figures/batch/amines_v14_compare.png; the overlay draws v14 for these three.
Isopropylamine and triethylamine keep a too-weak 2J band: not explained by N-H coupling.

## Residual field (v16) and non-sudden switching (owner's questions)

Stage 2 already has one decay rate per model line (`--family-edges lines`); no rate of the 2J-band lines of
isopropylamine or triethylamine is at a bound, so the weak 2J band is an area, not a width problem. A 13C
preparation weight does not change the 2J / J band area ratio (1.25 for every weight).

Residual field: `configs/batch_2026-10-10_v16.json` (`field: true`, start from the chosen pass).

| entry | previous | v16 objective | field transverse / z (nT) | verdict |
|---|---|---|---|---|
| isopropylamine | v2 s2 0.063 | **0.034** | 9.0 / 49.9 | **v16**: J band much better; 2J band still too weak |
| triethylamine | v1 s2 0.090 | **0.033** | 6.4 / 45.6 | **v16**: J band better; 2J band still too weak |
| tert-butylamine | v14 s2 0.110 | 0.109 | 26.8 / 40.7 | same, keep v14 |

The three samples give a consistent residual field of about 45-50 nT (mostly along z): this set was measured
in a residual field; the other amines should be refitted with `field: true` too.

Non-sudden switch (feasibility, not in the code): rho0 = sum gamma I_z evolved under H_J + Zeeman(B(t) z),
B(t) = 100 uT exp(-t / tau), then the zero-field amplitudes. 2J / J band ratio: tau 0 1.25, 0.1 ms 1.23,
0.3 ms 1.13, 1 ms 0.67 (signal 40 %), 3 ms 2.8-3.8 (signal 3 %), nearly identical for isopropylamine,
triethylamine and sec-butylamine. It lowers the ratio where the data need it higher, and it would act on
sec-butylamine alike, whose 2J band fits with the sudden drop: rejected, no change of the protocol.
Remaining candidate: a second species (triethylamine's 2J lines at 251 / 255 Hz coincide with diethylamine's,
a common impurity of triethylamine); test with a mixture model.

## Correction: the weak 2J band was a display artefact (later the same day)

Asked by the owner, a stage 2 with rates <= 3 1/s (runs/processed/batch_<id>_v16_rate3_s2, stopped) was worse
(isopropylamine 0.0375, triethylamine 0.0389). Its zoom on 225-265 Hz showed that the v16 fit already matched
the 2J band. Raw heights at 248-252 Hz (no AsLS): triethylamine data 0.0319, model 0.0334; isopropylamine data
0.0087, model 0.0084. In the full-range overlay the simulation's automatic AsLS (smoothness 0.25 Hz) took two
thirds of the model band (0.0116), because the choice only checked the five highest peaks (all in the J
band). scripts/overlay_panels.py auto_baseline now scores the height loss of every line cluster; with it both
2J bands agree in the overlay. The open point "2J band too weak" and the impurity hypothesis are withdrawn;
the non-sudden-switching test above stays valid as a test but answered a question that did not exist.

Conditional numerical results of the given structures.

Commit: (this commit)

## Why diethylamine is poor (owner's question)

- The data: SNR of the z5 average 17 (8001 scans), against 37-220 for the other amines with 2800-6700 scans.
  The raw measurement had 10000 scans; the collaborator's old tool kept 8001 (selection_info.txt: median
  residual 2.9e12, against 0.8-3.6e10 for the other measurements).
- The magnetometer output sits at a large, drifting offset (-17000 to -28000 counts) and reaches the ADC rail
  (-28590): 4543 of 8001 scans have a flat (clipped) segment, 3977 of them within 0.1-4.1 s (scan 0 from
  0.05 s, scan 1 from 8 s). New `scripts/clipped_scans.py` writes the keep list; `average_scans.py --keep`
  accepts its JSON.
- Averaging only the 4024 scans free of clipping in 0.1-4.1 s (~/research/zulf/data/processed/
  2025-10-10-diethylamine/unclipped-4s/) is worse: SNR 11.9 (window 0.3 1/s) and 6.6 (1.5 1/s), against 17.1 and
  13.9. Clipped segments are constant (no noise after the drift filter), and the clipped scans still carry the
  early signal. The noise of the unclipped scans (late noise 77) equals that of the other amines (76-78), so the
  signal per scan is about four times weaker.
- Likely cause: the sensor worked far from its zero (offset near the rail), with reduced, nonlinear
  sensitivity. A measurement problem; to be raised with the collaborator (re-measure with the sensor zeroed).

## Scan selection from the original 10000 diethylamine scans (2026-10-12)

The original run (Google Drive "Raw Data"/diethylamine_9999.zip, 10000 scans, 2025-10-10) is now stored at
`original/2025-10-10-diethylamine/diethylamine/` (20001 files, checked). Per-scan statistics
(~/research/zulf/data/processed/2025-10-10-diethylamine/original-10000/scan_stats.json): pinned fraction in
0.1-1.1 s, template amplitude in 110-270 Hz, band rms, robust z.

| scans | kept by the old tool | dropped by the old tool |
|---|---|---|
| pinned at the rail (> 90 %) | 3101, amplitude 0.0 +- 0.2 | 398 (a few with very large band power) |
| partly pinned | 869, amplitude 22.8 +- 6.8 | 167 |
| never pinned | 4031, amplitude 5.4 +- 0.5 | 1434, amplitude 6.6 +- 1.2 |

The old broadband residual metric is dominated by the magnetometer offset: it kept the pinned (signal-free)
scans and dropped good ones; it did remove a few scans with huge band power, which spoil a plain mean.

Line height over noise (window 1.5 1/s, record 4 s) at 124.9 / 251.1 / 255.6 Hz:

| average | scans | heights |
|---|---|---|
| all 10000 | 10000 | 3.2 / 1.8 / 0.9 |
| old selection (raw/, plain mean) | 8001 | 14.1 / 6.7 / 4.0 |
| old selection, z5 | 7532 | 7.7 / 4.0 / 2.2 |
| never pinned | 5465 | 3.3 / 2.6 / 1.2 |
| not pinned (> 90 %) | 6501 | 10.8 / 4.1 / 2.1 |
| not pinned, band z <= 3 / 5 / 8 / 15 | 5693 / 6011 / 6199 / 6313 | 10.0 / 7.1 / 2.8; **15.8 / 9.6 / 3.7**; 15.2 / 7.8 / 3.4; 4.5 / 1.3 / 0.9 |

New `scripts/select_scans.py` (pinned scans out, band 110-270 Hz robust z <= 5): 6011 scans, the best of these.
A clip-aware weighted mean was tried and dropped (lower line heights). Fit on it: configs v18 (running).
Open: why partly pinned scans carry a larger amplitude; whether select_scans also helps the other runs.
- Owner's check: band outliers out only, pinned scans kept (`select_scans.py --keep-pinned`, 9404 scans; 596
  band outliers, 106 of them pinned scans with very large band power): 3.5 / 2.9 / 1.2, far worse. Pinned scans
  are not harmless: many leave and return to the rail (pinned 90-100 %), and each step adds broadband power while
  they carry no signal. The pinned-scan cut stays.
- v17 (old 8001-scan plain mean, settings of v14): objective 0.692, 1J(C2,H) 146.9 Hz (far from about 133):
  not used.

## Fits on the new selection (v18, v19)

| version | data | objective | 1J C1 / C2 (Hz) | figure (runs/figures/batch/diethylamine_v18_compare.png) |
|---|---|---|---|---|
| v14 | z5 of the old selection (7532) | 0.535 | 124.81 / 132.29 | 125, 251, 255 Hz fitted; quiet near 200 Hz |
| v18 | select_scans band z5 (6011 of 10000) | 0.273 | 125.56 / 135.04 | strong model lines at 200-225 Hz (XA2 band) not in the data |
| v19 | as v18, free isotopologue amplitudes | 0.230 | 124.13 / 128.39 | worse: model lines at 190-205 Hz several times the data |

Lower objectives with worse pictures: with signal weighting (weight 0.2 away from the data cores) and this
noise level, model lines in regions without data lines cost little while the fit follows noise elsewhere. The
objective is not a reliable yardstick here; v14 stays in the overlay. A carbon weakened by 14N is not
supported (free amplitudes made the XA2 band stronger). Next: fit the new average with equal weighting or a
penalty for model lines absent from the data.
- v22 (new selection, `band_weighting: equal`, new option `--band-weighting`): objective 0.832 (not comparable with
  signal weighting), 1J 124.48 / 129.88 Hz; no spurious model lines (the XA2 region is quiet as in the data), J
  band and the 251 Hz line placed; 251 Hz at about 60 % of the data height, 256 Hz weak
  (runs/figures/batch/diethylamine_v22_compare.png). Next: v23 = v22 with a fitted residual field.

## Residual field along z only (v24, 2026-10-12)

Fits with a transverse-and-z field (v20 DMEDA / sec-butylamine, v21 sec-butylamine slow exchange, v23
diethylamine) took 2-5 h per stage: a transverse field breaks the total-M blocks of the eigen-decomposition.
Stopped by PID at 06:45 with: sec-butylamine v20 stage 2 0.092 after 3.7 h (v13 0.079, worse; field 1 / 45 nT),
v21 stage 1 0.108 after 5.4 h, DMEDA v20 stage 1 0.141 after 2.5 h (v13 stage 1 0.187), diethylamine v23
stage 1 0.949 after 3.3 h (v22 0.967). The fitted fields were mostly along z, so v24 fits a z field only
(`field_axes: z`, `field_start: 0,0.045`; new run_batch key) for DMEDA, sec-butylamine, ethylenediamine,
dipropylamine and diethylamine (v22 settings).
- v24 (z field only) was still about 20 times slower than zero field (ethylenediamine stage 1 78 min against 4)
  and did not help: ethylenediamine stage 2 0.071 at 33 min (v14 0.033; field 15 nT), sec-butylamine stage 1
  0.106 after 83 min (v13 0.097), DMEDA still in the component search after 2.4 h. Stopped by PID at 09:10.
  The 45-50 nT field appears only in the runs of early October (isopropylamine 09-30, tert-butylamine 10-03,
  triethylamine); ethylenediamine (09-27) gave 15 nT, sec-butylamine is from 2026-04.
- v25: diethylamine (10-10) with the v22 settings in a fixed z field of 46 nT (new run_batch key
  `fixed_field`; no field parameter).

## Diethylamine left out; feature labels (owner, 2026-10-12)

- Diethylamine is marked "experimental spectrum too poor" (magnetometer offset at the ADC rail, half of the scans
  pinned, weak signal) and left out of the amines overlay (scripts/overlay_panels.py EXCLUDED); data and fits are
  kept. v25 (fixed 46 nT field, equal weighting) was stopped at stage 2 (0.894 after 11 min).
- Every overlay panel now lists the non-default choices of its fit, read from the fit's command, result and
  series (`fit_features`): slow N-H exchange, residual field (fitted or fixed, transverse / z), free isotopologue
  amplitudes, background off, equal weighting, window.
- sec-butylamine 115-140 Hz (C2 XA band and the methyl J bands) poorly fitted: v26 slow N-H exchange at zero field
  (window 1.5 1/s), v27 the same with the default window 0.3 1/s (narrow lines below 1 Hz spacing).
- Owner: some small lines are in the model but merged into the big neighbour. Cause: the stage-2 `lines` families
  chain (each line within 0.25 Hz joins its neighbour's group), so DMEDA 131.25-131.98 Hz (30 lines) shared one
  rate, 132.04-133.98 Hz another, and each sec-butylamine isotopologue had one rate over 123.5-126.5 Hz (C4: FWHM
  1.9 Hz, wider than the 0.7 Hz splitting at 124.7 / 125.4 Hz). New `--line-group-span` (default 0 = unchanged;
  run_batch key `line_group_span`) caps a group's span; test tests/test_joint_series.py
  test_line_groups_stop_chaining_at_the_span_limit. Stage-2 refits from v13 stage 1 with span 0.4 Hz:
  runs/processed/batch_{n-n-dimethylethylenediamine,sec-butylamine}_v13span_s2.
- Pitfall: a fit's recorded command (monitor/status.json argv) holds `--family-edges` as the numbers actually
  used, not the keyword `lines`. The v13span stage-2 refits copied that command, so `--line-group-span` had no
  effect (same families as v13): that comparison is void. Copy a command with `--family-edges lines` restored.
- New options `--line-min-distance` (0.25), `--line-min-relative` (0.05), `--edge-merge` (0.3) (defaults unchanged).
- Owner: one rate per line, couplings fixed. v13fixJ: stage 2 from the v13 result, all couplings held at v13
  (`--free-couplings '^NONE$'`), `--family-edges lines --line-min-distance 0.02 --edge-merge 0.02
  --line-min-relative 0.02`: 475 (DMEDA) / 504 (sec-butylamine) families.
- v13fixJ (couplings fixed, one rate per line; 1901 / 2521 parameters, 5-30 min per step): after 2 steps
  objective 0.0659 / 0.0767 (v13 0.074 / 0.079), the merged doublets unchanged (runs/figures/batch/fixJ_compare.png).
  Check: the same couplings with every rate at 1.5 1/s (runs/figures/batch/narrow_rate_check.png). sec-butylamine
  gets a peak at 124.6 Hz but also structure absent from the data (125.3, 127.6 Hz); DMEDA has no line at the
  130.9 Hz data peak at all. The merged peaks are a line-position (small-coupling) problem, not a width problem.
  Stopped (with v27, window 0.3 1/s, stage 2 0.129, worse). Next: local scan of the small couplings at 115-140 Hz.

## Small-coupling scan (2026-10-12)

New `scripts/coupling_scan.py`: 1D scans (J0 +- 8 Hz, step 0.25) of every coupling below 20 Hz, then 2D grids
(step 0.5) of the three that lower the window residual most; rates, phase and field held at the fit, gains solved
again at every point. Results runs/processed/coupling_scans/{dmeda,secbutyl}_v13.json.

- DMEDA v13, window 122-142 Hz (residual 0.00141, objective 0.0738): best single change J(C3,HN2) 1.87 -> -2.63 Hz,
  window -1.4 %, objective 0.200; best pair -1.4 %. Every other coupling: < 0.5 %.
- sec-butylamine v13, window 114-140 Hz (0.00177, 0.0781): best single J(C2,HC4) 6.11 -> 3.11 Hz, window -3.7 %,
  objective 0.135; best pair -4.9 % with objective 1.77.
- No small coupling, alone or in pairs, separates the merged doublets without spoiling the rest: the v13
  couplings are the best of their neighbourhood. The misfit points to a missing ingredient (slow N-H exchange for
  sec-butylamine, v26 running; or a component the model lacks), not to wrong small couplings.
