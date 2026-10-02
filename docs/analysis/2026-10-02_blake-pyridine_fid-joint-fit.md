# Blake pyridine series: joint fit of the FID-processed spectra (2026-10-02)

- Data: raw FIDs `pyridine_{02,10,33,50,66,75,100}` (CSV from the user, 65515
  points at 3999.94 Hz, 16.4 s; uploads, not committed); x = pyridine mole
  fraction 0.02, 0.10, 0.33, 0.50, 0.66, 0.75, 1.00 (seven points; x 0.10 is
  new).
- Question: with one processing, one phase rule and one baseline rule for all
  spectra (2026-10-01_blake-pyridine_data-consistency.md showed that the
  x 0.66 band deficit came from the earlier processed spectra), how do the
  joint-fit couplings and the small-coupling jumps between x 0.66 and 0.75
  change?
- Code: zulf_processing (process_record via Acquisition, asls_baseline),
  scripts/fit_joint_series.py (`--peak-min-sigma` new, commit c028bac);
  scratchpad full_phase.py, sym_phase.py, sym_plot.py.

## Spectra

Crop from sample 210 (52.5 ms, after the 497 Hz switching ringing) to the
end, SG 201/2 on the FID, no apodization, zero fill 1 (grid 0.0613 Hz).
First-order phase from each FID's own switching edge (3.18-3.25 ms), zero-order
phase 179 deg (peak symmetry 174 deg plus the 5 deg bias of that criterion on
a simulated spectrum). Real part, AsLS baseline (smooth 1.5 Hz, p 0.01),
130-210 Hz. Figure:
runs/processed/joint_peakpen_smooth_ms48/sym_phased_vs_blake.png.

## Fit

    python scripts/fit_joint_series.py --series series_fid.json \
      --structure '{"motif": "pyridine ring", "one_bond": {"A2": 178, "A3": 162.5, "A4": 161.5}}' \
      --couplings "$(cat pyr_lit.json)" --prior-sigma-hh 0.5 --prior-sigma-ch 1.5 --prior-weight 10 \
      --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 4 \
      --from-joint runs/processed/joint_peakpen_smooth_ms48/fit.json --seeds seeds_fid.json \
      --starts 40 --spread 1.5 --prior-starts 8 --seed 31 --workers 4 \
      --out runs/processed/joint_fid_sym_ms48

Peak tops for the missing-peak rows: prominence >= 12 % of the maximum and
>= 4 noise sigma (2 sigma gave 49 tops on the noisy x 0.02 spectrum, mostly
noise; 4 sigma: 7, 7, 8, 9, 8, 7, 6). Seeds (9): the 6 best solutions of the
penalty run on the processed spectra and the 48-start best, 48-more best and
48-more data-only best; x 0.10 takes the seed values of x 0.02 (nearest
node). Started 2026-10-02 (one BLAS thread per worker).

Result (48 starts in 2105 s): best total 0.4490; relative data-core residuals
0.352 / 0.266 / 0.235 / 0.237 / 0.227 / 0.207 / 0.159 (x 0.02 ... 1.00).
Superseded before analysis: these spectra start at 52.5 ms, so every line
carries crop wings (first lobe about -0.12 of its height), and the positive-
peak AsLS of this input flattened them only partly, while the model had ideal
lines (record=None). See 2026-10-01_blake-pyridine_data-consistency.md
(crop wings, baseline). Not used further.

## Second run: crop 0.1 s, re-done baseline

Input (user choice): crop 0.1 s to the end, 0.3 1/s, phase0 176.9 deg at the
switching edge, AsLS 1.0 Hz p 0.01 (one setting for the series, chosen by a
scan against Blake's spectra; figure rebaseline_asls10_vs_blake.png), 125-215
Hz on a 0.0307 Hz grid, ideal-line model (no record). Peak tops for the
missing-peak rows at 2 sigma: 11, 12, 11, 10, 8, 9, 7.

    python scripts/fit_joint_series.py --series series_rb.json ... (as above) \
      --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 2 \
      --from-joint runs/processed/joint_peakpen_smooth_ms48/fit.json --seeds seeds_fid.json \
      --starts 40 --spread 1.5 --prior-starts 8 --seed 37 --workers 4 --out runs/processed/joint_fid_rb_ms48

Started 2026-10-02 03:20 UTC; 48 starts in 3564 s. Best total 0.4225, next
0.466 and 0.468 (only one solution within 3 %: not converged). Relative
data-core residuals 0.272 / 0.246 / 0.248 / 0.230 / 0.200 / 0.188 / 0.211
(x 0.02 ... 1.00). The small peak at main + 1.5 Hz is reproduced
(spectra.png); the broad rise at 177-182 Hz is under-fit.

Couplings against the solution on Blake's processed spectra (penalty run,
2026-10-01_blake-pyridine_peak-penalty.md; figure
runs/processed/joint_fid_rb_ms48/J_trends_vs_blake.png; scratchpad
jump_tab.py, jcmp_rb.py):

| | this run | Blake spectra |
|---|---|---|
| largest step > 1 Hz at 0.66 -> 0.75 | 0 couplings | 7 |
| largest step > 1 Hz at 0.50 -> 0.66 | 5 | 0 |
| largest step > 1 Hz at 0.75 -> 1.00 | 9 | 3 |
| 1J(C2,H2) x 0.02 -> 1.00 | 178.78 -> 177.03 | 178.79 -> 177.32 |
| 1J(C3,H3) | 164.12 -> 162.60 | 164.22 -> 162.73 |
| 1J(C4,H4) | 158.84 -> 152.36 | 162.09 -> 156.93 |

The jump of the small couplings did not go away; it moved one step down in x
(0.50 -> 0.66) with the same pattern: J(A2,HA3) 4.5 -> -3.0, J(A2,HA5)
-1.8 -> 2.4, J(HA3,HA5) -0.5 -> 1.5, J(A4,HA3) -0.3 -> 1.6 Hz. Steps of
several hertz between neighbouring points, at a place that depends on the
processing, look like a switch between two near-equivalent solution branches
of the joint fit rather than a property of the samples. 1J(C2,H2) and
1J(C3,H3) agree within 0.3 Hz; 1J(C4,H4) is about 3.5 Hz lower throughout.

The chained third run did not start: its seeds mixed 7-node solutions of
this run with 6-node reference solutions (IndexError). The three reference
seeds were interpolated onto the 7 nodes (np.interp in x; x 0.10 between
0.02 and 0.33) and the third run started by hand at 04:21 UTC.

## Third run: no baseline, model through the same processing

User decision: no AsLS; the model is rendered with the same crop, window, SG
and phase (series entries carry `record` and `phasing`). Input
series_fid03_4.json (now the 8.1 s set): crop 0.1-8.1 s (user: keep at least 8.1 s),
0.3 1/s, zero fill 3 (0.0417 Hz grid, like Blake's 0.0452 Hz),
phase0 176.9 deg at the switching edge, no baseline. Seeds: the 6 best of the
second run and the three reference solutions; same fit settings, 48 starts,
seed 41, out runs/processed/joint_fid_route2_ms48 (chain script starts it
when the second run ends). Figures of this run show data and model after
subtracting the same data-derived AsLS baseline (display only; the residual
is unchanged).

Record length vs points (single-thread timing): 4.1 s with zero fill 6 gives
1440 points in 140-200 Hz and 2.2 s per jacobian; the full 16.4 s record with
zero fill 1 gives 976 points and 11.5 s. Line amplitudes after 4 s are below
e^-8 of their start (rates 2-3 1/s), so the longer record adds noise only.
Chosen by the user: 0.1-8.1 s, zero fill 3: 1440 points, 4.3 s per jacobian.

Complex fit (user question; no AsLS means every processing step is linear and
the model goes through the same steps): `--real-only false` on
fit_joint_series / reliability_series (default true; user: a parameter, not a
series-file key) reads the values as complex. Input
series_fid03_8zf3c.json: the same 0.1-8.1 s spectra kept complex after the
176.9 deg phase (real part = absorption, used for the peak tops); the shared-
phase gain refines a residual zero-order phase per spectrum, the delay is free
near the edge. Timing single thread: jacobian 4.6 s (real-only 4.9 s);
data-core residuals with the processed-spectrum couplings 0.27-0.39 for both.
Peak tops at 2 sigma: 19 on x 0.02 (noise), so 3 sigma: 11, 14, 12, 12, 11,
11, 7. The third run uses this input (out runs/processed/joint_fid_route2c_ms48,
`--real-only false --peak-min-sigma 3`, seed 41, `--from-joint` the second
run).

Result: 48 starts in 11371 s (started 04:21 UTC). Best hard total 0.2956
(start 3, the fourth-best solution of the second run as seed), next 0.3006
(+1.7 %) and 0.3041 (+2.9 %): three solutions within 3 %. Relative data-core
residuals (complex) 0.236 / 0.238 / 0.213 / 0.210 / 0.207 / 0.207 / 0.151.
The script's spectra.png shows magnitudes for complex data; the real parts
with the data-derived AsLS baseline (1.0 Hz, p 0.01) subtracted from data and
model alike (display only) are in spectra_baselined.png (scratchpad
replot_route2c.py). The 177-182 Hz region and the small peak are reproduced;
the largest remaining misfit is 148-152 Hz at x 1.00.

Couplings (figure runs/processed/joint_fid_route2c_ms48/J_trends_compare.png,
three runs; scratchpad jump_c.py):

| | complex, no AsLS | real, AsLS | Blake spectra |
|---|---|---|---|
| largest step > 1 Hz at 0.50 -> 0.66 | 8 | 5 | 0 |
| at 0.66 -> 0.75 | 3 | 0 | 7 |
| at 0.75 -> 1.00 | 3 | 9 | 3 |
| 1J(C2,H2) x 0.02 -> 1.00 | 179.13 -> 177.06 | 178.78 -> 177.03 | 178.79 -> 177.32 |
| 1J(C3,H3) | 164.22 -> 162.84 | 164.12 -> 162.60 | 164.22 -> 162.73 |
| 1J(C4,H4) | 162.04 -> 156.88 | 158.84 -> 152.36 | 162.09 -> 156.93 |

1J(C2,H2) and 1J(C3,H3) agree within 0.4 Hz in all three fits. 1J(C4,H4) of
the complex fit equals the Blake-spectra fit within 0.1 Hz; the AsLS fit is
3-4.5 Hz lower, so that value (and the small couplings, which differ from the
complex fit by up to 8.7 Hz, e.g. J(A4,HA2), J(A2,HA3)) depends on the AsLS
route and is not used. The small-coupling jump is still there, spread over
the two steps 0.50 -> 0.66 -> 0.75: J(A2,HA3) 5.02 -> 5.02 -> -2.19,
J(A2,HA5) -3.34 -> -1.94 -> 1.87, J(A4,HA3) -1.61 -> 2.17 -> 2.17,
J(HA3,HA4) 9.49 -> 6.48 -> 6.08 Hz; 1J(C4,H4) also drops 2.8 Hz at
0.50 -> 0.66.

Model delay per spectrum (free, +-10 ms): -4.11, -2.18, -3.91, -3.54, -0.85,
+2.58, -0.13 ms (x 0.02 ... 1.00). The phasing put the switching edge at
-3.22 ms in this convention, which the four dilute spectra match; from x 0.66
on the delay moves 3-6 ms away, at the same two steps as the coupling jump. With
the couplings held, the residual is 2-6 % higher per ms away from the fitted
delay (scratchpad delay_scan.py), so the delay is determined, but it is not
the same across spectra although every FID was referenced to its own edge: a
free delay can trade phase slope against line positions.

## Fourth run: one model delay for the series

User choice: one shared delay instead of a free delay per spectrum. New option
`--shared phase_delay` (fit_joint_series, reliability_series): the listed
spectrum parameters are one value for the series, placed after the coupling
blocks in z; start = mean of the per-spectrum starts (here the third run's
delays, mean -1.71 ms); Jacobian column = sum over spectra (test against
finite differences in tests/test_joint_series.py).

    python scripts/fit_joint_series.py --series series_fid03_8zf3c.json --real-only false \
      --shared phase_delay ... (as the third run) --peak-min-sigma 3 \
      --from-joint runs/processed/joint_fid_route2c_ms48/fit.json --seeds seeds_shared.json \
      --starts 40 --spread 1.5 --prior-starts 8 --seed 43 --workers 4 --out runs/processed/joint_fid_shdelay_ms48

Seeds (9): the 6 best solutions of the third run and the three reference
solutions. Started 2026-10-02 08:20 UTC.

Result: 48 starts in 11679 s. Best hard total 0.2965 (start 2), next 0.3045
(+2.7 %), 0.3116. Shared delay -1.86 ms (start -1.71 ms; switching edge
-3.22 ms). Relative data-core residuals 0.241 / 0.240 / 0.217 / 0.212 /
0.210 / 0.212 / 0.136: against the per-spectrum delays +0.9 to +2.1 % on
x 0.02-0.75 and -9.8 % on x 1.00 (a different minimum there). Figures:
runs/processed/joint_fid_shdelay_ms48/spectra_baselined.png (display
baseline as above), J_trends_compare.png (shared delay, per-spectrum delay,
Blake spectra).

| | shared delay | delay per spectrum |
|---|---|---|
| largest step > 1 Hz at 0.50 -> 0.66 | 7 | 8 |
| at 0.66 -> 0.75 | 5 | 3 |
| at 0.75 -> 1.00 | 0 | 3 |
| 1J(C2,H2) x 0.02 -> 1.00 | 179.01 -> 177.60 | 179.13 -> 177.06 |
| 1J(C3,H3) | 164.23 -> 162.76 | 164.22 -> 162.84 |
| 1J(C4,H4) | 162.06 -> 157.29 | 162.04 -> 156.88 |

Couplings agree with the per-spectrum-delay fit within 1 Hz except J(A2,HA6)
and J(HA2,HA6) at x 1.00 (5.8 and 5.5 Hz; the other minimum there) and
J(A4,HA3) (2.8 Hz). The small-coupling jump at 0.50 -> 0.75 stays with one
delay for all spectra (J(A2,HA3) 4.78 -> -2.99, J(A2,HA5) -2.10 -> 1.96,
J(HA2,HA4) 2.56 -> 5.01, J(HA3,HA4) 9.45 -> 6.37 Hz): the per-spectrum delay
jump of the third run was a consequence, not the cause.

Choice for the fifth run (rule stated before the result, scratchpad
choose_delay.py): shared delay if its total <= 1.05 x the per-spectrum one, no
spectrum residual up by more than 10 %, delay within 1.5 ms of the edge.
Ratio 1.003, largest rise 2.1 %, delay 1.36 ms from the edge: shared delay.

## Fifth run: no monotone constraint

User question: drop the monotone constraint. New option `--shape free`
(fit_joint_series): every coupling has its own value at every concentration
(z block = J at the nodes); the spectra are then tied only through the priors
on the series averages (and `--shared`, not used here). Start perturbations
shift all nodes of a coupling by one offset; prior-drawn starts are linear
in x. Test: values = parameters, Jacobian against finite differences.
Settings as the fourth run (user: take the better delay model of runs three
and four; chosen automatically by the rule above: shared delay), so the only
change against the fourth run is the constraint. Seeds: seeds_free.json (the
fourth run's 6 best and the three references), `--from-joint` the fourth run,
seed 47, out runs/processed/joint_fid_free_ms48 (scratchpad chain_free2.sh);
started 2026-10-02 11:35 UTC. Smoke test (one start from the third run's
best, 2 evaluations): total 0.2792 against 0.2956, J(A2,HA3) 6.60, 6.64,
5.06, 4.82, 5.26, -2.25, -2.86 Hz (x 0.02 ... 1.00).

Result: 48 starts in 11396 s. Best hard total 0.2361 (-20 % against the
monotone fit, 0.2965), next 0.2443, 0.2452. Shared delay -2.57 ms. Relative
data-core residuals 0.217 / 0.194 / 0.198 / 0.190 / 0.179 / 0.205 / 0.135
(1-19 % below the monotone fit). Figures runs/processed/joint_fid_free_ms48/
spectra_baselined.png, J_trends_compare.png (free, monotone, Blake spectra).
15 of 19 couplings go up and down by more than 0.3 Hz (e.g. 1J(C2,H2) 177.91,
178.58, 177.46 Hz at x 0.66, 0.75, 1.00), but the jump stays: J(A2,HA3)
6.00 -> -4.00 (0.66 -> 0.75), J(HA3,HA4) 11.56 -> 5.92, J(A4,HA3)
-3.50 -> 2.85, 1J(C4,H4) 160.10 -> 157.32 Hz (0.50 -> 0.66).

### Two-branch test (scratchpad branch_test.py, branch_test.json)

Every spectrum alone (couplings and rates free, delay held at -2.57 ms, no
priors, no peak rows), started once from the free fit's couplings at x 0.50
(branch A) and once from those at x 0.75 (branch B):

| spectrum | A rel. residual | B rel. residual | B cost / A cost - 1 |
|---|---|---|---|
| x 0.02 | 0.224 | 0.233 | +6.9 % |
| x 0.10 | 0.213 | 0.267 | +68 % |
| x 0.33 | 0.197 | 0.236 | +47 % |
| x 0.50 | 0.189 | 0.226 | +47 % |
| x 0.66 | 0.139 | 0.199 | +97 % |
| x 0.75 | 0.135 | 0.205 | +117 % |
| x 1.00 | 0.113 | 0.128 | +14 % |

Branch A is lower for every spectrum, the concentrated ones included (x 0.75:
0.135 against 0.205 in the joint fits). In branch A the couplings change
smoothly with x, 1J(C4,H4) 162.0 -> 158.7 Hz and 1J(C2,H2) 179.2 -> 176.5 Hz,
no jump. The joint fits do not reach it because of the priors on the series
averages: branch A puts J(HA3,HA4) at 11.2-13.3 Hz (prior 7.7 +- 0.5: 9 sigma),
J(HA2,HA5) 3.4 (prior 0.9: 4.9 sigma), J(HA3,HA5) -0.7 (prior 1.4: -4.1
sigma); branch B lies on the other side. Splitting the series (dilute in A,
concentrated in B) brings the averages to the priors (J(HA3,HA4) mean 8.69 Hz)
at the cost of the concentrated spectra. The jump is made by the series-average
priors, not by the samples, the processing, the delay or the monotone
constraint.

Branch A is far from the literature small couplings (J(HA3,HA4) about 12 Hz
against about 7.7 Hz), so it is not a better assignment either: these spectra
fix some combinations of the small couplings, not the couplings one by one
(as found on Blake's spectra, 2026-09-30_blake-pyridine_reliability.md).

## Conclusion

The complex fit without AsLS fits all seven spectra with residuals 0.15-0.24
and agrees with the Blake-spectra fit on all three 1J; the AsLS route biased
1J(C4,H4) by about -3.5 Hz. The small-coupling jump does not come from the
processing of any single spectrum (its place moves between 0.50 -> 0.66 and
0.66 -> 0.75 with the processing); in the complex fit it coincides with a jump
of the free model delay.

## Reliability classes of the complex monotone fits

reliability_series.py now handles fits with missing-peak rows: their stored
scores contain the hard rows, so the compared score is the data residual plus
those rows (prior removed); the peak settings must be equal in all fits, and
the direct check of each fit's best includes the rows. (The `--spectra` refit
check still compares without them.)

    python scripts/reliability_series.py --series series_fid03_8zf3c.json --real-only false --shared phase_delay \
      --structure ... --couplings "$(cat pyr_lit.json)" --signal-threshold 2.5 --signal-taper 4.0 \
      --fit per-spectrum-delay=runs/processed/joint_fid_route2c_ms48/fit.json \
      --fit shared-delay=runs/processed/joint_fid_shdelay_ms48/fit.json --out runs/processed/reliability_fid_monotone

Best data score 0.2928 (per-spectrum delay) and 0.2929 (shared); 4 solutions
within 3 % (3 + 1). Reliable (spread <= 1.6 Hz): 1J(C3,H3) 0.5, 1J(C2,H2) 1.0,
2J(C3,H4) 1.0, 4J(C3,H6) 1.2, 4J(C2,H5) 1.4 Hz. Trend only: 3J(H2,H3),
3J(C4,H2), 1J(C4,H4) (2.7 Hz), 2J(C2,H3), 4J(H3,H5), 4J(H2,H4), 3J(H3,H4),
4J(H2,H6), 3J(C2,H6), 2J(C4,H3). Not determined: 5J(H2,H5), 2J(C3,H2),
3J(C3,H5), 3J(C2,H4). Figure (user: classes instead of the branch-switch
shading): runs/processed/joint_fid_shdelay_ms48/J_trends_compare_reliability.png
(panel title: n-bond type and ring positions, N = 1; scratchpad
jcmp_annot.py). Caveats: only 4 solutions in the set, all in the split
branches; "reliable" means the near-best solutions agree, not that the value
is right (4J(C2,H5) changes sign at 0.66 -> 0.75 in all of them). The free fit
and the branch-A run are not pooled here.

## Sixth run: monotone, all spectra started in branch A

User choice (option 2 of three): monotone joint fit as the fourth run (shared
delay, same priors and peak rows), every start in branch A. Seeds
(seeds_branchA.json): the branch-A couplings of the two-branch test (one per
spectrum), their isotonic fit (per coupling the better of rising and falling,
pool-adjacent-violators) and the x 0.50 values at every x; 37 perturbed starts
around them (spread 1.5 Hz, smaller than the 5-10 Hz between the branches),
no prior-drawn starts; spectrum parameters from the fifth run.

    python scripts/fit_joint_series.py --series series_fid03_8zf3c.json --real-only false --shared phase_delay \
      ... (as the fourth run) --from-joint runs/processed/joint_fid_free_ms48/fit.json \
      --seeds seeds_branchA.json --starts 40 --spread 1.5 --prior-starts 0 --seed 53 --workers 4 \
      --out runs/processed/joint_fid_branchA_ms48

Question: with the priors unchanged, does the fit stay in branch A (and how
does its total compare with 0.2965), or does it return to the split? Started
2026-10-02 15:57 UTC.

Result: 40 starts in 9749 s. Best total 0.2321 (start 0, the isotonic seed),
next 0.2384 (+2.7 %), 0.2496. Data part 0.2210 + prior 0.0110, against 0.2929
+ 0.0036 for the split solution of the fourth run: data -25 %, total -22 %
(also below the free-shape fit, 0.2361). Shared delay -4.45 ms. Relative
residuals 0.236 / 0.216 / 0.195 / 0.189 / 0.148 / 0.144 / 0.126 (x 0.66 and
0.75: 0.148 and 0.144 against 0.21). The fit stays in branch A: no step above
2 Hz; largest steps J(A2,HA3) 5.24 -> 7.05 and J(HA3,HA4) 11.80 -> 13.28
(0.50 -> 0.66), J(HA3,HA4) 9.82 -> 11.40 and J(A4,HA2) 4.68 -> 6.35
(0.02 -> 0.10). 1J fall smoothly: 1J(C2,H2) 179.22 -> 176.53, 1J(C3,H3)
164.23 -> 161.76, 1J(C4,H4) 161.96 -> 158.57 Hz. The earlier fits never
reached branch A from their starts; the priors did not keep the fit from it.

Reliability, three complex monotone fits pooled (runs/processed/
reliability_fid_all_monotone): only the two best branch-A solutions are within
3 %; 17 couplings reliable (spread 0.1-1.4 Hz), 3J(C2,H6) and 3J(C4,H2) not
determined (1.6, 1.8 Hz). With two solutions in the set this is a statement
about those two, not a budget. Figures runs/processed/joint_fid_branchA_ms48/
spectra_baselined.png, J_trends_compare_reliability.png.

Branch A is the best fit to the data but several proton-proton couplings are
far from the literature values of pyridine, which hardly depend on the
solvent: 3J(H3,H4) 9.8-13.6 Hz (7.7), 5J(H2,H5) 3.2-3.6 (0.9), 4J(H2,H6)
1.2-2.9 (-0.1), 4J(H3,H5) 0.7 to -1.4 (1.4); also 4J(C3,H6) 3.2-4.0 (-1.4).
These spectra are mostly sensitive to the C-H couplings of each 13C
isotopologue; the H-H couplings enter only in combination with them, so the
data can trade H-H against C-H values. Branch A is therefore the best
description of the data within this model, not an assignment of the
individual small couplings.

## Uncertainty budget (queued after the sixth run)

User: the fits look good; the uncertainties must be right. The linearised
errors of fit.json (0.01 Hz for 1J, 0.02-0.06 Hz otherwise) assume white,
independent residuals and one minimum. Measured against that: residual
correlation length 11-20 points (0.45-0.84 Hz; errors x 3.3-4.5); the J of the
three complex fits differ by 4-34 x their linearised errors; the free-shape
fit scatters about a quadratic in x by 0.19 (1J(C3,H3)), 0.38 (1J(C2,H2)),
0.73 (1J(C4,H4)) and 0.2-2 Hz (small couplings). Planned parts, per coupling
and concentration:

1. noise: linearised error times sqrt(correlation length);
2. near-equivalent solutions: spread of the pooled 3 % set of all monotone
   shared-delay fits (reliability_series);
3. processing: local refits (one start from the reference) of six processing
   variants (scratchpad make_variants.py): crop start 320 / 480 samples
   (80 / 120 ms), apodization 0 / 0.6 1/s, record 6 / 10 s; zero-order phase is
   not varied (the shared-phase gain absorbs it);
4. bias: four synthetic series (make_synthetic.py), truth = reference + one
   constant offset per coupling (N(0, 0.3 Hz) 1J, N(0, 0.8 Hz) others),
   spectrum = model(truth) + the real complex residual of the reference
   (identity check: zero offsets give the data exactly), refitted from the
   reference with the same settings and priors.

Reference = lower total of the fourth and sixth runs. Chain
run_uncertainty.sh (10 refits, 4 at a time) starts when the sixth run ends.

## Conclusion (after runs four and five)

The jump of the small couplings between x 0.50 and 0.75 comes from the Gaussian
priors on the series averages: every spectrum alone prefers one branch (A),
which lies several sigma from the literature centres, and the joint fit meets
the priors by moving the concentrated spectra to another branch. Neither the
per-spectrum delay nor the monotone constraint causes it. 1J values depend on
the branch by up to 1.5 Hz at high x.

## Open points

- One shared delay (fourth run): the jump stays; the delay does not cause it.
- Search not converged (three solutions within 3 %).
- Joint fit without the series-average priors (or with per-spectrum priors), and the monotone joint fit
  started from branch A for all spectra.

Commit: see git log (this file).
