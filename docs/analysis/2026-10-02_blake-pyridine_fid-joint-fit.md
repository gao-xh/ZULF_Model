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

## Conclusion

The complex fit without AsLS fits all seven spectra with residuals 0.15-0.24
and agrees with the Blake-spectra fit on all three 1J; the AsLS route biased
1J(C4,H4) by about -3.5 Hz. The small-coupling jump does not come from the
processing of any single spectrum (its place moves between 0.50 -> 0.66 and
0.66 -> 0.75 with the processing); in the complex fit it coincides with a jump
of the free model delay.

## Open points

- Fix the model delay at each FID's edge (or one shared delay) and refit, to
  see whether the jump of the small couplings goes with the delay.
- Search not converged (three solutions within 3 %).

Commit: see git log (this file).
