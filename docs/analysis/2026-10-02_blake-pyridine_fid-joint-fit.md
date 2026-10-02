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

Results: pending.

Commit: (pending)
