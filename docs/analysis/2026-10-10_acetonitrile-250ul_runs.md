# acetonitrile 250 uL, fringe field: the second acquisition and both together (2026-10-10)

- Data (drive only, `/Volumes/Extreme/research/zulf/data/original/2026-09-28-acetonitrile-fringe-field-250ul/`):
  - `ACN_fringe_field_ms1_LF_250uL/`: run 1, 5539 scans, 2026-09-28 16:53 to 2026-10-02 15:54 (received 2026-10-06);
  - `received-2026-10-10/ACN_fringe_field_ms1_LF_250uL/`: the same folder name received again with 13924 scans:
    scans 0-8384 are a new acquisition, 2026-10-03 13:39 to 2026-10-09 13:52 (run 2); scans 8385-13923 are run 1
    renumbered (byte-identical: old n = new n - 8385, checked on scans 0, 1, 2000, 5538 and by the time stamps:
    the only time jump is at 8385). Same sequence file (standard_zf_2000Hz_no_dead.seq), 2000 Hz.
- Question (Xuehan): fit the 250 uL acetonitrile, now with the second acquisition.
- Code: 3bd42cb.

## Setup

Averages (scripts/average_scans.py, from the drive folder; ~/research/zulf/data/processed/
2026-09-28-acetonitrile-fringe-field-250ul/):

| average | scans | kept (--exclude-z 5) |
|---|---|---|
| run 1, z5 (2026-10-06) | 5539 | 4907 |
| run 2, z5 (`run-2026-10-03/z5`, `--keep 0-8384`) | 8385 | 7324 |
| run 2, all scans (`run-2026-10-03/all-scans`) | 8385 | 8385 |
| both, z5 (`combined-13924/z5`) | 13924 | 12242 |

Whole-grid series (make_series_entry.py, defaults: crop 0.1 s, record 8 s, 0.3 1/s, zero fill 3, grid
20-380 Hz minus the mains and instrument lines, `--exclude 81.5,86` as for ethanol):
runs/series/acn250_{run1,run2,combined}_full.

Fits (the model of the best earlier acetonitrile fit acn_rerun_fam2, objective 0.0439 on two bands: chain
C1 (CH3) - C2, 1J(C1,H1) start 136.3 Hz, field fitted from 37 / 45 nT; now on the whole grid with one decay
rate per model line and rate bounds 0.1-40 1/s as for ethanol, D59 amendment):

    python scripts/fit_joint_series.py --series runs/series/acn250_<run>_full/series.json --real-only false \
        --shape free --exchange fast --range 20,380 \
        --structure '{"compound": "acetonitrile", "chain": {"groups": [["C1","C",3],["C2","C",0]], "bonds": [["C1","C2"]]}, "one_bond": {"C1": 136.3}}' \
        --rate-bounds 0.1,40 --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 \
        --peak-min-sigma 3 --starts 8 --workers 3 --max-nfev 300 --fit-field --field-start 0.037,0.045 \
        --family-edges lines --out runs/processed/acn250_<run>_full_fit

## Results

| | run 1 (z5, 4907 scans) | run 2 (z5, 7324 scans) | both (z5, 12242 scans) |
|---|---|---|---|
| 1J(C1,H1) (Hz) | 136.282 +- 0.001 | 136.289 +- 0.001 | 136.286 +- 0.001 |
| B transverse / z (nT) | 38.4 / 43.1 | 39.8 / 45.0 | 39.4 / 43.9 |
| \|B\| (nT) | 57.7 | 60.1 | 59.0 |
| delay (ms) | 3.33 | 3.30 | 3.31 |
| score / residual | 0.233 / 0.490 | 0.214 / 0.466 | 0.165 / 0.405 |
| best starts | 0.2331, 0.2545 x3 | 0.2143, 0.2450 x3 | 0.1649, 0.1917, 0.2012 |

Each fit 8 starts in about 6 min; one of 8 rates at a bound (run 2). The whole-grid score is not comparable with
the two-band score 0.0439 of acn_rerun_fam2 (different fitted region, most of the grid holds only noise).

## Conclusion

Conditional numerical result for CH3-CN with a uniform static field: the new acquisition (2026-10-03 to 10-09)
reproduces run 1: 1J(C,H) within 0.007 Hz (136.28-136.29 Hz) and the field within 2.4 nT in size and a few
nT per component (|B| 58-60 nT, transverse 38-40 nT, z 43-45 nT), the same delay (3.3 ms). The field estimate of
the earlier two-band fits (37 / 45 nT, |B| 59 nT) holds on the whole spectrum. Both acquisitions together:
1J(C,H) 136.286 Hz, |B| 59.0 nT (transverse 39.4, z 43.9 nT).

## Open points

- The field of the rubbing alcohols should be compared with these 59 nT (overnight batch: isopropanol, MDMP;
  ethanol stage 3 gave 63.6 nT, z only).
- Run 2 alone vs run 1: the field size differs by 2.4 nT; a drift between the acquisitions or the fit spread.

Commit: (this commit)
