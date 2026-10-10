# Overnight plan: four new data sets onto the drive, then known-structure fits (2026-10-10)

Owner's request (2026-10-10, night): sort the raw data of four downloaded zips (14N_amines, Urine_Metabolites,
real_world_samples, ACN_fringe_field_ms1_LF_250uL) and the note "Important Data Notes.docx" onto the external
drive, not the Mac; delete the zips; fit the data one by one; accepted the next morning.

Decisions (owner, answered before the run):

- Raw data lives only on the drive `/Volumes/Extreme/research/zulf/data/` (same layout as the Mac:
  `original/<measurement>/`, `raw/<measurement>/`), one tar archive per measurement (exFAT allocates 1 MiB per
  file; about 30 000 scan files per measurement), with a file list and the sha256 of the archive.
- The existing Mac originals (acetonitrile, ethanol, isopropanol) move to the drive the same way and are deleted
  from the Mac after a checked copy (sha256 of the archive and a file-by-file comparison of its listing).
  Processed averages and spectra (small, regenerable) stay on the Mac in `~/research/zulf/data/processed/`.
- The source still has the data, so the zips are deleted after their content is checked against the archives.
- Fits: the known structure of each compound, a quick first pass per compound so every one has a result by
  morning. The compound is named in the file / folder name (owner). Samples in `real_world_samples` and
  measurements whose name says "field" are fitted with a static field; all others at zero field (default).
- Drive copies: `~/code/extreme-drive/scripts/archive_measurement.py` (one uncompressed tar per measurement,
  listing with every file's sha256, the archive read back from the drive and compared file by file).

Notes from "Important Data Notes.docx" (owner's collaborator):

- `14N_amines` also holds the brominated compound methyl 2-bromopropionate.
- `real_world_samples`: the rubbing alcohols in brown polypropylene, unknown field, probably the acetonitrile
  field (about 59 nT); fit the field of both alcohols and compare.
- `real_world_samples/methyl_dimethylphosphonate`: taken in double acquisition mode: divide the frequency axis
  by 2 (process with half the recorded sampling rate); a strong 31P band at about 0-70 Hz, 13C lines weaker.

## Steps

1. Wait for the downloads; list every zip (measurements, scan counts, sampling rates, sequences).
2. Per measurement: name `<YYYY-MM-DD>-<compound>[-details]` (date of the first scan); tar the folder as given
   to `original/<measurement>/<folder>.tar` on the drive, file list `<folder>.files.txt`, sha256 into the drive
   README inventory and a "on the drive" line in `~/research/zulf/data/README.md`; read-only.
3. Duplicates of data already sorted (acetonitrile; ethanol and isopropanol) are compared by file list and sha256
   and kept once.
4. Per measurement, in a temporary folder on the Mac (deleted afterwards): `scripts/average_scans.py` (all scans
   and the z5 selection, halves) into `~/research/zulf/data/processed/<measurement>/`, then
   `scripts/make_series_entry.py` (whole grid) into `runs/series/<measurement>/`.
5. Known structure per compound (motif or chain specification, or a molecule from SMILES), quick fit:
   `fit_joint_series.py` whole grid, `--family-edges lines`, rate bounds 0.1-40 1/s, a few starts; field only
   where the rule says so. One analysis log per data set group in `docs/analysis/`, results table here.
6. Move the old Mac originals to the drive (step 2 checks), delete the Mac copies, delete the zips.
7. Commit and push the logs; the morning summary lists what is done, what failed and what to check.

## Changes during the night (owner)

- The drive was reformatted with 4 KiB clusters (owner, `newfs_exfat -R -b 4096`), so raw data is stored as plain
  folders (`~/code/extreme-drive/scripts/store_measurement.py`: the folder as received + `<folder>.files.txt` with
  every file's sha256, read back from the drive and checked), no tar. `~/AGENTS.md`, the drive README and the
  extreme-drive README say so. Speed of a 2 GiB file unchanged (write 892 / read 885 MB/s; 1 MiB clusters 918 / 888).
- Processed data must be on the drive too (owner): `~/code/extreme-drive/scripts/sync-to-drive.sh` copies
  `~/research` (with `processed/`) at the end.
- Acetonitrile 250 uL first (owner): docs/analysis/2026-10-10_acetonitrile-250ul_runs.md.

## What the downloads held

Safari unpacked the zips itself (no .zip left). 18 new measurements, about 33 GB:

| measurement (drive) | area | scans | rate |
|---|---|---|---|
| 2026-04-12-methyl-2-bromopropionate | original | 8296 | 4000 Hz |
| 2025-11-04-n-n-dimethylethylenediamine | original | 6385 | 4000 Hz |
| 2025-10-10-diethylamine | raw (selected scans) | 8001 | 4000 Hz |
| 2025-11-07-dipropylamine | raw (selected scans) | 5707 | 4000 Hz |
| 2025-09-27-ethylenediamine | raw (selected scans) | 5991 | 4000 Hz |
| 2025-09-30-isopropylamine | raw (selected scans) | 5677 | 4000 Hz |
| 2026-04-06-sec-butylamine | original | 6681 | 4000 Hz |
| 2025-10-03-tert-butylamine | raw (selected scans) | 6713 | 4000 Hz |
| 2026-04-09-tetramethylpiperidine | original | 5741 | 4000 Hz |
| 2025-09-25-triethylamine | original | 2817 | 4000 Hz |
| 2025-12-09-pyridine-dilutions-and-pyridinium | original (02, 10, 33, 66, 100 mol %, pyridinium 6 M) + raw (50, 75 mol %) | 89575 | 4000 Hz |
| 2026-03-16-l-lactic-acid | raw (selected scans) | 7688 | 4000 Hz |
| 2026-03-10-l-alanine | raw (selected scans) | 12124 | 4000 Hz |
| 2026-03-20-glycerol-40wt-h2o | original | 7747 | 4000 Hz |
| 2026-04-01-serine | original | 7711 | 4000 Hz |
| 2026-07-30-methyl-dimethylphosphonate-double-acq | original | 15962 | 8333 Hz (recorded; 4167 Hz real) |
| ethanol, isopropanol (real_world_samples) | duplicates of the drive copies, checked file by file | | |
| acetonitrile 250 uL (received again) | original/2026-09-28-.../received-2026-10-10/: scans 0-8384 a new run (2026-10-03 to 10-09), 8385-13923 the old run renumbered | 13924 | 2000 Hz |

"raw (selected scans)": folders with `data/` (the kept scans, renumbered) and `selection_info.txt` of the old
ScanSelector (e.g. diethylamine 8001 of 10000 kept). The extra files in tetramethylpiperidine are named
`5_4_26_sec_butylamine_*` (probably misnamed spectra of another tool; kept as received).

## Fits

Batch: `scripts/run_batch.py configs/batch_2026-10-10.json --parallel 3 --workers 3` (results table:
docs/analysis/2026-10-10_batch_results.md). Structures: the confirmed specifications of configs/confirmed_samples.json
(triethylamine, ethylenediamine, lactic acid, alanine, pyridine and pyridinium), the others from SMILES
(zulf_hypothesis.molecule; generic 2J/3J, 1J guesses), methyl dimethylphosphonate as a spin system (31P with the
P-CH3 and two O-CH3 protons, 13C at either carbon; literature-range starts). Not fitted: tetramethylpiperidine
(18-19 protons per isotopologue; the simulation stops, needs a reduced model).

## Progress

- 2026-10-10 ~01:00: all measurements except the pyridine series on the drive and checked; Mac originals and
  Downloads folders removed after the check.
- Acetonitrile 250 uL (two acquisitions) fitted first (owner): docs/analysis/2026-10-10_acetonitrile-250ul_runs.md.
- Acetonitrile LF_1 (downloaded during the night) fitted and compared with 250 uL:
  docs/analysis/2026-10-10_acetonitrile-lf1.md. Finding: with crop 0.1 s an early transient after the field
  switch (not removed by the drift filter) biases the fit; with the high SNR of LF_1 the fit imitated it with a
  117 nT field. With crop 0.3 s both fit to the noise (1J 136.28 Hz, |B| 64 / 61 nT; decay times per line in the
  log). The batch keeps crop 0.1 s (the plan); a crop-0.3 check of the high-SNR entries is an open point.
- ~02:00 session restart: the batch, the intake and the LF_1 fit had stopped. Batch restarted (resumable); the
  intake had stopped inside 02_mol_pyridine: the partial drive copy was removed and the intake restarted at 02:46.
- ~03:30: the isopropylamine stage 2 fit failed (broad model, the sharp 133-135 Hz lines missed;
  runs/figures/batch/isopropylamine.png). Cause: the SMILES structures started every sp3 1J at 125 Hz, but CH
  next to N is 133.5 Hz (2026-10-03 analysis); diethylamine's 2J -10 Hz likely the same. New 1J starts from the
  neighbours (D61 amendment, 8e713a7). The batch was stopped (running then: diethylamine s2, sec-butylamine s2,
  dipropylamine s1, all redone) and restarted with `--parallel 2 --workers 2 --skip <the SMILES entries>`; the
  SMILES entries with changed starts run as a second pass: `scripts/run_batch.py configs/batch_2026-10-10_v2.json
  --parallel 2 --workers 2 --tag _v2 --table docs/analysis/2026-10-10_batch_results_v2.md` (fits
  runs/processed/batch_<id>_v2_s1|s2; the first-pass fits stay for comparison). At most 8 fit processes.
- Diagnostic figures of finished entries (experiment, the stage-2 fit's own model, residual):
  runs/figures/batch/<id><tag>.png (scripts/batch_figures.py). The first version silently drew the Studio quick
  look instead of the fit's model when the exact curve had not been computed; that made good fits (alanine) look
  broad. Fixed (it computes fit_model_curve first and skips an entry without it); all figures redrawn.
- ~04:40 assessment from the corrected figures: tert-butylamine, triethylamine, ethylenediamine good; alanine and
  lactic acid follow their lines (scores high from the noise of the data); isopropylamine first pass misses the
  134.8 Hz line (1J start), the second pass fits (stage 2 0.0626; couplings within 0.1-0.2 Hz of the 2026-10-03
  analysis: 1J 133.47 / 124.37, 2J -4.12 / -1.47, 3J 5.23 / 6.00 Hz); diethylamine (second pass 0.559) misses the
  125 Hz line, low SNR with strong drift below 100 Hz.
- Load: at most about 9.5 of 10 cores (batch 3 x 3 workers + one 3-worker fit), no thermal or performance warning
  (pmset); no extra parallel jobs added.
