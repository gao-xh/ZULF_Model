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

(filled in during the night)
