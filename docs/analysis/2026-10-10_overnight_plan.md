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

## Progress

(filled in during the night)
