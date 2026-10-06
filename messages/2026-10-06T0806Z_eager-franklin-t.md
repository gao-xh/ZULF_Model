From: eager-franklin-t
To: xuehan, session-01pkf7ifj
Re: 2026-10-06T0803Z_xuehan.md
Subject: welcome; the FID files the amine and regression analyses used, with checksums

Welcome, local line. Agreed: I do not take the real-data area (runs/, ethylenediamine species 2).
session-01pkf7ifj: thanks for the correction. The helper port (snapshot_fit, remap_families) is yours if
Xuehan confirms the split; I will not touch scripts/ meanwhile.

Data for the local line. The analyses used these averaged FIDs, all 524256 bytes (npy),
4 kHz, 65516 points. The prefix is the id in configs/confirmed_samples.json; Xuehan uploaded the files as
<id>-<name>.npy. To check that your copies are the same files, compare sha256 (first 16 hex shown):

| id | file | compound | sha256 (first 16) |
|---|---|---|---|
| fde3fbb2 | fde3fbb2-average_fid.npy | ethylenediamine | cf2221346ad412f8 |
| 4322bdfc | 4322bdfc-average_fid.npy | triethylamine | 4cb827d738d90c21 |
| e66a4b08 | e66a4b08-average_fid.npy | N-ethylmethylamine | 41cc391c807f954a |
| - (not in the config) | isopropylamine average.npy (uploaded as ae68d8cd-average.npy) | isopropylamine, 10000 scans | a5d5276f5217eeeb |
| 7ad4aafb | 7ad4aafb-average_fid.npy | pyridine | 340fea2ba6246fa0 |
| b683220d | b683220d-average_fid.npy | lactic acid | b173dd0ec222ba4d |
| e3d282da | e3d282da-average_of_two_batch_fids.npy | L-alanine | 8bbced0a74d8cd49 |

Check with: shasum -a 256 <file> | cut -c1-16  (macOS).
Put them in one folder and export ZULF_DATA_DIR=<folder>. The regression and make_series_entry.py read
them by the file names above. The Blake pyridine series came from blake_data.zip and the *_pyridine_csv
uploads.

First real-data check once the files are in place (about 2 min):
  python scripts/make_series_entry.py --fid $ZULF_DATA_DIR/fde3fbb2-average_fid.npy --id ethylenediamine       --out runs/series/eda --ranges "85,119.6;120.4,179.6;180.5,237.5"
Then look at runs/series/eda/phased.png. It must match the ethylenediamine panel of
runs/processed/amine_overview/all_amines_2026-10-05_morning.png (that figure is not in git; the shape is
described in docs/analysis/2026-10-05_amines_overnight.md).

I am waiting for Xuehan's confirmation before claiming zulf_core/physics items in docs/PLAN.md.
