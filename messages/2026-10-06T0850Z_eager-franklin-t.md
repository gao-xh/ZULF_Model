From: eager-franklin-t
To: xuehan
Re: 2026-10-06T0821Z_eager-franklin-t.md
Subject: all cloud data sent to Xuehan as zulf_cloud_data.zip (session chat), not through git

On Xuehan's request I sent every data file this cloud session used, plus all fit results, to Xuehan in
the session chat. The file is zulf_cloud_data.zip, 21 MB, sha256 first 16: e67673d8ed08d1fb. Experimental data stays
out of git, including this branch (AGENTS.md).

Contents (layout of the local data README; see MANIFEST.md and SHA256SUMS inside):
- raw/<compound>-<id>/: the 8 averaged FIDs: fde3fbb2, 4322bdfc, e66a4b08, 7ad4aafb, b683220d, e3d282da,
  7dc9a043 (blind, no signal), isopropylamine ae68d8cd-average.npy.
- raw/pyridine-series-blake/: fids/<x>_pyridine_csv (Blake's raw FIDs as CSV, x = 02..100) and spectra/
  (Blake's processed spectra 6_3_26_pyridine_<x>, 6_14_26_pyridinium_6M).
- original/isopropylamine-ae68d8cd/0.dat, 0.ini: one scan file and its settings.
- fits/<run>/: fit.json, RUN_LOG.md and J_table.csv of all 129 runs; fits/_series/: the amine series
  files (ranges, phase).
Check after unpacking: cd zulf_cloud_data && shasum -a 256 -c SHA256SUMS
Then move raw/* and original/* into ~/research/zulf/data/ and the runs into runs/processed/.
The RUN_LOG commands carry cloud paths: replace --series with a rebuilt series file (make_series_entry.py)
and --from-joint with the fit.json from fits/.
