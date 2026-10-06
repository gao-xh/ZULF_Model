# acetonitrile (fringe field, 250 uL): first processing and known-structure fit with a fitted field (2026-10-06)

- Data: NMRduino run `ACN_fringe_field_ms1_LF_250uL` (5539 scans 0-5538, 2026-09-28 16:53 to 2026-10-02
  15:54; sequence standard_zf_2000Hz_no_dead.seq; **2000 Hz**, 65516 decoded points = 32.8 s). Stored locally
  (never committed) as `~/research/zulf/data/original/2026-09-28-acetonitrile-fringe-field-250ul/`; zip sha256
  43bae6be26f0932d. Averages under `~/research/zulf/data/processed/2026-09-28-acetonitrile-fringe-field-250ul/`.
- Question: first processing and fit of a new sample. Xuehan: the measurement was made in a field (fringe
  field), so the fitter needs a field parameter (default stays zero field).
- Code: 68a505f at the start (averaging, zero-field fits); 4fbac0a for the field fits (D47,
  `fit_joint_series --fit-field`). Run logs: `runs/processed/acn_*/RUN_LOG.md`.

## Setup

1. Averaging (`scripts/average_scans.py`, new): every scan decoded with `zulf_core.io.decode_dat`. Per-scan
   deviation = RMS of (scan - mean) on 0.1-4.1 s after removing each one's mean and linear trend; robust z
   (median, 1.4826 MAD). Median deviation 93 ADC, late noise 113 ADC (drifts slowly by +-1 % over the run);
   632 scans (11 %) have z > 5 (deviation 100-1000 ADC, scattered over the run). Two averages: all 5539 scans
   (`all-scans/`) and the 4907 scans with z <= 5 (`z5/`), each with even / odd half averages.

       python scripts/average_scans.py ORIGINAL/ACN_fringe_field_ms1_LF_250uL PROCESSED/all-scans
       python scripts/average_scans.py ORIGINAL/ACN_fringe_field_ms1_LF_250uL PROCESSED/z5 --exclude-z 5

   FID diagnostics (both): plateau end 3.5 ms, ringing end 43.5 ms.
2. Overview (magnitude, SG 201, crop 0.1 s, 8 s, 0.3 1/s): 13CH3 band at 136.30 Hz (207 sigma) with side lines
   134.45 / 137.1-138.2 Hz; 2J band split 271.2-271.4 / 273.7-274.5 Hz with a weak centre 272.5-272.8 Hz.
   Instrument lines 922.97 Hz (very strong), 294.0 Hz. Unassigned: a dispersive feature at 83-85 Hz, a sharp
   line at 154.1 Hz, a line at about 627 Hz, a broad hump at 20-22 Hz (low-frequency background region).
3. Processed spectrum (z5 average; halves and all-scans the same way):

       python scripts/make_series_entry.py --fid PROCESSED/z5/average_fid.npy --id acetonitrile \
           --out runs/series/acn --sampling-rate 2000 --ranges "125,150;255,292"

   Crop 0.1 s, record 8 s, 0.3 1/s, zero fill 3. The phase calibration of the config is for the 4 kHz
   sequence; on this 2 kHz run the calibrated spectrum is far from absorptive (2J band dispersive). The complex
   fit fits the zero-order phase (shared-phase gain) and the delay, so this does not enter the results; the
   fitted delay is -4.1 ms (edge 3.785 ms).
4. Model: 13CH3-C (chain C1 with 3 H, C2 without H), 14N not modelled (decoupled / no 14N spin), fast exchange
   irrelevant. In 125-292 Hz only the 13C@C1 isotopologue has lines (XA3: J and 2J). Common options:

       python scripts/fit_joint_series.py --series runs/series/acn/series.json --real-only false --shape free \
           --exchange fast --range 125,292 \
           --structure '{"compound": "acetonitrile", "chain": {"groups": [["C1","C",3],["C2","C",0]],
                         "bonds": [["C1","C2"]]}, "one_bond": {"C1": 136.3}}' \
           --rate-bounds 0.2,15 --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 \
           --peak-min-sigma 3 --starts 8 --workers 4 --max-nfev 300 --out runs/processed/NAME [options below]

   | Run | Options |
   |---|---|
   | acn_ch3_fast | (zero field, one rate) |
   | acn_field_<perp>_<z> | `--fit-field --field-start perp,z` for six starts |
   | acn_field_fam2 | `--fit-field --family-edges 135.5,137.2,200 --from-joint runs/processed/acn_field_0.005_0.08/fit.json` |
   | acn_field_fam4 | as fam2 with edges 135.5,137.2,200,272.1,273.1 |
   | acn_{even,odd,allscans}_field_fam2 | fam2 settings on the half averages and on the all-scans average |

## Results

Zero-field model (acn_ch3_fast): objective 0.5816, data residual 0.735 (all 8 starts equal); 1J 136.308 Hz;
the model puts one line at 272.6 Hz where the data have a split pair, and its J line is too low. Halves: 1J
136.325 / 136.295 Hz, residual 0.74 each.

Field grid on the magnitude (13CH3 sticks, Lorentzian 0.35 Hz, J 136.3 Hz): best near a longitudinal field of
0.08 uT (relative magnitude residual 0.667 -> 0.490); its 2J outer lines 271.33 / 273.88 Hz match the data.
A purely transverse field splits the J line into an equal doublet (not observed).

Fitted field, one decay rate:

| Field start (perp, z) uT | Objective | Residual | 1J (Hz) | B_perp (uT) | Bz (uT) |
|---|---|---|---|---|---|
| 0.005, 0.08 | **0.1992** | 0.439 | 136.280 | 0.0301 | 0.0530 |
| 0.05, 0.05 | 0.1992 | 0.439 | 136.280 | 0.0301 | 0.0530 |
| 0.02, 0.02 | 0.1979 | 0.441 | 136.284 | 0.0304 | 0.0496 |
| 0.02, 0.08 | 0.2016 | 0.441 | 136.280 | 0.0295 | 0.0525 |
| 0.04, 0.12 | 0.3963 | 0.611 | 136.241 | 0.011 | 0.125 |
| 0.08, 0.02 | 0.3114 | 0.547 | 133.816 | 0.158 | 0.797 |

The field reproduces the side lines of the J band (134.65 / 137.9 Hz) and the 2J triplet (271.5 / 272.6 /
273.5 Hz). Remaining misfit: line widths (sharp J centre, broader side lines and 2J lines).

Fitted field with rate families (acn_field_fam2, edges 135.5, 137.2, 200 Hz): objective **0.0439**, data
residual **0.211**, all 8 starts equal, no boundary hits; the end component search finds nothing better.

| Average | Objective | Residual | 1J (Hz) | B_perp (uT) | Bz (uT) | abs(B) (uT) | Rates (1/s): <135.5 / centre / 137.2-200 / 2J |
|---|---|---|---|---|---|---|---|
| z5 (4907 scans) | 0.0439 | 0.211 | 136.282 | 0.0371 | 0.0453 | 0.0586 | 0.78 / 0.43 / 2.76 / 2.16 |
| even half | 0.0706 | 0.273 | 136.286 | 0.0372 | 0.0431 | 0.0570 | 1.07 / 0.45 / 3.76 / 2.18 |
| odd half | 0.0592 | 0.248 | 136.279 | 0.0378 | 0.0460 | 0.0595 | 0.66 / 0.42 / 2.41 / 2.16 |
| all 5539 scans | 0.0435 | 0.211 | 136.281 | 0.0395 | 0.0445 | 0.0595 | 0.78 / 0.44 / 2.57 / 2.27 |

Splitting the 2J band into three more families (acn_field_fam4) changes the objective by 1 % (0.0435): not
kept. Linearised 1J error 0.001 Hz (known to be 10-30 x too small, WORKFLOW section 9); the half-to-half
spread is 0.007 Hz and the field model changes 1J by 0.026 Hz against the zero-field fit. The delay is -4.1 ms
in every field fit.

Figures:
- runs/acn/overview.png (magnitude 0-1000 Hz, all-scans vs z5)
- runs/series/acn/phased.png (calibrated phase, fit ranges)
- runs/acn/zero_vs_field.png, runs/acn/field_families.png (data, models and residuals, zooms 130-142 and
  266-280 Hz)
- ~/research/zulf/data/processed/2026-09-28-acetonitrile-fringe-field-250ul/{all-scans,z5}/scans.png
  (per-scan deviation and late noise)

## Conclusion

Conditional numerical result for one 13CH3 model with a uniform static field and frequency rate families:

- A static field of about 0.06 uT (B_perp about 0.037-0.040, Bz about 0.043-0.046 uT) explains the split J and
  2J bands; without it the 13CH3 model fails (objective 0.58 against 0.044). Both halves give the same field
  within 3 nT and 1J within 0.007 Hz.
- 1J(13C,1H) = 136.28 Hz (spread over halves and averages 136.279-136.286 Hz).
- The side lines decay faster than the J centre (0.7-3.8 against 0.43 1/s). This is the pattern of an
  inhomogeneous field (field-sensitive lines broaden more; the J centre moves only at second order with a
  longitudinal field). The rate families stand in for a field distribution, which the model does not have.
- Excluding the 11 % disturbed scans changes nothing measurable (0.0439 against 0.0435).
- B_perp and Bz trade against each other between fits (0.030 / 0.053 with one rate, 0.037 / 0.045 with rate
  families); the magnitude is better defined than the direction.

## Open points

- A field distribution (gradient over the sample) instead of rate families, and whether the field direction is
  known from the setup (fringe field of which magnet, orientation of the sample).
- The 2 kHz sequence needs its own phase calibration (configs/confirmed_samples.json has the 4 kHz one).
- Unassigned lines: 83-85 Hz (dispersive), 154.1 Hz, about 627 Hz; check them in a no-signal run of the same
  setup. The 13C@C2 (nitrile) isotopologue (2J(C,H) about -10 Hz: lines near 10 and 20 Hz) lies in the
  low-frequency background region; the processed grid starts at 20 Hz.
- 14N: the nitrile N is not modelled (no 14N spin); the 14N relaxation model of eager-franklin-t
  (physics.exchange) is zero-field only at first.
- What "ms1" and "LF" in the run name mean (asked).

Commit: (this file's commit)
