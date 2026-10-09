# ethanol (70 % rubbing alcohol, ULF run): first processing and known-structure fit (2026-10-08)

- Data: NMRduino run `70per_ethyl_alcohol_rubbing_alcohol_ULF` (8869 scans 0-8868, 2026-09-13 18:11 to
  2026-09-17 11:33; sequence standard_zf_8333Hz_no_dead.seq; **8333 Hz**, 65516 decoded points = 7.86 s).
  Stored locally (never committed) as `~/research/zulf/data/original/2026-09-13-ethanol-70pct-rubbing-alcohol-ulf/`;
  averages under `~/research/zulf/data/processed/` with the same name.
- Question: first processing and fit of a new sample (Xuehan: "run the ethyl in Downloads").
- Code: 15a7cbc at the start; sampling-rate fix in this commit. Run logs: `runs/processed/ethanol_*/RUN_LOG.md`.

## Setup

1. Averaging: `scripts/average_scans.py` on the run (all 8869 scans) and with `--exclude-z 5` (7630 kept). The
   per-scan deviation steps from about 110 to about 180 ADC at scan 6000 (late noise drops at the same scan):
   the setup changed there; 1232 of the 1239 excluded scans are after scan 5990.
2. Pitfall found: `analyze_sample.py` processed every FID at the config rate (4000 Hz). At 8333 Hz all
   frequencies came out x0.48 (the 923 Hz instrument line looked like a 443 Hz line, the 2J band like 120 Hz),
   and the default window (4000 samples) was 0.48 s. Fixed: the rate comes from `scans.json` / `.ini` next to
   the FID (`zulf_processing.find_sampling_rate`) and the sample-count defaults are rescaled to the same times
   (`scale_sample_settings`); `make_series_entry.py` reads the rate the same way and has `--grid`.
3. Processed spectrum (z5 and both halves; crop 0.1 s, record 7.5 s, 0.3 1/s, zero fill 3):

       python scripts/make_series_entry.py --fid PROCESSED/z5/average_fid.npy --id ethanol \
           --out runs/series/ethanol --record 7.5

   Bands: 120.6-136.1 (J band, 13CH3), 198.2-220.8 (3/2 J, 13CH2), 244.1-260.2 Hz (2J, 13CH3); a dispersive
   feature at 83 Hz (also in the acetonitrile run, instrument) is left out (`--range 120,261`). The 4 kHz phase
   calibration does not hold for the 8 kHz sequence; phase and delay are fitted.
4. Model: chain C1 (3 H) - C2 (2 H), natural abundance, OH in fast exchange (not modelled):

       python scripts/fit_joint_series.py --series runs/series/ethanol/series.json --real-only false --shape free \
           --exchange fast --range 120,261 \
           --structure '{"compound": "ethanol", "chain": {"groups": [["C1","C",3],["C2","C",2]],
                         "bonds": [["C1","C2"]]}, "one_bond": {"C1": 126.0, "C2": 141.0}}' \
           --rate-bounds 0.2,15 --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 \
           --peak-min-sigma 3 --starts 8 --workers 4 --max-nfev 300 --out runs/processed/NAME [options]

   | Run | Options |
   |---|---|
   | ethanol_zf | zero field |
   | ethanol_field | `--fit-field --field-start 0.02,0.05` |
   | ethanol_field_fam_{z5,even,odd} | as ethanol_field, `--starts 16 --family-edges 170,230` |

## Results

| Run | Objective | Residual | 1J(C1,H1) | 1J(C2,H2) | 2J(C1,H2) | 2J(C2,H1) | 3J(H,H) | B_t / Bz (uT) | delay (ms) |
|---|---|---|---|---|---|---|---|---|---|
| zero field | 0.276 | 0.508 | 125.93 | 141.26 | -5.39 | -5.18 | 9.03 | - | -2.2 |
| field | 0.170 | 0.416 | 125.63 | 140.88 | -5.70 | -5.22 | 9.45 | 0.019 / 0.029 | -1.8 |
| field, families, z5 | 0.169 | 0.416 | 125.62 | 140.87 | -5.70 | -5.21 | 9.45 | 0.019 / 0.029 | -1.8 |
| same, even half | 0.137 | 0.371 | 125.61 | 140.84 | -5.73 | -5.22 | 9.46 | 0.019 / 0.028 | -1.9 |
| same, odd half | 0.150 | 0.387 | 125.58 | 140.67 | -5.72 | -4.74 | 9.43 | 0.020 / 0.029 | -1.7 |

Couplings in Hz. The starts are multimodal (16 starts: objectives 0.169 to 0.461; 3 of 16 reach the best).
Rate families change nothing (0.1698 -> 0.1690). Halves agree within 0.05 Hz on 1J(C1,H1) and 3J(H,H), 0.2 Hz
on 1J(C2,H2), 0.5 Hz on 2J(C2,H1).

Figures: runs/series/ethanol/phased.png, runs/processed/ethanol_*/spectra.png,
~/research/zulf/data/processed/2026-09-13-ethanol-70pct-rubbing-alcohol-ulf/all-scans/scans.png.

### Field compared with the acetonitrile run (Xuehan: the field should be the same)

Same settings as ethanol_field, 16 starts:

| Run | Field (uT) | Objective | Residual | 3J(H,H) |
|---|---|---|---|---|
| ethanol_acnfield: field held at the acetonitrile value | B_t 0.0371, Bz 0.0453 (fixed) | 0.283 | 0.530 | 9.85 |
| ethanol_field_acnstart: fitted, started there | B_t 0.0192, Bz 0.0289 | 0.170 | 0.416 | 9.45 |
| ethanol_field_hh7: fitted, 3J(H,H) prior 7 +- 0.5 Hz | B_t 0.0192, Bz 0.0289 | 0.170 | 0.416 | 9.45 |
| ethanol_field_oh: OH kept (slow exchange) | B_t 0.000, Bz 0.051 | 0.297 | 0.510 | 7.49 |

The acetonitrile field (|B| 59 nT) fits clearly worse than the fitted 35 nT (objective 0.283 against 0.170);
started at the acetonitrile value the fit returns to the same 35 nT, and both halves give it within 1 nT. The two
runs differ in date (2026-09-13 against 09-28), sequence (8333 Hz against 2000 Hz) and run name (ULF against
LF). Caveat: the ethanol fit is not noise-limited (residual 0.42 against a noise-only estimate of about 0.15;
acetonitrile 0.21 against 0.20), so a model error may bias the field. The OH in slow exchange fits worse; fast
exchange stays.

### Field refit (Xuehan: fit the field again)

1. Start grid (ethanol_fgrid_<t>_<z>, 8 starts each, transverse 0.005 / 0.04 / 0.08 x z 0.02 / 0.06 / 0.12 uT):
   5 of 9 return the same best solution (objective 0.1698, B_t 0.0192, Bz 0.0289 uT, delay -1.81 ms,
   3J(H,H) 9.45 Hz). Other minima: B_t 0, Bz 0.067 uT (0.229; delay -3.76 ms, 3J(H,H) 7.14 Hz); B_t 0.08-0.10
   uT (0.297-0.323).
2. Delay held near the edge (`--phase-delay-bounds=-4.2,-3.0`, ethanol_fedge_*, four field starts incl. the
   acetonitrile value): every start gives B_t 0.0198, Bz 0.0316 uT (|B| 37 nT), objective 0.211, 3J(H,H)
   8.99 Hz; the delay sits on the bound -3.0 ms. Halves: B_t 0.0194 / 0.0203, Bz 0.0292 / 0.0311 uT.

The field is about 35-37 nT (B_t about 0.020, Bz about 0.030 uT) whether the delay is free or held near the
edge, in both halves and from every start near it; the acetonitrile value (59 nT) is not a minimum here. The
Bz-only 67 nT minimum with the physical 3J(H,H) is worse in both settings. The fit pushes the delay away from
the edge, so the delay / phase of the 8 kHz sequence and 3J(H,H) stay open.

### Acetonitrile rerun with the same field start grid (Xuehan: run acetonitrile again)

Current code (sampling rate 2000 Hz now read from scans.json; edge 3.785 ms as before), series
`runs/series/acn_rerun`, the one-rate settings of docs/analysis/2026-10-06_acetonitrile_field-fit.md, the same
9 field starts as ethanol (runs/processed/acn_rerun_*): zero field 0.5816 (unchanged); best 0.1979-0.1992 from
4 of 9 starts, B_t 0.030, Bz 0.050-0.053 uT (|B| 58 nT), 1J 136.28 Hz, delay -4.1 ms (edge + 0.3 ms); other
minima Bz 0.11-0.125 uT (0.396) and B_t 0.17 / Bz 0.79 uT (0.312). This reproduces the 2026-10-06 one-rate
result. With the same model (one rate, free delay) ethanol gives B_t 0.019, Bz 0.029 uT (|B| 35 nT): the two
runs differ by about 23 nT in |B|, mostly in Bz.

### Rate families for all three samples, and isopropanol (Xuehan: add rate families to every fit)

Isopropanol: run `70per_isopropyl_alcohol_rubbing_alcohol_ULF` (9301 scans 0-9300, 2026-09-17 12:06 to
2026-09-20 18:29, 8333 Hz, same sequence as ethanol); `original/` and `processed/` folder
`2026-09-17-isopropanol-70pct-rubbing-alcohol-ulf`, z5 keeps 8577 scans. Series `runs/series/isopropanol`
(record 7.5 s, ranges 63.8-81.5, 86-99, 115-119.7, 120.5-154.1, 238-258 Hz; the 83 Hz feature left out); model
`{"motif": "isopropyl", "one_bond": {"C1": 141.0, "C2": 126.0}}`, fast OH exchange, `--range 60,260`.

Rate families: acetonitrile with the edges of 2026-10-06 (135.5, 137.2, 200 Hz, from the one-rate field fit);
ethanol and isopropanol with `--family-edges auto` from their best one-rate field fit (ethanol 18 edges,
isopropanol 9; the edges used are in each RUN_LOG.md).

| Run | Objective | Residual | B transverse (nT) | B z (nT) | abs(B) (nT) | delay (ms) | Couplings (Hz) |
|---|---|---|---|---|---|---|---|
| acn_rerun_fam2 | 0.0439 | 0.211 | 37.1 | 45.3 | 59 | -4.12 | 1J 136.28 (reproduces 2026-10-06) |
| ethanol_fam_auto | 0.108 | 0.334 | 22.3 | 34.2 | 41 | -2.46 | 1J 125.70 / 140.93, 2J -5.63 / -4.42, 3J(H,H) 9.05 |
| ethanol_fam_auto_edge (delay bounds -4.2,-3.0 ms) | 0.126 | 0.357 | 19.1 | 47.2 | 51 | -3.0 (bound) | 1J 125.46 / 140.86, 2J -5.24 / -4.95, 3J(H,H) 9.11 |
| ipa_zf (one rate) | 0.505 | 0.663 | - | - | - | -3.73 | |
| ipa_field (one rate, best of 3 starts) | 0.431 | 0.632 | 36.3 | 37.7 | 52 | -1.76 | |
| ipa_fedge (one rate, delay bounds) | 0.440 | 0.639 | 36.0 | 38.4 | 53 | -3.0 (bound) | |
| ipa_fam_auto | 0.232 | 0.471 | 40.5 | 35.7 | 54 | -1.93 | 1J(CH) 142.02, 1J(CH3) 126.51, 2J -4.15 / -3.64, 3J(H,H) 5.91, J(C2,HC3) 3.35 |

With rate families the fields are 41-51 nT (ethanol), 54 nT (isopropanol) and 59 nT (acetonitrile); the
ethanol value depends on the delay (41 nT free, 51 nT near the edge), so the 8 kHz delay calibration is the
open point for the field. Ethanol 3J(H,H) stays at 9.1 Hz with families. Isopropanol: 1J(CH3) 126.5 Hz,
1J(CH) 142.0 Hz; misfit left at 120-122 Hz and 253-258 Hz.

Figures (scripts/paper_figure.py; the detail panels now share one scale, a804beb / 462d052 / 4197b97 fixed the
spectra.png and paper_figure display): runs/figures/{acn_full,ipa_full,ethanol_full,ethanol_fam_full}.png.

### Field direction scan (Xuehan: does the direction matter?)

Whole-grid series `runs/series/ethanol_full`, automatic rate families of ethanol_full_fam_0.04_0.04 (best whole-
spectrum fit: B transverse 6.6, B z 52.2 nT, objective 0.209), |B| held at 52.7 nT, angle theta from z held,
couplings, rates, phase and delay fitted (2 starts from that fit's couplings; runs/processed/ethanol_dir_<theta>):

| theta (deg) | Objective | Residual | delay (ms) | 1J(C1,H1) | 1J(C2,H2) | 2J(C1,H2) | 2J(C2,H1) | 3J(H,H) |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.209 | 0.447 | -3.72 | 125.453 | 141.853 | -2.943 | -6.048 | 7.688 |
| 15 | **0.2005** | 0.439 | -3.69 | 125.457 | 141.876 | -2.926 | -6.003 | 7.662 |
| 30 | 0.232 | 0.476 | -3.72 | 125.453 | 141.848 | -2.943 | -5.944 | 7.688 |
| 45 | 0.307 | 0.556 | -3.75 | 125.507 | 142.065 | -2.817 | -5.211 | 7.713 |
| 60 | 0.423 | 0.645 | -4.48 | 124.199 | 140.804 | +3.225 | +2.173 | 8.435 |
| 75 | 0.614 | 0.748 | -3.93 | 124.079 | 135.375 | +3.710 | +3.676 | 7.724 |
| 90 | 0.657 | 0.764 | -3.82 | 124.036 | 135.604 | +3.725 | +3.337 | 7.685 |

With the whole spectrum the direction is constrained: the field lies within about 30 deg of the detection axis
(best near 15 deg); 45 deg and more fit clearly worse, and from 60 deg the couplings jump to another solution
(2J positive, 1J(C2,H2) 135-141 Hz). Within 0-30 deg the couplings move by at most 0.03 Hz (1J), 0.1 Hz (2J) and
0.03 Hz (3J): there the direction does not change the reported couplings. The earlier trade-off between the
components (segment fits) came from the unfitted gaps, not from the physics.

## Conclusion

Conditional numerical result for one ethyl model with a uniform static field:

- A small field (about 0.035 uT, B_t 0.019, Bz 0.029) improves the fit clearly (objective 0.28 -> 0.17); it is
  the same in both halves.
- 1J(C1,H1) = 125.6 Hz and 1J(C2,H2) = 140.7-140.9 Hz are reproducible; 2J(C,H) about -5.2 to -5.7 Hz.
- 3J(H,H) = 9.4-9.5 Hz is reproducible but above the usual 7 Hz of freely rotating ethyl groups
  (docs/REFERENCES.md); with residual 0.42 and broad lines this is not trusted yet: it may absorb a model
  error (field distribution, linewidths, the phase / delay of the uncalibrated 8 kHz sequence).
- The fitted delay (-1.8 ms) is 1.6 ms from the switching edge (3.42 ms); on the 4 kHz sequence fit and edge
  agree within 0.05 ms. The 8 kHz sequence needs its own delay / phase calibration.

## Open points

- 3J(H,H): fit with an H-H prior (`--prior-sigma-hh`) and compare objectives; profile the objective over 3J.
- Delay: fix the delay at the edge and compare; calibrate the 8 kHz sequence on a known sample.
- Scans after 6000 (setup change): fit the all-scans and the scans >= 6000 averages separately.
- What changed at scan 6000; what "ULF" means for this run (field value known?).
- The 13C@C2 band (198-221 Hz) is broad; the 2J band (244-260 Hz) structure is not reproduced.

Commit: (this commit)
