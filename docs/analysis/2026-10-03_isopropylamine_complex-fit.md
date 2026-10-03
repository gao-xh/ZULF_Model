# Isopropylamine: complex known-structure fit with the series fitter (2026-10-03)

- Data: 10000-scan average `/home/user/zulf_data/isopropylamine/average.npy`
  with `0.ini` (NMRduino, 4000 Hz, 65516 points, 16.4 s; sequence
  standard_zf_4000Hz_no_dead.seq; not committed). Isopropylamine is the
  development molecule, never a blind test.
- Question (user): fit isopropylamine with the current fitter (complex data,
  model rendered through the same processing, missing-peak penalty, analytic
  Jacobian, multi-start), and see whether the NH2 group leaves a trace (fast
  against slow N-H exchange).
- Code: scripts/fit_joint_series.py (new: one spectrum with `--shape free`,
  per-entry fit `ranges`, `--exchange slow|fast`; same options in
  reliability_series.py); scratchpad ipa/explore.py, scan.py, apod_view.py,
  make_series.py.

## Processing (checklist of skills/zulf-fid-processing)

1. Sampling rate 4000 Hz (ini); 65516 points.
2. Diagnostics: switching edge 3.497 ms (half height), ADC plateau to
   2.75 ms, ringing to 54.75 ms, first point anomalous.
3. Crop-start scan (record 8 s, no window, zero fill 4, magnitude of the
   strongest line at 133.3 Hz): FWHM 1.09 Hz at 30 ms (ringing, peak/noise 4),
   0.53 Hz at 50 / 75 / 100 ms, 0.56 at 150 ms, 0.62-0.66 at 200-500 ms.
   Chosen: 0.1 s (plateau, 45 ms after the ringing).
4. Window (crop 0.1 s): noise/max 0.0084 / 0.0050 / 0.0043 / 0.0040 at
   0 / 0.3 / 0.6 / 1 1/s; the magnitude FWHM of the 133.3 Hz line grows from
   0.53 to 1.0 Hz at 0.3 1/s because its shoulder merges (apod_view.png).
   Chosen: 0.3 1/s (as for the pyridine series), record 0.1-8.1 s, zero
   fill 3 (0.0417 Hz grid).
5. Phase: instrument calibration (configs/confirmed_samples.json, same
   sequence): delay = edge - 0.033 ms = 3.464 ms, phase0 176.3 deg. The low
   band is absorptive; the 245-257 Hz band alternates in sign (phased.png).
   The complex fit refines a residual zero-order phase in the gain.
6. Bands: 115-160 Hz and 235-262 Hz above noise; mains harmonics at 60.06
   and 240.06 Hz. Fit ranges 85-119.6, 120.4-165, 232-239.6, 240.5-265 Hz
   (120 and 240 Hz left out; 85-115 Hz kept for a possible 15NH2 line near
   1.5 |1J(N,H)|).

## Fits

    python scripts/fit_joint_series.py --series series_ipa.json --real-only false --shape free \
      --exchange fast|slow --range 85,265 \
      --structure '{"motif": "(CH3)2CH-NH2", "one_bond": {"C1": 133, "C2": 125, "N1": -65}}' \
      --signal-threshold 2.5 --signal-taper 4.0 --peak-penalty 5 --peak-smooth 0.03 --peak-min-sigma 3 \
      --starts 64 --spread 2.0 --seed 71 --workers 2 --out runs/processed/ipa_fast|ipa_slow

Fast exchange: the NH2 protons dropped (13C at C1 and at the two equivalent
methyl carbons; 6 free couplings). Slow exchange: NH2 kept (also the 15N
isotopologue; 12 free couplings). No priors. Peak tops for the missing-peak
rows: 16. Started 2026-10-03 04:37 UTC.

## Results

| fit | starts | best total | relative residual | next best |
|---|---|---|---|---|
| fast exchange | 64 (618 s) | 0.1456 (2 starts) | 0.376 | 0.3178 |
| slow exchange, generic starts | 64 (1791 s) | 0.1828 | 0.426 | 0.2158 |
| slow exchange, nested start (the fast solution, N-H couplings 0) | 1 (66 s) | **0.1271** | 0.353 | - |

The slow model contains the fast one (N-H couplings 0 decouple the NH2 protons
from the 13C isotopologues), so the generic slow starts missed that basin: their
N-H couplings started at the motif values (3J(H1,N-H) 5.5, 2J/3J(C,N-H)
-4.5 / 4.5 Hz).

Couplings (Hz):

| coupling | fast | slow (nested) |
|---|---|---|
| 1J(C1,H1) (CH) | 133.45 | 133.43 |
| 1J(C2,H2) (CH3) | 124.44 | 124.42 |
| 3J(H1,H2) | 6.11 | 6.18 |
| 2J(C1,H2) | -4.21 | -4.28 |
| 2J(C2,H1) | -1.69 | -1.76 |
| 3J(C2,H3) (to the other methyl) | 5.18 | 5.20 |
| 2J(C1,HN) | - | 0.83 |
| 3J(H1,HN) | - | -0.71 |
| 3J(C2,HN) | - | 0.56 |
| 1J(15N,HN) | - | -64.89 |
| 2J(15N,H1) / 3J(15N,H2) | - | -6.39 / -2.18 |

Decay rates: fast 1.86 (13C at C1) and 3.08 1/s (methyl 13C); slow 1.53, 1.95
and 7.3 1/s (15N). Model delay -3.62 / -3.57 ms (edge -3.46 ms).

All couplings of the carbon skeleton agree within 0.1 Hz between the two
regimes. The slow fit is 13 % lower with 6 more couplings, but uses the freedom
for couplings of the NH2 protons below 1 Hz (slow exchange would give
3J(H1,N-H) of about 5-7 Hz) and for a broad 15N line (7.3 1/s) near 98-100 Hz
that the data do not show clearly (85-112 Hz panel). Remaining misfits are the
same in both: the lines at 121.6 / 122.3 and 133.3 Hz are sharper and taller
in the data, 131.8 Hz is too high in the model, the 251 Hz line is narrower in
the data; one decay rate per isotopologue cannot give lines of different
width.

Figures: runs/processed/ipa_fast/spectrum_fit.png,
runs/processed/ipa_slow_nested/spectrum_fit.png (real part with a display
AsLS baseline subtracted from data and model, imaginary part, 85-112 Hz);
run logs runs/processed/ipa_{fast,slow,slow_nested}/RUN_LOG.md.

## Conclusion

The carbon-skeleton couplings of isopropylamine are the same in both exchange
regimes: 1J(CH) 133.4, 1J(CH3) 124.4, 3J(H,H) 6.1-6.2 Hz. No resolved NH2
structure: no clear 15N line and N-H couplings under 1 Hz, consistent with
fast (or intermediate) N-H exchange in the neat amine. The 13 % gain of the
slow model is not evidence of slow exchange; it can absorb line-shape errors.

## Open points

- Line widths: lines of different width within one isotopologue (an intermediate
  N-H exchange rate would broaden lines unevenly; or a Gaussian/Voigt part).
- Exchange rate fitted for the NH2 group (Liouville model, D45).
- Uncertainty budget as for the pyridine series.

Commit: see git log (this file).
