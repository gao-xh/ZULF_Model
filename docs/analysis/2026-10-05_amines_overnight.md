# Amines, overnight optimization (2026-10-05 / 06)

Continuation of docs/analysis/2026-10-04_amines_complex-fit.md (see its "Overnight plan"). Data, processing and
fitter as there.

## Step 1: ethylenediamine, natural-abundance 15N lines (user: "try 15N")

Script: scratchpad amines/scan_15n.py.

- **Model.** The slow AA'BB' model is rebuilt at every 1J(15N,H) (start and bounds there), with the couplings and
  rate families of eda_aabb. The N-H couplings to carbon and to the CH protons are 0, so the 13C isotopologue is
  exactly the eda_aabb fit and only the 15N@N1 component changes. Its rate families are refit at every value.
- **Scan.** 55-135 Hz in 2 Hz steps. The strongest 15N line sits at 1.5 x 1J (|amplitude| about 0.01).
- **Baseline.** eda_aabb without 15N-H lines: 0.0449.

| 1J(15N,H) (Hz) | 15N line (Hz) | Whole objective |
|---|---|---|
| 65 | 97.5 | 0.0463 |
| 79 | 118.5 | 0.0454 |
| 119 | 178.5 | 0.0443 |
| 125 | 187.5 | 0.0652 |
| **127** | **190.5** | **0.0528** |
| 129 | 193.5 | 0.0639 |
| 133 | 199.5 | 0.0436 |

A 15N line placed on the pair (1J about 127 Hz) makes the fit worse: it is broad, and a single line where the data
show two sharp ones. No value improves the fit by more than 3 %, and the two small dips (119 and 133 Hz) are
outside any physical 1J(15N,H) of an amine. The 190.6 / 192.1 Hz pair is not the 15N isotopologue.
Figure: runs/processed/eda_aabb/scan_15n.png.

## Step 2: ethylenediamine, delay prior

eda_aabb_d holds the delay to the instrument prior (-4.6..-2.6 ms) and ends at the -2.6 ms bound.

| Run | Objective | Delay (ms) | Couplings (Hz) |
|---|---|---|---|
| eda_aabb (delay free) | 0.0449 | +0.58 | 1J 131.16, 2J -1.77, geminal -13.79, J 1.48, J' 16.04 |
| eda_aabb_d (prior) | 0.0452 | -2.6 (bound) | 1J 131.09, 2J -1.73, geminal -13.41, J 1.38, J' 16.48 |

The couplings do not depend on the delay. The 208-238 Hz undulation is the same in both, so the delay is not its
cause. The likely cause is the crop wings of the broad lines in the four coarse rate families, which have little
weight outside the line cores. Queued: eda_aabb_fam with 12 family edges (186-212 Hz), from eda_aabb.
Figure: runs/processed/amine_overview/ethylenediamine_delay.png.

## Running at 08:50 UTC

| Run | What it tests |
|---|---|
| ipa_fast_fam_cs | isopropylamine with the component search |
| nema_fast_fam_cs | N-ethylmethylamine, component search from the new basin |
| tea_fam121 | triethylamine with narrow rate families at the 121 Hz CH3 lines (120.8 / 121.2 / 121.6) |
| eda_2sp | ethylenediamine plus a second AA'BB' unit with a free ratio |

A job queue (scratchpad amines/queue.sh, queue.txt) starts the next script whenever fewer than four fits run.
