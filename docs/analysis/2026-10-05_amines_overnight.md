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

## 09:44 UTC check-in

The container restarted around 09:4x and every job was lost; the monitor records survived.

**Component search reruns** (from the monitor records; the end passes were not reached):

| Sample | Best objective | Previous best |
|---|---|---|
| isopropylamine (ipa_fast_fam_cs) | 0.0475 | 0.0475 |
| N-ethylmethylamine (nema_fast_fam_cs) | 0.0501 | 0.0501 |

The start passes found no better basin for any component. These fits are at their best for the current models.

### Step 3: ethylenediamine plus a second ethylene unit (eda_2sp)

Setup: two AA'BB' units with a free ratio (`--structure` list; 20 seeds over the second unit's 1J, J and J').
Best start before the restart: objective **0.0099**, against 0.0449 for one species (-78 %).

| Coupling (Hz) | Species 1 | Species 2 |
|---|---|---|
| 1J(C,H) | 131.32 | 126.75 |
| 2J(C,H) | -1.94 | -2.03 |
| geminal | -14.43 | -9.11 |
| J (same side) | 1.74 | 5.10 |
| J' (across) | 15.05 | 12.02 |

Gain ratio species 2 / species 1: 0.17.

What species 2 explains:
- the 190.6 / 192.1 Hz pair;
- the 197.8 Hz shoulder;
- the 208-238 Hz undulation, which is now fitted (its residual drops from about 0.005 to about 0.001).

What remains: the 199.5 / 200.2 Hz dip and peak structure.

Species 2 is a conditional result: a second, less abundant ethylene unit with smaller 1J, larger geminal and less
anti-dominated vicinal couplings. Candidates are a minor species in the sample (ethylenediamine carbamate from CO2,
or a protonated or hydrogen-bonded form) or a slowly exchanging second conformer population. Note that a carbamate
or a mono-protonated ethylenediamine would have two inequivalent CH2 groups, so a symmetric unit is only a first
model.

Figures:
- runs/processed/amine_overview/ethylenediamine_2sp.png
- runs/processed/amine_overview/ethylenediamine_2sp_components.png

Restarted: eda_2sp_refit (from this start, 6 starts, component search) and tea_fam121.

## 10:45-11:25 UTC

The container restarted again at 10:44 (uptime 0). The queue was restarted from the monitor snapshots.

### Step 4: ethylenediamine, two AA'BB' units, refit (eda_2sp_refit2)

Setup: from the 0.0095 snapshot of eda_2sp_refit; 6 starts, spread 0.3, max_nfev 300, component search start and
end. RUN_LOG: runs/processed/eda_2sp_refit2/RUN_LOG.md.

| Start | Objective |
|---|---|
| 0 (centre) | 0.00949 |
| **2** | **0.00676** (stopped at max_nfev 300) |
| 5 | 0.00798 |
| 3 | 0.00942 |
| 1 | 0.01004 |
| 4 | 0.04388 |

The component search (start and end) found no better basin. The data-region residual went from 0.189 (one unit) to
0.081.

| Coupling (Hz) | Species 1 | Species 2 |
|---|---|---|
| 1J(C,H) | 131.33 | 126.68 |
| 2J(C,H) | -1.98 | -3.22 |
| geminal | -11.74 | **-18.33** |
| J (same side) | 1.88 | 5.88 |
| J' (across) | 14.79 | 15.45 |

Gain ratio species 2 / species 1: 0.26.

Caveats (this is a conditional result):
- **Geminal of species 2.** -18.3 Hz is outside the usual sp3 CH2 range (about -10 to -15 Hz). It is probably
  compensating for something the model lacks, such as an asymmetric species.
- **Model-only line at 177.5 Hz.** The model puts a line there where the data show none (see the paper figure). The
  joint fitter weights model-only lines at only 0.2 (signal weighting from data peaks alone). This motivated
  step 5.

What remains: the 197.5-200.5 Hz structure (residual about +-0.01 on 0.1).

Figures:
- runs/processed/amine_overview/ethylenediamine_2sp_refit2.png (one unit vs two)
- runs/processed/amine_overview/ethylenediamine_2sp_refit2_components.png
- runs/processed/eda_2sp_refit2/paper_figure.png

### Step 5: model-line passes in the joint fitter (commit 0ecce84)

Change: `fit_joint_series --model-line-passes N`, which mirrors refine's `signal_model_passes`.
- After the fit, the best model's incoherent line envelope above 2 noise sigma joins the signal-weighting cores
  (`MixtureForward.add_signal_cores`).
- The best solution is then refit. This repeats until no point is added.
- Every solution is rescored under the final weights. The record keeps the objective under the original weights
  for comparison with older runs.
- Test: `SignalWeightingTests.test_added_signal_cores_match_a_forward_built_with_them`. The weights match a forward
  built with `signal_extra_hz` at the same points.

At the eda_2sp_refit2 optimum, pass 1 adds 1226 core points:
- the 164-180 and 180.5-189 Hz wings;
- small regions inside 190-204 Hz;
- 205-224.5 Hz.

The objective under the new weights at that same vector is 0.0110. Run: eda_2sp_mlp (one start, from the
eda_2sp_refit2 optimum, 3 passes).

### Step 6: ethylenediamine, AA'BB' plus an asymmetric X-CH2-CH2-Y (ABCD) unit (eda_abcd)

Motivation: a carbamate (H2N-CH2-CH2-NH-COO-) or a mono-protonated ethylenediamine has two inequivalent CH2 groups.

Setup:
- New motif "X-CH2-CH2-Y (ABCD)" (commit 518d2c8): one 1J per carbon, one geminal per CH2, only the a/b mirror.
- 15 seeds: 1J(C2,H) in 126.8-143 Hz times three (J, J') pairs. Species 1 and the rates come from the 0.0095
  snapshot.
- `--from-joint` now accepts models with couplings the previous fit lacks.

Running.

### Step 5 result: eda_2sp_mlp (model-line passes only)

- Pass 1 adds 1226 core points; pass 2 adds none.
- Objective (new weights) 0.01103 -> 0.01073. Under the original weights: 0.00694 (against 0.00676).
- The couplings barely move (species-2 geminal -18.07).
- **The 177.5 Hz line stays at full weight**, so the fit "needs" it in this model.

Line list: it is species 2's line at 177.54 Hz, relative amplitude 0.38 of that species' strongest line.

Diagnosis: this line shares the lowest rate family (below 195.3 Hz) with the sharp 190.6 / 192.1 Hz pair, so it has
to be as sharp as they are.

Figure: runs/processed/amine_overview/ethylenediamine_2sp_mlp.png.

### Step 7: rate-family edge at 185 Hz (eda_2sp_fam185)

Setup:
- Edges 185, 195.3, 199.5, 202.4.
- Started from the refit2 optimum, with the old rate of the lowest family copied into both new families.
- 4 starts, spread 0.15, component search at the end, 2 model-line passes.

Result:
- Best start 0.00663, then 0.00661 after the end pass.
- After the model-line pass: 0.00709 (new weights), **0.00662 under the original weights**. This is the new best
  for ethylenediamine; the data-region residual is 0.082.
- The 177.5 Hz line is gone: species 2 below 185 Hz now decays at 9.0/s, against 2.0/s for its 185-195 Hz lines.
- Couplings as in refit2:

| Coupling (Hz) | Species 1 | Species 2 |
|---|---|---|
| 1J(C,H) | 131.34 | 126.68 |
| 2J(C,H) | -1.99 | -3.22 |
| geminal | -11.84 | -18.00 |
| J | 1.91 | 5.85 |
| J' | 14.70 | 15.46 |

Boundary hits: species 1 185-195 Hz and species 2 195.3-199.5 Hz decay rates sit at the 15/s upper bound.

Figure: runs/processed/amine_overview/ethylenediamine_2sp_fam185.png.

**Band diagnosis of 196.5-200.8 Hz** (scripts/band_diagnosis.py; figure runs/processed/eda_2sp_fam185/band_197_200.png):
- The band holds 58.8 % of the cost; relative residual 0.077.
- Every lever is weak (largest band gain 1.7 %), so no free parameter acts on this misfit and the model lacks
  something there.
- Candidates:
  - separate rates for the 197.8 and 198.7 Hz lines of species 1, which share one family; queued as eda_2sp_fam198
    with an edge at 198.3 Hz;
  - a line-shape effect of the strongest lines;
  - an asymmetric main species (eda_abcd tests an asymmetric second species).

### Step 8: rate-family edge at 198.3 Hz (eda_2sp_fam198)

Setup:
- Edges 185, 195.3, 198.3, 199.5, 202.4, from eda_2sp_fam185; rates were remapped (each new family takes the rate
  of the old family containing it).
- 3 starts, spread 0.1, 2 model-line passes.

Result:
- Starts 0.00520, 0.00521. After the model-line pass: 0.00554 (new weights), **0.00509 under the original weights**.
- Data-region residual 0.072.
- **The couplings do not move** (all within 0.06 Hz of fam185); only the decay rates of the 197.8 and 198.7 Hz
  lines separate.

Boundary hits: species 2's two families between 195.3 and 199.5 Hz are at the 15/s bound. Species 2 broadens its
lines there out of the way, which suggests its model has more line intensity in that band than the data allow.

Figure: runs/processed/amine_overview/ethylenediamine_2sp_fam198.png. The legends of these figures give the
fit.json score, which is under the final (model-line) weights.

### Step 6 result: AA'BB' plus ABCD (eda_abcd)

- Best 0.00691, against 0.00676 for the symmetric second unit with the same family edges. The extra freedom gives
  no gain, so there is no evidence for an asymmetric species 2.
- Species 2 couplings:
  - 1J(C1,H) 126.93 and 1J(C2,H) 140.15;
  - geminals -17.9 and -15.9;
  - vicinal 12.0 and 7.2;
  - 2J(C2,H) +0.75.
- The symmetric model is kept.

### Ethylenediamine summary so far

| Run | Model | Objective (original weights) | Residual |
|---|---|---|---|
| eda_aabb | one AA'BB' unit | 0.0449 | 0.189 |
| eda_2sp_refit2 | two units | 0.0068 | 0.081 |
| eda_abcd | AA'BB' + ABCD | 0.0069 | 0.081 |
| eda_2sp_fam185 | two units, edge 185 | 0.0066 | 0.082 |
| eda_2sp_fam198 | two units, edges 185 + 198.3 | **0.0051** | 0.072 |

### Narrow rate families for the other amines

tea_fam121 (edges at 120.8, 121.2 and 121.6 Hz around the sharp 121.08 Hz CH3 line) is at 0.0383 mid-fit,
against 0.049 before. The same is queued for:
- isopropylamine: ipa_fam122, edges 121.4, 121.95 and 122.6 around the 121.67 / 122.29 Hz lines;
- N-ethylmethylamine: nema_fam121, edges 121.0 and 121.65 around 121.33 Hz.

Both start from their best fits with remapped rates (scratchpad amines/remap_families.py), 4 starts, component
search at the end and 2 model-line passes.

### Step 9: isopropylamine, the 118.5-126 Hz band

ipa_fam122 (narrow rate families at 121.4 / 121.95 / 122.6 Hz):
- Start 0 converged in 2 evaluations at the old objective 0.04753, so the run was stopped.
- Reason, from the line list: the model has **no line at 121.67 Hz**. Its 13C@C2 lines are at 122.13, 122.15 and
  122.36 Hz; the data peaks are at 121.67 (0.020) and 122.29 (0.036). Rates cannot fix a missing line.

Band diagnosis of 118.5-126 Hz (figure runs/processed/ipa_fast_fam2_4j/band_118_126.png):
- The band holds 60.6 % of the IPA cost; relative residual 0.367.
- Only conflicts: J(C2,HC1) wants -0.18 Hz and J(HC1,HC2) -0.12 Hz, but the rest of the spectrum resists.
- The whole-objective refit of the two levers gains nothing.

Second species (same motif, free ratio):
- The fit from 12 seeds stays at 0.04753. A copy of species 1 has zero gain and no gradient.
- Direct grid (scratchpad amines/scan_ipa_p2.py): species-2 1J(C2,H) 122-127 Hz, 3J(C2,HC3) 3.5 / 5.07 / 6.5 and
  J(HC1,HC2) 5 / 5.94 / 7, with the gains solved.
- Best 0.04685 (-1.4 %) at 1J 125.25, 3J 6.5.

A second isopropylamine-like species does not supply the 121.67 Hz line. Open: an N-H coupled (slow-exchange)
methyl band, or a different impurity.

Slow-exchange N-H test for isopropylamine (user's hypothesis that the N-H protons couple):
- Model: `--exchange slow` with the fast-fit couplings and rates.
- Grid: J(C2,HN1) -1..4.5, J(HC1,HN1) 0 / 3 / 6, J(C1,HN1) -3 / 0 / 3, with the gains solved.
- Best 0.0743, at all N-H couplings 0, against 0.0475 for fast exchange. The natural-abundance 15N-H lines are not
  in the data, and any 3J(C2,HN) splitting of the methyl lines makes the fit worse (0.5 Hz: 0.086; 1 Hz: 0.155).
- As for N-ethylmethylamine, the data favour fast N-H exchange.

### Step 10: N-ethylmethylamine and triethylamine with narrow rate families (running)

- nema_fam121 (edges 121.0 and 121.65 Hz added): best start 0.0467 mid-fit, against 0.0501.
- tea_fam121: four starts finished, best 0.03834 (previous best 0.049); the end pass is running.

Band diagnosis of the nema_fam121 snapshot (0.0467):

| Band | Share of the cost | Relative residual | Verdict |
|---|---|---|---|
| 13C@C3 (N-CH3) 127.6-134.8 Hz | 32 % | 0.215 | all levers weak: the model lacks something |
| 13C@C1 (C-CH3) 120.2-126.0 Hz | 28 % | 0.237 | conflicts with other bands (rate families) |

A pattern across the amines: after the fixes above, the largest remaining misfits are the methyl 13C bands near
1J (120-135 Hz). In isopropylamine this is the missing 121.67 Hz line; in N-ethylmethylamine, the N-CH3 band.
For the N-bonded carbons this points to the 14N item in the PLAN (spin-1 14N with quadrupolar relaxation).

### Step 10 result: triethylamine (tea_fam121)

- Edges 120.8, 121.2 and 121.6 Hz were added. The rates restarted, and the start component search brought the
  objective from 0.471 back to 0.070.
- Four starts, best **0.03834** (against 0.04897 for tea_fast_fam_cs, -22 %). The end component search found no
  better basin.
- Data-region residual 0.196 (was 0.215). The 121.08 Hz sharp CH3 line is now fitted with its own decay rate.

| Coupling (Hz) | tea_fast_fam_cs | tea_fam121 |
|---|---|---|
| 1J(C1,H1) (CH2) | 130.754 | 130.685 |
| 2J(C1,H2) | -4.659 | -4.864 |
| J(C1,HX) | -2.147 | -1.722 |
| 3J(H1,H2) | 7.200 | 7.168 |
| 1J(C2,H2) (CH3) | 125.003 | 125.004 |
| 2J(C2,H1) | -3.059 | -3.072 |
| J(C2,HX) | -0.499 | -0.459 |

The 1J values move by less than 0.07 Hz. The small couplings move by up to 0.4 Hz (J(C1,HX)), which is more than
their formal std (0.05). That formal std understates the model-dependence.

Left over:
- the 235 / 237.5 Hz small data peaks, which have no model line;
- the 124 Hz peak height;
- the 245-247 Hz shoulder.

Figures:
- runs/processed/amine_overview/triethylamine_fam121.png
- runs/processed/tea_fam121/paper_figure.png

### N-ethylmethylamine result (nema_fam121_final)

- The container restarted at 13:47. The best monitor vector (0.04579) was refit for one start, 80 evaluations,
  into runs/processed/nema_fam121_final.
- Objective **0.0458** (against 0.0501, -9 %); data-region residual 0.215 (was 0.224).
- Every coupling moves by less than 0.05 Hz.
- The 121.3 Hz data peak is still not fitted. As in isopropylamine, the model has no line there: a position
  problem, not a rate problem.

Figures:
- runs/processed/amine_overview/N-ethylmethylamine_fam121.png
- runs/processed/nema_fam121_final/paper_figure.png

## Morning summary (2026-10-05, about 14:00 UTC)

| Sample | Best run | Objective, start of night -> now | Data residual | What changed |
|---|---|---|---|---|
| ethylenediamine | eda_2sp_fam198 | 0.0449 -> **0.0051** | 0.189 -> 0.072 | second AA'BB' unit (ratio 0.26); rate edges 185 / 198.3 Hz |
| triethylamine | tea_fam121 | 0.0490 -> **0.0383** | 0.215 -> 0.196 | narrow rate families at the 121.08 Hz line |
| N-ethylmethylamine | nema_fam121_final | 0.0501 -> **0.0458** | 0.224 -> 0.215 | narrow rate families at 121.33 Hz |
| isopropylamine | ipa_fast_fam2_4j | 0.0475 -> 0.0475 | 0.217 | none (all tests negative) |

The ethylenediamine objective is under the original weights; the fit.json score after the model-line pass is
0.0055.

Overview figure: runs/processed/amine_overview/all_amines_2026-10-05_morning.png.

Paper figures:
- runs/processed/eda_2sp_fam198/paper_figure.png
- runs/processed/tea_fam121/paper_figure.png
- runs/processed/nema_fam121_final/paper_figure.png
- runs/processed/ipa_fast_fam2_4j/paper_figure.png

Negative results:
- No 15N isotopologue in ethylenediamine.
- No asymmetric (ABCD) second species in ethylenediamine.
- Isopropylamine has no second same-molecule species, no slow N-H exchange, and narrow families do not help.
- The NEMA N-H exchange rate test (earlier) found k -> fast.

Conditional points to keep in view:
- **Ethylenediamine species 2.** The second unit (1J 126.7 Hz, geminal -18.0 Hz, J 5.85, J' 15.5) is a numerical
  result, not an identified compound. The geminal is outside the usual range. Two of its decay-rate families sit
  at the 15/s bound, a sign that it has more line intensity at 195-199.5 Hz than the data allow.
- **Rate-family edges are nuisance parameters chosen by hand.** They lowered the objectives without moving the
  couplings by more than 0.07 Hz (1J) or 0.4 Hz (small couplings).
- **Common open misfit.** After these fixes, the largest misfit of three amines is a methyl 13C band near 1J with a
  data line the model does not have: isopropylamine 121.67 Hz, N-ethylmethylamine 121.33 Hz (plus the N-CH3 band),
  and the triethylamine 235 / 237.5 Hz small lines. Next candidates are 14N (spin 1, quadrupolar relaxation)
  for the N-bonded carbons and a common minor species (for example the CO2 carbamate) for all samples.

## Sign test of 2J(C,H) (user question: can 2J(C,H) be negative?)

Script: scratchpad amines/flip_2j.py.
- From the best fit, the sign of the two-bond C-H couplings is flipped and everything is refit (least squares,
  analytic Jacobian, the problem's bounds, 40 evaluations, one start). No prior on the couplings.

| Sample | 2J in the fit (Hz) | Fit | Flipped start | Flipped refit |
|---|---|---|---|---|
| triethylamine (tea_fam121) | J(C1,H2) -4.86, J(C2,H1) -3.07 | 0.0383 | 1.12 | 0.527 (2J +2.8 / +3.2, 3J(H,H) runs to 12.1) |
| ethylenediamine (eda_2sp_fam198) | species 1 J(C1,HC2a) -1.99 | 0.0051 | 1.03 | 0.198 (2J at the +3.7 bound, vicinal J distorted) |

Negative 2J (relative to positive 1J) is strongly favoured: the flipped fits are 14x and 40x worse, and the rest of
the network distorts to compensate.

Limits:
- the positive side is bounded at +3.7 Hz by the template bounds;
- one start and 40 evaluations, so the positive region is not searched exhaustively.

Literature on the sign (web search, 2026-10-05):
- Two-bond 13C-1H couplings are usually small and negative, about -4 to -6 Hz. They can be positive: aldehyde
  protons about +25 Hz, terminal ethynyl +40 to +50 Hz (general reviews; P. E. Hansen, Prog. NMR Spectrosc. 14,
  175 (1981), doi:10.1016/0079-6565(81)80001-5).
- Closest analogue to the amines: ethanol at zero field, 12CH3-13CH2-OH, with 2J(C,H) = -4.6 Hz, 1J 140.4 Hz and
  3J(H,H) 7.1 Hz in the simulations that matched the data. Source: "Optical detection of NMR J-spectra at zero
  magnetic field", arXiv:0901.4069. The CH3-13CH2-N couplings fitted here are -4.86 (triethylamine) and -4.50
  (N-ethylmethylamine).
- Orientation dependence: an electronegative substituent (O) on the 13C anti to the coupled proton contributes
  positively, gauche negatively (Schwarcz, Cyr, Perlin, Can. J. Chem. 1975, doi:10.1139/v75-262).
- No amine-specific 2J(C,H) values were found. DFT prediction of zero-field multiplets for natural-abundance 13C
  (Andrews et al., arXiv:2604.26071) is a route to independent values.

## What is ethylenediamine species 2? (user question, 2026-10-05)

In eda_2sp_fam198, species 2 has gain ratio 0.261 against species 1, and about 25 % of the integrated
|component| signal.

Its parameters do not fit the obvious candidates:
- **1J(C,H) 126.7 Hz is low for a CH2 bonded to N.** Species 1 is at 131.3 Hz. A protonated (ammonium) or
  carbamate N-CH2 should have a larger 1J, not a smaller one.
- **The geminal coupling of -18.0 Hz is outside the usual sp3 range**, and two of its rate families sit at the
  bound.

Candidates and their status:
- **Ethylenediamine carbamate / ammonium salt from air CO2.** Chemically plausible: ethylenediamine takes up CO2.
  But the carbamate is asymmetric. The ABCD model neither requires nor excludes asymmetry (0.0069 against 0.0068),
  and the low 1J argues against an N-CH2 next to a carbamate or ammonium N.
- **Protonated or hydrated forms in fast exchange** average with the free base, so they cannot give a separate
  spectrum.
- **A second conformer is excluded as a separate species.** gauche / anti rotation is far faster than the
  millisecond scale, so a slowly exchanging conformer is not plausible. Correction to step 3 of this log, which
  listed it as a candidate.
- **A stand-in for missing physics** (for example 14N coupling to the N-bonded 13C) is possible. The unusual
  geminal and the bound hits support this reading as much as a real compound.

Decisive checks:
- high-field 1H / 13C NMR of the same sample (a species at about 25 % would be obvious);
- a fresh sample kept under inert gas against an air-exposed one;
- deliberately adding CO2 (or a known carbamate) to see whether the 190.6 / 192.1 Hz lines grow.

Figure of the two species separately: runs/processed/amine_overview/ethylenediamine_two_species.png (scratchpad
amines/eda_species_split.py).
- In eda_2sp_fam198, species 2 has lines at 190.54, 192.00, 193.71 and 200.25 Hz. The 200.25 Hz sharp line is
  species 2's; the data-minus-species-1 trace matches it.
- Correction to step 3: in the final fit, species 2 has no line at the 197.8 Hz shoulder.
- What species 1 leaves at 197-200 Hz is an oscillating residual that species 2 does not explain.
