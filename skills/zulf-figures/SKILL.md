---
name: zulf-figures
description: Make publication-style figures of ZULF fits - whole spectrum plus detail bands, experiment and simulation overlaid, isotopologue bars and structure insets, high-resolution display processing and a display-only baseline correction that removes the rolling baseline without lowering the peaks. Use when a fit result is to be shown to people (slides, paper, colleagues).
---

# Figures of ZULF fits

Script: `scripts/paper_figure.py` (the fit's own fit_joint_series options + `--fit RUN/fit.json --fid FID`);
helpers `scripts/figure_tools.py` (palette, RDKit structure insets); baseline `zulf_processing.anchor_spline_baseline`.
Worked example: isopropylamine, docs/analysis/2026-10-03_isopropylamine_complex-fit.md ("Paper-style figure").

    python scripts/paper_figure.py <fit options> --fit runs/processed/NAME/fit.json --fid DATA/average_fid.npy \
      --segments "104,158;222,276" --gains 1,2.5 --insets insets.json --colors '{"13C@C1": "#2f8f5b"}' \
      --title NAME --figure runs/processed/NAME/paper_figure.png

## Layout (what people liked)

- (a) whole range (5-300 Hz), experiment above simulation, nothing cut; dashed boxes mark the detail bands.
- (b) the line bands with a broken axis, experiment (black, thin, on top) and simulation (purple, thicker,
  alpha 0.8) overlaid on one baseline; a weak band scaled (x2.5) and labelled; legend at the right.
- Isotopologue bars above the line clusters (one colour per 13C isotopologue, fixed order from
  `figure_tools.PALETTE`, checked with the dataviz validator) and RDKit structure insets in the empty parts of the
  panels: the 13C and its protons in the isotopologue colour, the coupled carbon-bound protons grey, N/O-bound
  protons plain (decoupled by fast exchange). Insets need RDKit (`pip install rdkit`, extra `figures`).
- Power-line harmonics (60.06 Hz x n) are kept and marked with a small triangle, never cut.
- 1J on the structure: an inset entry with `"bond_note": "J(C2,HC2)"` and `--uncertainty
  RUN/coupling_diagram.json` marks one 13C-H bond in the isotopologue colour and writes 1J = J +- sigma under
  the drawing with a leader line (RDKit's own bond notes are too small to read).

## Display processing (not the fit's)

- A figure of a fit draws the fit's own model: in StudioSession call `fit_model_curve()` before `simulate()`, and
  check that the result has "fit_re"; without it `simulate()` returns only the quick-look simulation (other
  rates and gains), which made good fits look broad (2026-10-10, scripts/batch_figures.py).
- Display baseline per experiment (owner, 2026-10-10): tune the AsLS (or spline) baseline separately for each
  spectrum (SNR, line widths and spacings differ) and write its parameters on the figure; never reuse one set of
  parameters for different experiments. Subtract the same data baseline from the simulation drawn over it. With
  close lines (about 1 Hz) the crop wings between them cannot be removed by a smooth baseline without cutting the
  lines; say what is left (acetonitrile 2 mL: smoothness 0.25 Hz, p 0.05; 250 uL: pre-smooth 0.15 Hz, 0.3 Hz,
  p 0.05; a coarse smoothness makes broad humps beside the lines).
- An automatic AsLS choice must check every line cluster, not only the highest peaks: a flexible baseline
  (0.25 Hz) took two thirds of the height of the weaker, broader 2J band of triethylamine / isopropylamine in
  the simulation, which looked like a misfit for a day (2026-10-11; raw model and data agreed within 5 %).
  When a band looks misfitted in a baseline-corrected overlay, compare raw heights before changing the fit.
- Humps beside the lines (crop wings, broad structure; in data and model alike) stay with a plain AsLS. The
  overlay also tries a masked Whittaker smoother (line cores +-0.8 / 1.2 Hz weight 0, smoothness 1-3 Hz) and
  scores the residual structure beside the lines; the kind (AsLS or masked) is chosen for experiment and
  simulation together (different kinds made the heights disagree), the parameters separately (2026-10-12).
  With dense multiplets little is left between the lines and AsLS stays.
- Experiment vs simulation: take both routes and keep the better (owner, 2026-10-13). Common = the experiment's
  baseline subtracted from both (the drawn difference equals the true residual; the simulation can pick up the
  data's baseline wiggles); separate = each curve its own baseline (clean simulation, but different baselines can
  create a false difference: serine v11 triplet 25 % low, rms 0.055 -> 0.084). overlay_panels keeps the route with
  the smaller drawn difference within 6 Hz of the lines and prints it in the footer.
- A display baseline ("model", D62) built from a wrong fit draws structure into the data trace (LF_1 wings); look at
  the plain real part too before reading features off a corrected figure.

- Thinner lines: the whole record (crop 0.1 s to the end, 16.2 s), window 0.1 1/s (0 is sharper but noisy in weak
  bands), zero fill 4; phase from the instrument calibration.
- The model is rendered through the same processing with the fitted couplings, rates and delay; gains solved on
  the fit ranges and kept for the whole range (`predict(fixed_gains=...)`), so the wide view is a prediction.
- Remove the residual zero-order phase of the gains (modulo 180 deg) from data and model alike.

## Display baseline (same steps for data and model)

The rolling baseline (period about 1 / (crop + delay), about 10 Hz for a 0.1 s crop) and the negative crop wings
between close lines are processing effects of the record start (skills/zulf-phasing), not signal.

1. Anchor spline: anchors more than 1 Hz from every model line above 20 Hz, away from the power-line harmonics and
   from narrow data peaks near the lines (+-5 Hz); knots every 1.5 Hz; protected runs longer than two knot spacings
   bridged by a straight line; the noise-free model with `k_sigma=inf` (robust rejection removes its anchors).
2. Only within 2.5 Hz of the model lines (1 Hz taper): AsLS 1.5 Hz, which lifts the valleys between lines.

Pitfalls met on the way (isopropylamine): AsLS alone leaves the rolling baseline (wiggles at 20-60 Hz);
a spline with 6 Hz knots leaves 2-5 Hz ripples; a spline through anchors inside a cluster or with wide protection
pulls the peaks down (negative valleys); a cubic piece across a wide cluster overshoots (bridge it); protecting
every narrow data peak keeps the low-frequency wiggles. Write "display baseline corrected" in the caption
(`OUT.caption.txt` is written with every setting); the fit is never affected.

## Couplings on the molecule

`scripts/coupling_diagram.py --fit RUN/fit.json --smiles SMILES --atoms '{"C1": 1, "HC1": "H@1", ...}'
--variants OTHER/fit.json ... --noise-scale sqrt(L)`: one curve per coupling between its two groups (heavy atoms at
their RDKit positions, each proton group one node in the widest gap between its carbon's bonds), width ~ sqrt|J|,
blue J > 0, orange J < 0, solid / dashed / dotted for reliable / trend / not determined, labels J +- sigma with
sigma = sqrt((sqrt(L) x linearised)^2 + (half range over the model variants)^2) (L: integrated residual
correlation length in points; 6.6 for isopropylamine). The linearised error alone is 10-30 x too small; the spread
over refits with other model settings (rate families, extra couplings, exchange) dominates. Writes OUT.png and
OUT.json (J, sigma, parts, class).

## Automatic display baseline (owner, 2026-10-17)

One route for every figure: `overlay_panels.display_baseline(run, view)` (overlay panels, showcase figures).
The fitted lines narrower than 3 Hz FWHM (rate + window) are masked over +-(2 FWHM + 0.3 Hz), overlapping masks
join (a dense band is masked whole); a Whittaker smoother runs through the rest, smoothness by blocked
cross-validation with a floor of 1.5 Hz; every masked run is bridged by a straight line. Then the common or the
separate route for the fit, the smaller drawn difference kept. Check new spectra with
`python scripts/baseline_check.py OUT.png RUN:lo,hi ...` before using a figure.
Why each piece: a fixed +-0.6 Hz mask let the smoother dip under the isopropylamine 2J band (lines 3-5 Hz wide);
one smoother across wide masks dipped 0.3 of the signal under the lactic-acid bands; broad model "lines"
(10-40 1/s) masked triethylamine 150-220 Hz whole; a smoothness below 1.5 Hz followed the noise and the narrow
structure the fit misses (glycerol, DMEDA 30-100 Hz) and would hide it. Removing the processing distortion by a
model rendered without drift filter was tried and rejected (dispersive tails). The manual BASELINE_OVERRIDES
are no longer used.

