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

## Display processing (not the fit's)

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
