---
name: zulf-spectrum-interpretation
description: Read a processed ZULF J-spectrum qualitatively before fitting - which spin groups produce which lines (XH, XH2, XH3 at J, 3J/2, 2J), relative strengths of 13C and 15N lines, what fine structure and band width imply, and how natural-abundance isotopologues and their abundance ratios constrain the molecule. Use to turn a line list into testable hypotheses.
---

# Interpreting a ZULF J-spectrum

Reference tables: `references/line_patterns.md` (patterns, amplitudes,
typical coupling ranges) and `references/isotopologues.md` (natural-abundance
bookkeeping). Check any pattern claim with `compute_transitions` on a small
system rather than from memory.

## From lines to groups

- One line near J: 13C-H (XH). One line near 3J/2: XH2. J and 2J: XH3, the 2J
  line stronger than J when the methyl is isolated.
- A J line without the matching 2J line argues against an isolated methyl;
  check the 2J position against instrument lines first.
- Width and fine structure: a spread of tens of Hz needs several coupled
  protons; regular splittings of about 1-2 Hz point to small long-range
  couplings; strongly coupled proton networks (J_HH comparable to long-range
  J_CH) turn clean J / 2J lines into multiplets and can weaken a methyl's 2J.
- Place the band on the 1J(CH) scale (sp3 C-C, C-O/N, aromatic, aldehyde);
  these ranges generate hypotheses, they do not decide them.

## Heteronucleus strength

- Amplitude scales with (gamma_H - gamma_X)^2: 15N-H about 2.2x 13C-H per
  molecule; with 0.37 % vs 1.1 % abundance a 15N site gives about 0.75x one
  13C site.
- Exchangeable N-H and O-H protons usually decouple (fast exchange in neutral
  water): no 1J(N,H) lines; 15N then shows only small couplings at low
  frequency. Slow exchange (acidic, non-aqueous) would show 15NH3+ near
  73-75 Hz and its 2J near 146-150 Hz.
- Protonation state: amino acids in water near neutral pH are zwitterions
  (-NH3+, -COO-); alanine pKa about 2.3 / 9.7, so the neutral -NH2 fraction at
  pH 6 is about 2e-4. -NH2 dominates only above pH about 10 (15NH2 would give
  one line near 3/2 abs(1J), 1J about -65 Hz). Choose the group in a model from
  the solvent and pH, and ask for them.

## Natural abundance

- Each carbon position gives its own single-13C isotopologue; equivalent
  positions add up. All isotopologues share the proton couplings.
- Fitted amplitudes are abundance ratios - a structural test (one CH and one
  CH3 carbon: 1:1; isopropyl methyl:CH 2:1; 15N:13C about 0.35).
- Missing bands matter: no further band means no further protonated carbon.
  Non-protonated carbons (C=O, COOH) give only low-frequency lines.
- Chirality is invisible: enantiomers have identical J-spectra.

## Output of this step

A short list of hypotheses, each with the observation it explains, the line
it predicts that is not yet checked, and the abundance ratio it requires.
