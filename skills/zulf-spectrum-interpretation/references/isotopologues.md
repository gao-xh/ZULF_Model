# Natural-abundance isotopologue bookkeeping

- Natural abundances: 13C 1.07 %, 15N 0.37 %. At natural abundance each
  molecule carries at most one rare spin; double labels are negligible.
- A molecule with carbons C1..Cn gives n single-13C isotopologues. Equivalent
  carbons add: two equivalent methyls make that isotopologue twice as
  abundant.
- All isotopologues of one molecule share every proton-proton coupling: tie
  them in `RefineSettings.ties` (group names differ per component because the
  heteronucleus position changes the group order; map them explicitly).
- The solver renders each component per molecule, so fitted amplitudes
  (abs of `gains`) compare abundances directly. Expected ratio for two
  distinct protonated carbons: 1 : 1; for 15N vs one 13C site: about 0.35.
- Carbons without protons (C=O, COOH) give only small-coupling lines, usually
  below 20 Hz; they are hidden by SG baseline removal unless the low-frequency
  region is processed with an exponential baseline model instead.
- A missing band is information: no band means no protonated carbon of that
  type in the molecule (within SNR).

Example (sample e3d282da): CH carbon 1J 145.4 Hz and CH3 carbon 1J 129.8 Hz,
2J(C,H) -4.1 / -4.5 Hz, 3J(HH) 6.6 Hz shared; amplitudes 1 : 0.77 (complex
fit); isopropyl alternative fitted 0.55 where 2.0 is required -> rejected.
