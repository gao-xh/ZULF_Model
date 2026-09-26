# Line patterns at zero field (sudden drop, gamma-weighted)

Amplitudes below were computed with `compute_transitions` for isolated groups
(heteronucleus X coupled to n equivalent protons), 1J(CH) = 140 Hz,
1J(NH) = -90 Hz:

| Group | Lines (relative amplitude) |
|---|---|
| 13C-H | J (127) |
| 13CH3 | J (63), 2J (79) |
| 15N-H | abs(J) (275) |
| 15NH3 | abs(J) (137), 2 abs(J) (172) |

- XH2 gives one line near 3J/2.
- Amplitude scales roughly with (gamma_H - gamma_X)^2. gamma: 1H 42.58,
  13C 10.71, 15N -4.32 MHz/T, so 15N-H is about 2.2x 13C-H per molecule.
- Further couplings split and shift these lines by fractions of the small
  couplings; a strongly coupled proton network (J_HH comparable to the
  long-range J_CH) produces multiplets rather than clean J / 2J lines, and can
  weaken the 2J line of a methyl.
- Typical 1J(CH): sp3 C-C 125-135 Hz, sp3 C-O / C-N 140-150, aromatic
  155-180, aldehyde 170-230, alkyne about 250. Values here are general
  knowledge, not measurements; refine them.
- Typical 1J(15N,H): NH3+ in amino acids about -73 to -75 Hz; visible only
  with slow exchange.
