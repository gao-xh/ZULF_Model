# Blake pyridine series: field-switch sequence check (2026-09-29)

- Data: none (simulation of the acquisition protocol); the spectra concerned are
  the Blake pyridine series and pyridinium (processed real parts, not committed).
- Question: does the Ajoy-lab sequence produce the state our `sudden_drop` model
  assumes (sum gamma I along the detection axis, sudden removal)?
- Code: 4bd904e (`scripts/field_switch_check.py`); 777f859 (pulse test).

## Setup

Sources: Andrews et al., PNAS Nexus 4, pgaf187 (2025): 9.4 T prepolarisation,
shuttle in GF1 = 80 uT (x), GF1 ramped off over 30 ms onto GF2 = 1 uT (-y, OPM
axis), GF2 switched off suddenly, no pulse. arXiv:2604.26071 (search snippet
only): a 90 deg 1H pulse is mentioned.
Checks: (1) a DC pulse (90 deg 1H, 22.6 deg 13C) or a pi pulse on 13C vs the
sudden drop, pyridine-2-13C; (2) brute-force 64-dim propagation of the full
sequence for pyridine-2-13C (57 lines) and -4-13C (18 lines), hold at 1 uT from
0 to 300 ms and the dephased long-hold limit; line amplitudes compared with the
ideal sudden drop after a free complex scale and time shift.

    python scripts/field_switch_check.py C2
    python scripts/field_switch_check.py C4

## Results

- Pulses: identical frequencies and normalised amplitudes (overlap 1.0), scale
  only (0.31 with sign -1; 1.67). Test
  `test_dc_pulse_only_rescales_two_nucleus_types` (fails by 9 % with 15N).
- Sequence, relative residual of the amplitudes vs sudden drop:

| hold at 1 uT | C2 | C4 |
|---|---|---|
| 0 ms | 0.019 | 0.018 |
| 1-30 ms | 0.03-0.15 | 0.07-0.12 |
| 100 ms (coherent) | 0.33 | 0.37 |
| 300 ms (coherent) | 0.53 | 0.59 |
| long, dephased | 0.092 | 0.076 |

- Only ~26 % of the magnetisation ends along y (last part of the rotation not
  adiabatic): an overall scale only.

## Conclusion

A pulse is harmless for 1H/13C isotopologues. Time spent at 1 uT before the
switch-off changes intensities inside the multiplets by several to tens of
percent, enough for the small couplings to shift by Hz. The hold time is not
given in the paper.

## Open points

- Ask for the sequence timings (GF1 ramp vs arrival, GF1-off to GF2-off delay,
  GF2 fall time) and whether the June 2026 data used a pulse.
- Or fit a preparation hold time t_hold.

Commit: 4bd904e
