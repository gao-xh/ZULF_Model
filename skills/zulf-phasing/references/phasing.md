# Phasing notes

Conventions
- `render.phasing.phase_correct(values, f, phase0, delay, acquisition)`
  multiplies by exp(-i (phase0 + 2 pi f (delay + reference))), where the
  reference is the crop start. Solver parameters (shared gain phase and
  `phase_delay`) use the same convention, so a complex fit's phase can be
  applied directly.

Options, in order of preference
1. Complex fitting with phase and delay fitted (no phasing needed).
2. Model-free phase: `scripts/phase_coherence.py`. For each trial delay, the
   phases of the strongest peak points (local maxima of abs(S) above a noise
   multiple) are corrected; the doubled angle removes the sign ambiguity of
   ZULF lines; coherence R = abs(sum w exp(2 i theta)) / sum w, with w the
   squared peak magnitude. Take the maximum; pick the overall sign so the
   strongest peak is positive.
3. `estimate_phase` (library): may stop at its +-10 ms delay bound; check.

Ambiguities
- With all peaks inside a narrow band (e.g. 129-149 Hz), phase0 and delay
  trade off; the coherence plateau can span 10 ms. Different settings then
  give nearly the same real spectrum inside the band and differ outside it
  (for example the shape of a weak 259 Hz line).
- A 2 pi wrap: a delay change of 1/f at the band centre is invisible there.

Ripple after phasing
- Anything smooth at the record start becomes a ripple of period
  1/(crop + delay) after first-order correction (0.05 s crop -> about 20 Hz;
  0.10 s -> about 10 Hz; a no-signal dataset shows it too).
- The solver's background columns follow the correction for phased spectra,
  so the background is smooth in the record frame (commit 0aeb959). Before
  that fix, real-only fits used narrow lines to chase the ripple.
- For display, the ripple is expected; do not interpret it as lines.
