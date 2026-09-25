# Physical and numerical conventions

These conventions are frozen for generator version `gen-v1`. Changing any of
them requires a `docs/DECISIONS.md` entry, a new generator version and new
tests.

## Units

| Quantity | Unit |
| --- | --- |
| J couplings, frequencies, linewidth-free transition positions | Hz |
| Decay rates | 1/s (amplitude decays as `exp(-R t)`) |
| Time | s |
| Magnetic field (future extension) | microtesla |
| Gyromagnetic ratio | `gamma / (2 pi)` in Hz per microtesla |

## Hamiltonian

At zero field the Hamiltonian in frequency units (Hz) is

    H = sum_{i<j} J_ij  I_i . I_j

Energies from `eigh(H)` are therefore in Hz and transition frequencies are
energy differences in Hz. Time evolution uses `exp(-i 2 pi H t)`. The sign of
`J` enters the Hamiltonian directly. The sign of gamma never changes `J`.

## Nuclei

Registry defaults (Hz per microtesla): 1H 42.577, 13C 10.708, 15N -4.316,
19F 40.078, 31P 17.235. 15N has a negative gamma and spin quantum number +1/2.
Only 1H, 13C and 15N are in the v1 generator scope. 14N and 2H are not
modelled as spins in v1; their possible influence is a documented modelling
limitation, not an assumption that it is zero.

## Magnetic equivalence

A set of spins of the same nucleus is magnetically equivalent when every spin
in it has identical couplings to every spin outside it. Equivalent spins are
combined into collective spins. Couplings inside a group commute with every
observable used here and do not change transition frequencies; they are
unobservable and not fitted.

## Observation protocol (v1)

1. Prepolarization produces `rho0 proportional to sum_i p_i I_z,i` with
   `p_i = gamma_i` (high-field thermal state; the identity part is dropped).
2. The field is removed suddenly; no pulse is applied (protocol `sudden_drop`).
3. The detected signal is `s(t) = Tr[rho(t) D]` with `D = sum_i d_i I_z,i` and
   `d_i = gamma_i` (magnetometer sensitive to magnetization).
4. The absolute scale is per molecule: traces are divided by the full Hilbert
   dimension, so isotopologue contributions can be combined with abundances.

Because the zero-field Hamiltonian is isotropic, a z or x axis gives the same
spectrum. This matches `ZULF_Analysis_Tools` (x preparation and detection).
Ideal instantaneous rotations can be inserted after the drop through
`Protocol.pulses`; they are an extension point and not part of v1 data.

A residual static field during evolution is optional: `Protocol.field_ut =
(Bx, By, Bz)` in microtesla adds `-sum_n gamma_n B . I_n` (Hz) to the J
Hamiltonian, where z is the preparation and detection axis. The default is
exact zero field. The isotropy argument above no longer holds with a field:
a longitudinal field shifts lines only at second order, a transverse field
splits them at first order (for a 13C-1H pair by `(gamma_H + gamma_C) B / 2`).

For molecules with only two nucleus types and any longitudinal preparation
`sum_n w_n I_z,n` with z detection, the line shape of each component is
independent of the weights `w_n` up to one overall scale, because the total
`F_z` commutes with the zero-field Hamiltonian. Preparation weights or a
single ideal pulse therefore cannot change relative intensities inside one
1H-13C isotopologue; they only rescale it.

## Global sign

Negating every J leaves the zero-field spectrum unchanged for real preparation
and detection operators. Comparisons allow one global sign flip; canonical
forms make the largest observable heteronuclear coupling positive.

## Transition lists

Using eigenvectors with ascending energies, each pair `a < b` with
`f = E_b - E_a > 0` has complex amplitude

    A_ab = 2 rho_ab D_ba

and the noiseless signal is `s(t) = dc + sum_k Re(A_k exp(2 pi i f_k t))`.
For the default protocol all quantities are real. Transitions with identical
frequency (to `merge_tolerance_hz`, default 1e-7 Hz) are merged by summing
amplitudes; the merged frequency is the amplitude-magnitude weighted mean.
Transitions are never removed because of small weight unless a caller passes
an explicit, recorded relative threshold.

## Acquisition and preprocessing (single operator)

Every processing step is optional and off by default (the pure route). An
experimental recipe turns steps on explicitly in `configs/acquisition_*.json`.
A second pure route, `ContinuousRenderer`, gives infinite-record Lorentzian
spectra with no sampling at all.

The FID is a real sequence `x[m]`, `m = 0..N-1`, sampled at `fs`, recorded at
time `t_m = time_origin_s + m / fs` (default origin 0).

1. Optional Savitzky-Golay baseline subtraction on the full record with mirror
   edges: `y = x - savgol_filter(x, sg_window, sg_order, mode="mirror")`.
   Window 0 (default) disables this step.
2. Optional crop to samples `start_sample <= m < stop_sample`;
   `n = stop - start` (default: the full record).
3. Optional removal of the retained mean (default false).
4. Spectrum at any frequency `f`:
   `X(f) = (1/n) sum_{m=start}^{stop-1} y[m] exp(-2 pi i f (m - start) / fs)`.

Native bins are `f = k fs / n`. A zero-fill factor `z` gives `f = k fs/(n z)`;
this is an exact evaluation of the same finite record, not added information.
No apodization is applied. The first decoded point is retained.

Simulated spectra go through exactly this operator. Two exact backends exist:
the analytic backend evaluates it in closed form for damped oscillations
(mirror-edge SG corrections, mean removal); the time backend synthesizes the
FID with a type-1 NUFFT (relative error about 1e-12) and calls
`process_record` and `evaluate_spectrum`. The time backend is chosen
automatically when the grid is FFT-compatible or Gaussian broadening is used.
Mixtures and nuisances are summed as one FID and processed once.

## Model input features

Complex spectra on the model grid are divided by a robust scale (default the
root-mean-square of the magnitude over the grid) and presented as channels
`[real, imag]` with an optional `magnitude` channel. The scale is returned with
the features so predictions of relative contribution can be mapped back.

## Line broadening

Each component (or each transition family inside a component) has one
effective exponential rate R (homogeneous, Lorentzian, absorption FWHM R/pi)
and an optional Gaussian frequency spread with standard deviation sigma in Hz
(inhomogeneous: field inhomogeneity, scan-to-scan jitter). The FID envelope is
`exp(-R t - (2 pi sigma t)^2 / 2)`, giving Voigt lines. These are
phenomenological effective parameters, not a microscopic relaxation model.

## Spin ordering and equality of interpretations

Relabelling equivalent or identical-nucleus spins does not change the physics.
Comparisons of systems always search over permutations that preserve isotopes;
canonical orderings are conveniences for training targets only.
