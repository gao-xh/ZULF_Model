# Decision log

Each entry records what was decided, why, and what would change it.
Open questions are listed at the end and mirrored in `docs/PLAN.md`.

## D1. New repository, no runtime dependency on earlier projects (2026-09-25)

Physics, rendering and refinement are rewritten here, borrowing ideas from
`ZULF_Analysis_Tools` and `ZULF_NMR_Suite`. Importing them would couple this
package to molecule-specific code (`jfit.py` parameter names) and to their
storage layout. Cross-checks against those implementations belong in optional
tests that are skipped when the other repository is absent.

## D2. One problem specification object (2026-09-25)

All problem dimensions (allowed spin counts, maximum spins per component,
allowed nuclei, maximum components, maximum equivalence-group size, J bin
layout, spectral grid) live in `zulf_model.spec.ProblemSpec`, loaded from
`configs/problem_v1.json`. Generator, codec, models, trainer and solver read
only this object. Extending to other spin counts or nuclei is a configuration
change followed by retraining; code must not contain literal spin counts.

## D3. Frequency units in the Hamiltonian (2026-09-25)

The Hamiltonian is built in Hz so eigenvalue differences are transition
frequencies directly. Propagators use `exp(-i 2 pi H t)`.

## D4. Default observation protocol `sudden_drop` (2026-09-25)

gamma-weighted z preparation and detection, no pulse, per-molecule trace
normalization. Equivalent to the x-axis convention of `ZULF_Analysis_Tools`.
Must be confirmed against the laboratory sequence (open question Q1).

## D5. Absolute per-molecule amplitudes for mixtures (2026-09-25)

Component amplitudes are not renormalized to unit sum. Relative contributions
are abundance times per-molecule signal, then multiplied by a random
experimental response in the generator and by fitted gains in the solver.
`ZULF_Analysis_Tools.jfit` normalizes weights per isotopomer, which discards
the physical relation between isotopologues; `nitrogen.py` keeps absolute
weights. This package follows the latter.

## D6. Exact analytic rendering through the acquisition operator (2026-09-25)

The renderer evaluates finite-record DTFTs of damped oscillations analytically,
with SG mirror-edge and mean-removal corrections, so data generation and
refinement do not depend on time-domain synthesis. Time-domain synthesis
remains the test reference.

## D7. Rule-derived J uses bond-path classes averaged over graph symmetry (2026-09-25)

Couplings are sampled per atom pair by path length and hybridization, then
averaged over pairs with identical Weisfeiler-Lehman color classes of the
isotope-labelled graph. This produces magnetically equivalent rotor groups
(CH3, freely rotating CH2) and respects molecular symmetry broken by 13C or 15N
labels. Diastereotopic distinctions are not modelled in `gen-v1`.

## D8. Exchangeable protons dropped by default (2026-09-25)

N-H and O-H protons are removed from spin systems by default (fast exchange
decoupling), matching the carbon-bound isopropylamine skeleton model. The
generator option `exchangeable_protons` can keep them.

## D9. Two target representations (2026-09-25)

The CNN baseline predicts padded set targets (component slots, spin slots,
J matrices with masks). The Transformer decodes a grammar-constrained token
sequence over magnetic-equivalence groups. Both decode to `Interpretation`, so
evaluation and solver are model-agnostic.

## D10. Device policy for CUDA and Apple MPS (2026-09-25)

Simulation and rendering run in NumPy float64 on the CPU (MPS lacks float64);
networks train in float32 through `zulf_model.device`, which selects
cuda, mps or cpu. AMP defaults to bf16/fp16 on CUDA and off on MPS.
Checkpoints store CPU tensors and load with `map_location`.

## D11. Processing is optional; the pure route is the default (2026-09-25)

`Acquisition()` defaults to no crop, no SG, no mean removal and zero time
origin. `ContinuousRenderer` gives infinite-record Lorentzian spectra with no
sampling at all. Experimental recipes enable steps explicitly; training can
use `processing.mode` = pure, fixed or randomized.

## D12. Global J sign is unobservable (2026-09-25)

Negating every J negates all energies; with real preparation and detection
operators frequencies and amplitudes are unchanged (tested). Matching allows
one global sign flip; canonical forms make the largest observable
heteronuclear coupling positive. Relative signs remain observable (tested).
Pulsed protocols with complex operators must revisit this.

## D13. Signal-free components are rejected by the generator (2026-09-25)

Components without an observable heteronuclear coupling (for example only
protons) give no zero-field signal and are never emitted.

## D14. Matched linewidth continuation in the solver (2026-09-25)

Narrow zero-field lines make the J objective highly nonconvex (a 1 Hz error
moves a 2J line by 2 Hz, several linewidths). The solver refines through a
schedule of extra exponential apodization (default 10, 3, 1, 0 per second)
applied identically to data (re-processed from the stored FID) and model
(through the same acquisition operator). Only the final, unbroadened level
produces results. Broadening only the model was tested and fails (tests keep
the matched version). Without a stored FID the schedule is skipped and
flagged.

## D15. Real-signal gain response (2026-09-25)

A real FID's spectrum is not complex-linear in the component gain because of
the negative-frequency mirror term. Gains are always represented by two real
columns (response to 1 and to i); the shared-phase model uses
cos(phi) col0 + sin(phi) col1.

## D16. One tool registry for agents, CLI and MCP (2026-09-25)

Claude and GPT agents use the same typed tools as the command line. Tools
return JSON and file paths, never large arrays, and restate the interpretation
rules. Long-running tools run as background jobs. The MCP server and API
exports are generated from the registry so they cannot drift apart.

## D17. Processing recipes are proposed per dataset (2026-09-25)

Crop, SG window and nuisance terms differ between experiments.
`diagnostics.diagnose_fid` reports the first-point anomaly, saturation plateau,
ringing end (local envelope decay rate) and baseline fits, and lists candidate
recipes. None is chosen automatically; recipes are compared on held-out
acquisitions, and parameters that move between recipes are reported as
processing-sensitive.

## D18. Global pattern search supplies refinement starts (2026-09-25)

On the isopropylamine data every local refinement from the earlier default J
values converged to a wrong basin: the strongest simulated methine line sat
near 137 Hz while the data has it at 133.4 Hz, and no continuation level
bridged that gap. `solver.search` adds a phase-insensitive objective: the
magnitude spectrum of the processed record after a cosine start ramp and end
taper (which removes the broadband tail of an abruptly cropped record) is
compared with, per component, the magnitude of the complex sum of exact
transitions through the same taper and decay. Keeping the complex sum inside a
component is essential: antiphase neighbours inside a multiplet cancel, and a
sum of line magnitudes left a 30 to 50 percent residual even at the true
parameters (1 percent with the complex sum). Gains are nonnegative, each band
has a linear baseline, bands are weighted equally. Differential evolution
searches couplings and rates within the policy bounds; distinct low-cost
members (couplings at least 0.5 Hz apart) are polished and handed to the
complex refinement as starts. The search objective is never used for ranking.

## D19. Optional residual field during evolution (2026-09-25)

A low-cost shield leaves a residual field, and its effect on line positions
and intensities is a candidate explanation for model misfit. It is exposed as
`Protocol.field_ut` (default zero), computed exactly within the existing
sector decomposition and checked against full-space propagation. It is not
fitted by default; any use must be reported with the result.

## D20. Phase-corrected input as an option (2026-09-25)

Experimental spectra are phased by hand (zero- and first-order) and read as
the real part. `ProblemSpec.grid.phasing = "corrected"` reproduces this: each
training spectrum is multiplied by the exact inverse of its global phase,
phase delay and crop reference, then by a residual error drawn from
`ProcessingConfig.residual_phase0_range_rad` and `residual_delay_range_s`, so
the model tolerates imperfect manual phasing. Inference applies the same
operator with the operator's phase0 and delay (`render.phasing`); the crop
reference is added from the acquisition, so the values do not depend on the
crop. The negative-frequency mirror term of a real FID keeps its conjugate
phase, so correction is exact only up to that term (about 0.1 percent of the
peak in the tests). The default stays "none" (real and imaginary channels with
random phase). Channel count of the encoder follows the spec.

Automatic phasing (`render.phasing.estimate_phase`, `phasing = {"auto": true}`
at inference) handles lines of either sign: after a running-median background
removal in the crop frame (where the broadband tail of a cropped record is
smooth), each strong peak's phase is taken from a symmetric complex window sum
at a parabolically refined centre, and phase0 plus delay are fitted modulo pi.
Synthetic recovery: phase within 0.07 rad, delay within 0.05 ms. Background
removal and phasing must happen in this order: after first-order correction the
crop tail becomes a sinusoid with period 1 / crop time.

## D21. Free component ratios; targets use rendered weights (2026-09-25)

Relaxation during transfer, polarization and detection change isotopologue
ratios, so `PerturbationConfig.component_ratio_mode = "free"` draws each
component weight from `free_weight_log10_range` independently of abundance.
In every mode the contribution targets (set head and tokens, and component
order) are the rendered weights, abundance times response gain, not the
nominal abundances, so the network is not asked to predict a quantity the
spectrum does not show.

## D22. Pre-rendered shards next to live rendering (2026-09-25)

Rendering an experimental-route sample costs about 0.16 s on one core, too
slow to feed a GPU from a few workers. `training.prerender` writes shards
(float16 features, padded targets, JSON interpretations, manifest with the
spec digest) with seeds (seed, shard index), resumable and parallel.
`data.prerendered` (or a curriculum stage override `prerendered`) trains from
them; mismatched spec digests are refused. Live rendering stays the reference
and the default.

## D23. Set head predicts equivalence groups, not spins (2026-09-25)

The first CPU verification run (4-6 spins, 1000 steps) showed that the
spin-level set head never formed magnetic-equivalence groups (every candidate
had singleton groups, so structure coverage stayed at 0.008), mis-sized
components, and regressed J to averages because spin slots of the same nucleus
swap order between similar samples. The head now predicts, per component slot,
a (nucleus, group size) class for each canonical group slot plus group-pair J,
exactly the representation of the token grammar. Decoding is a small dynamic
program: groups packed at the front, nucleus order non-decreasing, total spins
in `spin_counts`, at least two nuclei; the system is built with
`SpinSystem.from_group_couplings`, so equivalence is exact. Codec set targets
gain `group_class`, `group_mask` and `group_couplings` (spin-level keys are
kept for analysis); older shards must be re-rendered.

## D24. The solver fits the first-order phase by default (2026-09-25)

The simulated recovery benchmark failed on 4- and 5-group systems even from
the true couplings (relative residual 0.24, a weak coupling pushed from
-0.39 Hz to its bound) because observations carry a random delay (first-order
phase, up to 2 ms) while the solver fitted only a shared zero-order phase; the
couplings absorbed the phase error. With the delay fitted, the same system is
recovered to 0.007 Hz from 0.5 and 2 Hz starts and the delay to 0.923 ms
(true 0.918 ms). Real data need it too (isopropylamine: about 6 ms).
`ParameterPolicy.fit_phase_delay` now defaults to True (bounds +/-10 ms).

## Open questions

- Q1. Exact laboratory preparation, pulse and detection sequence.
- Q2. Acquisition parameters for v1 (sampling rate, record length, crop and SG
  ranges) from the isopropylamine dataset. `configs/acquisition_v1.json` holds
  placeholders marked `frozen: false`.
- Q3. Literature ranges and signs for 1J(15N,1H), 2J/3J(15N,1H) and 1J(13C,15N).
- Q4. Operational definition of structure match, observational equivalence and
  unidentifiable samples (current implementation: `evaluation.matching`).
- Q5. Solver convergence tolerance after permutation matching.
- Q6. Component upper limit and mixture ratio ranges for v1 training.
- Q7. GPU model and memory for training.
- Q8. Whether solver tools are exposed through an MCP server.
- Q9. Blind-test molecule list.
