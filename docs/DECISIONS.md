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
