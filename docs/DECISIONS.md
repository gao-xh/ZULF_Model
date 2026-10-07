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

## D25. Refinement against experiments always tries opposite signs (2026-09-25)

Only the global sign of all J is unobservable. Relative signs are observable
whenever the coupling network connects the spins, even through a single
coupling of a few hertz (13C-15N-H-H test system: flipping only 1J(NH)
changes the spectrum by 19-26 percent; with no couplings between the CH and
NH parts it changes nothing). The evidence sits in small splittings that the
network misreads (CPU verification: 15N-1H strong median error 115 Hz, a sign
error) and that a local fit cannot cross, because a coupling whose interval
excludes zero never changes sign. `RefineSettings.sign_variants` therefore
adds, for every candidate, variants with one large coupling flipped and with
all couplings of one heteronuclear group flipped (`solver.variants`), removes
variants equivalent under the global sign, refines all of them and keeps every
branch. Rule for experimental work: enable it whenever J are refined against
an experimental spectrum, and report which sign variant won and by how much
(held-out ranking when available). Matching allows an independent sign per
connected block of the true coupling network (`evaluation.matching.
coupling_blocks`), since only disconnected blocks have unobservable relative
signs; the relabelling is also searched on |J|.

## D26. Simulated annealing as an alternative global search (2026-09-25)

`SearchSettings.method = "dual_annealing"` runs scipy's generalized simulated
annealing on the pattern objective (independent seeded runs until enough
distinct minima or the budget), next to the default differential evolution.
The solver's deterministic counterpart of annealing is the matched linewidth
continuation (D14): data and model are broadened identically and the
broadening is lowered step by step. Which global method recovers more
systems is to be measured with `scripts/solver_recovery.py`.

## D27. Generator output must not depend on the Python hash seed (2026-09-25)

`MoleculeGraph.wl_colors` used Python's `hash()`, which is randomized per
process for strings. Colour numbering, skeleton hashes and therefore family
ids and the train/val/test split differed between processes, and the same
seed produced different molecules. Pre-rendering with several worker
processes mixed split definitions, so the training-time validation set of the
CPU verification runs overlapped with training families: the logged coverage
(CNN+Transformer structure@10 0.45, J@10 0.30) is inflated. Re-evaluating
the same checkpoint in a fresh process gave structure@10 0.15-0.19 and J@10
about 0.01-0.02 (still not clean, since that process's split also overlaps
the workers'). Fix: a hashlib-based stable hash; a regression test runs the
generator under three PYTHONHASHSEED values and requires identical samples
and splits. All verification numbers recorded before this fix are
superseded; relative comparisons made inside one process (set versus
sequence model, beam versus solver rerank) keep their direction.

## D28. Guarded linewidth continuation (2026-09-26)

On a crowded 6-spin system (79 transitions, 0.16 Hz lines) refinement started
at the true couplings ended 5.9 Hz away with relative residual 0.49 when the
default continuation (extra broadening 10, 3, 1, 0 per second) was used, but
stayed at the truth (0.005 Hz, residual 0.138 = noise floor, 7 times faster)
without it: broadening merges neighbouring lines and moves the coarse
optimum out of the correct basin. Each start now also runs a direct
full-resolution fit (`RefineSettings.guard_continuation`, default True) and
the better full-resolution score of the two paths wins. Continuation remains
useful for starts several hertz off in sparse spectra (D14).

## D29. Core and model packages (2026-09-26)

The refinement and comparison side is used on experimental data without GPUs
and changes at a different pace from the learning side, but both must use one
forward model. The repository therefore holds two packages: `zulf_core`
(NumPy/SciPy: physics, rendering, solver, matching, I/O, diagnostics) and
`zulf_model` (PyTorch: spec, generator, codec, synthesis, models, training,
fine-tuning, proposers, agent tools). The dependency is one-way and tested.
Code is never copied between them. Moving `zulf_core` to its own repository
later only changes packaging, not imports.

## D30. Refinement on processed spectra (2026-09-26)

Spectra are often processed elsewhere (cropping, baseline, manual phasing) and
only the spectrum is kept. `ObservedSpectrum.from_spectrum(frequencies, values,
ranges, record, phasing, real_only)` accepts such a spectrum;
`io.load_spectrum_table` reads .npy/.npz/CSV tables and the agent tool
`refine_candidates` takes `spectrum` instead of `source`. The model is made to
match the data instead of undoing the processing: with `record` (sampling
rate, points, crop, SG) the finite-record lineshape is exact, without it ideal
Lorentzian lines are used; `phasing` (phase0_rad, delay_s in the
render.phasing convention) multiplies every model column by the same
correction phasor, and a real (absorption) spectrum is compared on the real
part only. Complex input and a record reproduce the FID route. Global pattern
search needs the FID and is skipped with the flag
`search_unavailable_without_fid`. For real-only systems the background
projection uses an SVD basis because the imaginary background columns vanish;
complex systems keep QR (a change there moved a weakly identified nuisance
rate, so it was reverted).

## D31. Analytic variable-projection Jacobian (2026-09-26)

Finite-difference Jacobians cost one full model evaluation (eigendecomposition
of every sector block plus rendering) per free parameter per iteration, and
they inherit any roughness of the reduced objective. `refine` now passes an
analytic Jacobian to the trust-region solver (`RefineSettings.jacobian`,
"kaufman" by default, "analytic" for the exact form, "finite_difference"):

* Couplings: `physics.derivatives.transition_derivatives` differentiates the
  signal Tr(U rho U^+ D) with the Daleckii-Krein form, so degenerate levels
  need no special handling (rotations inside a degenerate cluster do not
  change the signal). Each transition gets (dA, T): the derivative signal is
  Re((dA + 2 pi i t T) exp(2 pi i f t)).
* Rates and the phase delay: closed forms on the same transitions.
* Rendering: `Renderer.render_pair_directions` synthesizes all derivative
  signals with one batched NUFFT per rate group and pushes them through the
  same acquisition operator (SG, crop, mean, apodization, phase reference);
  `ContinuousRenderer` has the Lorentzian closed form.
* Gaussian widths and nuisance parameters: central differences of their
  columns only (no eigendecomposition). The continuous route with Gaussian
  widths falls back to forward differences and flags it.
* Elimination of the linear variables (amplitudes, shared phase, background
  and nuisance amplitudes) is exact: dl/dx = -H^+ (A^T J + dA^T r) with
  H = A^T A + S, where S carries the residual-weighted second derivatives in
  the shared phase. Without S and the dA^T r term this is Kaufman's
  approximation, exact only at zero residual.

The Jacobian matches central differences to 1e-6 off the optimum. Checking
this exposed that the shared phase with a clamped (zero) amplitude was found
by golden section only to about 1e-8 rad, which made the objective rough at
that level; it is now polished with the analytic phase derivative on the
active set.

Benchmark (`scripts/solver_recovery.py`, 8 generator systems of 4 to 6 spins,
8192 points, SNR 100 to 300, two starts per offset, one thread per run,
success = reaching the optimum found from the truth within 0.1 Hz):

| Jacobian | offset | reached | total time | median time | model evaluations |
|---|---|---|---|---|---|
| exact analytic | 0.5 Hz | 16/16 | 28 s | 1.4 s | 841 |
| Kaufman | 0.5 Hz | 16/16 | 26 s | 1.5 s | 794 |
| finite difference | 0.5 Hz | 16/16 | 59 s | 2.7 s | 7561 |
| exact analytic | 2 Hz | 11/16 | 162 s | 7.8 s | 4907 |
| Kaufman | 2 Hz | 11/16 | 126 s | 6.9 s | 3967 |
| finite difference | 2 Hz | 11/16 | 309 s | 15.2 s | 40075 |

All three reach the same points; the analytic forms are 2 to 2.5 times
faster. Kaufman's form is marginally faster than the exact one here and is
the default. The 2 Hz misses (3 systems) are basin failures shared by every
mode and need global search starts, not a better local step. Two systems
end 0.40 and 0.16 Hz from the truth from every start, including the truth
itself: at this SNR the data optimum is that far from the truth, so the
recovery benchmark now scores against the optimum refined from the truth
(`basin_of_attraction(reference="refined_truth")`, the script default)
and reports both distances.

## D32. Signal-focused weighting with model cores (2026-09-26)

`band_weighting="signal"` (d5543fc, 196e14f, 9e46dbd, 6595bad) scales the
residual by the noise sigma (MAD of the narrow excess) and weights points by
their distance to peak cores: 1 at a core, Gaussian fall-off (2 Hz) to 0.2.
Data cores are narrow features above 4 sigma. A data-only mask let a model
put lines where the data show none at low cost (e3d282da: 263-265 Hz), so
the model's own signal joins the cores and the fit is repeated until the
cores are stable. The model signal is read from the transition lists, not
peak-picked from the rendered spectrum: first per-line heights (missed
broad overlapping lines and a cancellation hole at 139 Hz), now the
incoherent envelope sum_k abs(g_c a_k) P_R(f - f_k) above 2 sigma, which
covers isolated lines of either sign, broad components and cancellation
holes. Cost: the envelope's tails cover most of each band, so the focus on
peaks is weaker there. Models are compared on `data_region_residual` (data
cores only), the same region for every model.

## D33. Fixed amplitude ratios (2026-09-26)

`RefineSettings.amplitude_ratios` (2a16dfe) merges components into one
column sum_c r_c col_c before the linear solve; the Jacobian merges the
derivative columns the same way (tested against finite differences). Used
for isotopologue abundances of a pure sample. With free amplitudes and free
per-component rates one component turned into broad background (e3d282da);
fixed ratios plus a tied rate prevent it. Free ratios stay the first test:
on b683220d they came out at the natural 1 : 1 by themselves.

## D34. Programmatic hypotheses (2026-09-26)

The blind analyses built isotopologue sets and their ties by hand; group
indices differ between isotopologues, which made the ties error-prone.
`zulf_core.hypothesis` describes fragments by labels (sites, proton groups,
couplings, symmetry) and builds the set automatically: one component per
symmetry-distinct labelled site and isotope, proton groups merged where the
labelled site's stabiliser makes them equivalent, couplings tied by symmetry
orbit and across isotopologues, natural-abundance ratios, and components
without couplings or without lines in the fitted ranges omitted with a
reason. Named couplings ('J(Ca,Ha)') replace solver names in reports.
Checks, extension moves, templates and reference couplings are registries so
that further cases add entries rather than code paths. Replay on b683220d:
builder-built H7 and H8 reproduce the hand-built predictions to 2e-13; the
checks flag the H7 rate asymmetry and the isopropyl abundance contradiction.
Stage 2 (band inventory, group candidates, fragment enumeration) and stage 3
(budgeted search over moves, knowledge matching) are planned in PLAN.md.

## D35. Hypothesis proposals from the spectrum; neural models as hint providers (2026-09-26)

Stage 2 of D34 (`propose_hypotheses`): band inventory (lines above 4 sigma of
narrow excess, instrument lines set aside, a noise floor of 1e-3 of the peak
for near noise-free input); group candidates per band from isolated X-Hn
patterns (XH at J, XH2 at 3J/2, XH3 at J and 2J), scored by partner lines
that must reach half their height predicted from the anchor; a line explains
another band within a multiplet width (10 Hz) only if its predicted height
reaches 0.4 of that band's peak; pattern priors (15N 0.6); minimal covers of
the strong bands (weak bands, SNR < 8, may stay open at a 0.1 score penalty);
all trees, isolated parts and doubled methyl leaves as topologies; couplings
beyond 1J from a registered prior. Both confirmed samples rank the correct
bonded CH3-CH first (lactic acid with no open band; alanine with the 74 Hz
band open), with isopropyl as the next alternative.

Neural models enter as `HintProvider`s (`ProposerHints` adapts any proposer;
`zulf_model.evaluation.model_hints` loads a checkpoint). Hints add a bounded
bonus (0.1) to agreeing candidates, add hint-only candidates at score 0, and
pass interpretations on as extra candidates; they never enter the final
ranking, which stage 3 takes from refinement and checks. The insight report
lists per band whether hints agree, hints that predict lines where the data
show none, and a per-source summary. On the two samples: the transformer
checkpoints give no one-bond couplings; the set model gives 7 group hints
per sample, 1 agreeing with the data (147.6 and 131.9 Hz) and 6 elsewhere.

## D36. Minor isotopologues and labelling schemes (2026-09-26)

Multiply labelled isotopologues (13C-13C, 13C-15N) are about 1 % of their
parents at natural abundance. They are weighted by their abundance, not by a
tuning factor: a `Labeling` (natural by default; `Labeling.enriched` for
uniform or site-specific enrichment) gives exact label-set probabilities
(labelled sites a, unlabelled sites 1 - a). Components below
`primary_fraction` (0.05) of the strongest set are minor: they follow the
primary component sharing most labelled sites through the solver's
`amplitude_map` (components x free amplitudes; generalises the fixed ratios
of D33 and also fixes ratios per part of a combined model) and a tied decay
rate, so they add no free parameters and cannot absorb background. Sets
whose relative amplitude is below `min_ratio` are omitted with the reason;
the proposal pipeline sets it to 2 / peak SNR (lines must be able to reach
2 sigma), which drops the 13C-13C set on the current data (0.0108 < 0.050)
and keeps it at higher SNR or in enriched samples. (Changed 2026-09-27 at
the user's request: the pipeline default is single-label natural abundance;
doubly labelled sets need `Labeling.natural(max_labels=2)`. At SNR 266 the
doubles made pyridine / benzene motifs ten-component models and the search
very slow.) Chemically but not
magnetically equivalent label sets (both methyl carbons of isopropyl) keep
separate spins with tied symmetric couplings. Site-specific enrichment that
breaks a fragment symmetry is rejected.

## D37. zulf_hypothesis as a separate package (2026-09-26)

User: keep the hypothesis system apart from the refiner and the models. The
code of D34-D36 moved from `zulf_core.hypothesis` to the top-level package
`zulf_hypothesis`. It uses zulf_core only through public APIs (`refine`,
`RefineSettings`, `MixtureForward` for predictions, spectra, physics) and
never changes them; the neural-model adapter moved from
`zulf_model.evaluation` to `zulf_hypothesis.adapters.model_hints`, which
imports zulf_model lazily. `tests/test_spinsystem.py` checks the boundaries:
zulf_core imports neither zulf_model nor zulf_hypothesis; zulf_hypothesis
imports neither torch nor zulf_model; zulf_model does not import
zulf_hypothesis. The solver's `amplitude_map` (D36) stays in zulf_core as a
general linear amplitude constraint with its own tests. Earlier entries
that name `zulf_core.hypothesis` refer to this package.

## D38. Stage 3: budgeted hypothesis search (2026-09-26)

`zulf_hypothesis.search_hypotheses` refines the top stage-2 proposals (and,
optionally, hinted interpretations) in two variants: free amplitudes and
rates, and abundance ratios with one shared rate. Jobs run in parallel
processes through the public solver API. Every result is scored on one
yardstick (`scoring.py`): the data cores of signal weighting (the same points
for every hypothesis), a per-component complex noise from the quiet parts of
the spectrum, zero-fill decimation, and BIC = chi2 + k ln N with k the free
nonlinear parameters plus amplitude degrees of freedom. The checks run on
each result (the other variant as peer); a new "misfit" finding (reduced
chi2 > 3) and the rate asymmetry trigger extension moves for the best
clean hypotheses. Extensions are built with exchangeable protons kept (the
move hypothesises slow exchange), compared with the same variant of their
parent and accepted only if the BIC improves by 6. The best hypothesis is
the lowest-BIC one without warnings (else the lowest overall). Known
compounds are matched by mapping template fragments onto the hypothesis
fragment (sub-fragments allowed) and translating the couplings. On a
synthetic CH(XH)-CH3 spectrum the loop proposes CH3-CH, extends it by the
coupled proton (Delta BIC about -22000) and ranks the extension first.

## D39. Motif library scan and recognition benchmark (2026-09-26)

Stage-2 group candidates assume isolated X-Hn patterns; coupled CH2 groups
and aromatic rings break that (recognition tests: ethyl, propane-like,
CH2-CH2, benzene and pyridine were missed). `zulf_hypothesis.motifs` adds
whole-motif proposals: registered motifs (methyl, CH-CH3, ethyl, isopropyl,
CH2(CH3)2, CH2-CH2, CH3-C-CH3, benzene ring, pyridine ring; amines with slow
N-H exchange: CH3-NH3+, CH3-NH2, CH3CH2-NH3+, (CH3)2CH-NH3+, (CH3)2CH-NH2;
fast-exchange amines are covered by the carbon motifs because the N-H protons
decouple) are built with
generic couplings (sp3, aromatic, or amines with slow N-H exchange; not
compound-specific), their 1J values
taken from band positions through the X-Hn line ratios; every combination is
evaluated with a closed-form linear solve (free complex amplitude per
isotopologue, background) on the common yardstick and the best per motif
kept. The pipeline lists motif proposals separately (`motif_proposals`);
the search refines the top motifs next to the group proposals.
`zulf_hypothesis.benchmark` simulates registered truth cases with their own
couplings (natural abundance, SNR 60) and reports the rank of the first
proposal with the same topology: the motif scan ranks the truth first in 10
of 12 cases (methyl, CH-CH3, isopropyl, ethyl, CH2(CH3)2, CH2-CH2, benzene,
pyridine, CH3-15NH3+, (CH3)2CH-NH3+); CH3-15NH2 is second behind an
isolated methyl (weak 15N line against the BIC penalty) and CH3-C(O)-CH3
third (practically identical to an isolated methyl). The group path alone
ranks the truth first in 2 of 12. Update (fde3fbb2, ethylenediamine): case
13 "ethylenediamine (fast exchange)" is found as CH2-CH2 at rank 1 (11 of 13);
a slow-exchange H2N-CH2-CH2-NH2 motif was added. Update (4322bdfc,
triethylamine): motif "N-ethyl (Et3N)" (reduced fragment, other N-CH2 protons
as one group) and case 14 "triethylamine (reduced)", found at rank 2 behind
an isolated methyl (weak 13CH2 isotopologue): 11 of 14 at rank 1.

## D40. Quasi-likelihood ranking and a margin for clean hypotheses (2026-09-27)

On real spectra every hypothesis leaves structure (reduced chi2 about 14
for the CH-CH3 samples, about 380 for pyridine), which inflates all chi2
differences. Scores are now chi2 / c_hat + k ln N with c_hat the smallest
reduced chi2 of the evaluated set (at least 1), recomputed after every round;
extensions are accepted on this scale at the end. The best hypothesis is the
lowest score, unless one without warnings lies within `clean_margin` (10),
which is then preferred. Before the change the clean-first rule chose a
benzene ring 5e4 BIC worse than the (warned) pyridine fits, and alanine
accepted a coupled-proton extension that is 1.5 worse on the quasi scale.
The motif quick scan now fixes the isotopologue ratios (one molecule), so a
15N isotopologue cannot take over a 13C band (amine motifs had ranked first
on the lactic-acid and alanine data); the synthetic benchmark is unchanged
(10 / 12 rank 1).

## D41. Extensions start warm and cold (2026-09-27)

Extensions (nested models from moves) are refined from several starts: the
parent's refined couplings (warm), the proposal's original values (cold),
and perturbations of both (`extension_starts` 8, `extension_spread_hz`
1.5; 4 starts were not enough on lactic acid), plus the distinct starts of
the solver's global pattern search over the small couplings and rates with
the one-bond couplings held (`extension_global_search`, 240 s budget). A warm start alone inherits the parent's compensations (on lactic
acid: methyl lines broadened in place of the missing coupling) and stayed in
the parent's basin (the extension was rejected); the cold start had reached
the better basin (chi2 1153 -> 528). Motif and group proposals with the same
component systems are refined once.

## D42. Motif screen, proton-count move, demotion of collapsed fits (2026-09-27)

Triethylamine (4322bdfc) was missed in three ways, each now changed:
- The quick scan picks 1J from band positions with generic couplings; ethyl
  ranked 6th and only the first `top_motifs` (3) were refined. The search
  now gives the first `motif_screen` (8) scan proposals a short refinement
  (fixed-abundance variant when it constrains anything; 3 starts, spread
  2 Hz, continuation 10/3/1/0 1/s, 60 s) and refines the best `top_motifs`
  by score in full; the screen is logged (`step: screen`). A one-start
  screen without continuation chose CH3-NH3+ over the truth on synthetic
  CH-CH3 data (1J starts 2-3 Hz off), hence the broadened continuation.
- New move `change_proton_count`: one more or one fewer proton (1-3) on a
  site and its symmetry orbit, couplings kept (warm start), ring sites and
  exchangeable groups left alone; triggered by misfit, abundance,
  collapsed_component or rate_asymmetry. CH-CH3 -> ethyl is one step.
- New check `collapsed_component` (warn): an isotopologue whose amplitude
  per unit of natural abundance is below 0.2 of the strongest. Fits with a
  finding in `SearchSettings.demote` (default this one) rank after all
  others and cannot be `best`; the table marks them (`demoted`) and their
  delta can be negative. The free CH-CH3 + HX fit on 4322bdfc (CH component
  0.061) had won by 60 quasi-BIC units. D40's clean-margin rule applies
  among the fits that are not demoted. Extensions that failed the acceptance
test rank after the other fits and cannot be `best` either (alanine rerun:
the rejected coupled-proton extension scored 2.3 below its parent, under
`accept_delta` 6, and had become `best`).

## D43. Automatic variants, fit_structure, fit phasing, exact short-record rendering (2026-09-27)

From e66a4b08 (N-ethylmethylamine), where every improvement of the known-
structure fit had been proposed by hand:
- Variant "ratios" (natural ratios, a free rate per isotopologue) next to free
  and fixed; extensions refined in the parent's best variant only.
- `fit_structure(fragment, observed)`: every exchange regime of N/O/S-H
  protons (dropped = fast, kept = slow), all variants, thorough round-0 starts
  (perturbed plus global pattern search), triggered moves, one yardstick; a thin
  entry point over `search_hypotheses`, the solver unchanged.
- Moves: change_proton_count rescales 1J to keep the main X-Hn line (XH at J,
  XH2 at 1.5 J, XH3 at J); free_remote_couplings (couplings built as 0);
  gaussian_line_shape (one Voigt width shared by all components).
- Phasing: the minimum-entropy and per-line model-free estimators both fail on
  dense J multiplets (lines of either sign closer than one width): noise-free
  rendering of the e66a4b08 fit gave single-line phase errors of 10-40 deg,
  and on the real spectrum the delay was 4.7 ms off (reference: complex fits,
  -3.4 to -3.5 ms). The phase-first route therefore phases the spectrum with
  the phase of the best complex fit (all lines modelled jointly) and refits
  the real part; `estimate_phase_lines` stays as a coarse tool for resolved
  spectra with its limits documented.
- Solver speed: models are rendered over `Acquisition.local_record()` (stop +
  h + 1 samples) because the SG baseline is local; retained samples agree
  with the full-record rendering to 4e-12, one e66a4b08 refinement went from
  37.7 s to 10.7 s (same 91 evaluations).
- Addendum (same session): extensions start nested at their parent's fit
  (held couplings stay held in warm starts; new N-H couplings start near
  0.2 Hz; the Gaussian width at 0.02 Hz); from generic starts extensions had
  ended worse than their parents. fit_structure fits the fast-exchange
  structure first and tries the slow-exchange and protonated (ammonium) forms
  as warm-started model moves, variants fixed and ratios, only
  structure-preserving moves (`SearchSettings.moves`). On e66a4b08 the
  protonated form with free remote couplings reached chi2 17456 against 23868
  for fast + remote J, a tie on the quasi-BIC scale (+2.8). The motif scan
  refines 1J of its 8 leading motifs locally (+-4 Hz). Intermediate exchange
  was not implemented (both limits fitted, neither explained the misfit).
- Blind search: motif screen in the ratios variant, two extension rounds,
  extensions compared with the nearest kept ancestor, identical extensions
  not refitted, single-component models fitted as free. Regression over the
  six confirmed samples: right skeleton first in 6 of 6 (one round: 5 of 6).

## D44. Experimental processing as its own package; phase per dataset (2026-09-27)

User: phase correction is data processing, the phase comes from the
instrument and environment, and must be recalibrated for every dataset
(processing parameters such as SG window and crop may differ per spectrum).
- New package `zulf_processing`: raw FID and settings, diagnostics, a
  per-dataset `ProcessingPlan` (crop start after the switching ringing, other
  parameters default for now, every choice with its reason), per-dataset
  phase, `process_dataset`. The processing operator stays single
  (`zulf_core.render.acquisition`); phase convention stays in
  `zulf_core.render.phasing`. The model-free estimators added in D43
  (per-line fits, switching edge) moved there from zulf_core.
- Delay per dataset from the raw FID: the field switch-off edge (half height,
  3.41-3.51 ms in the confirmed datasets) is the start of the zero-field
  evolution; it matches the complex-fit delay within 0.05 ms on three samples,
  lactic acid fits equally well with it, and fitted delays of narrow-band
  spectra absorb model error (pyridine: +1.0 or -3.9 ms by start).
- Zero-order phase per dataset: global search over 0-180 deg together with
  the delay in +-0.5 ms around the edge, then Nelder-Mead fine-tune, with a
  registered criterion (derivative entropy, per-line coherence). Without the
  edge, bands near harmonics make delays about 1 / (2 x line spacing) apart
  equivalent (4 ms on a synthetic J / 2J spectrum). Accuracy on the confirmed
  samples: scripts/validate_processing.py (ANALYSIS_LOG); the criteria are
  the next thing to optimise.
- Addendum (validation): no model-free criterion reaches the calibration.
  Instrument calibration + per-dataset edge (`calibrated_phase`: delay =
  edge + offset, phase0 fixed per instrument/sequence, config
  "phase_calibration") is within 14 deg of each sample's own complex fit at
  125-250 Hz (leave-one-out, 5 samples; ethylenediamine reference itself
  ambiguous); entropy / line criteria err 2-86 deg. So: delay measured per
  dataset, phase0 from the instrument calibration by default; the criteria
  stay available (`phase_criterion`) and are validated by the same script.
  The calibration must be refreshed when the sequence or hardware changes.
- Addendum (speed): `diagnose_raw` skips the multi-exponential fits of
  `zulf_core.diagnose_fid` (546 s per 64k-point FID, unused here); processing
  of one dataset now takes 0.1 s.

## D45. Chemical exchange as a Liouville-space model with a fitted rate (2026-09-29)

Exchangeable protons (N-H, O-H) were modelled only in the two limits (fast:
dropped; slow: static spins, D43). `zulf_core.physics.exchange` adds the
intermediate regime: each exchanging spin x is replaced by an uncorrelated
solvent spin at rate k_x (1/s), d delta/dt = -2 pi i [H, delta] + sum_x k_x
(Tr_x(delta) (x) 1_x / d_x + eps_x I_z,x - delta), high-temperature and
linear in polarization; eps_x is the incoming spin's deviation weight
(default: its preparation weight, the solvent was prepolarized alike and its
polarization is static at zero field; 0 for an unpolarized pool).

- Restriction: zero-field H and the exchange map commute with rotations and
  conserve M, and preparation, source and detection are rank-1, q = 0, so
  the dynamics live in the rank-1 subspace of q = 0 operators (Casimir
  eigenvalue 2): 1001 dimensions for 7 spin-1/2 (3432 for 8) instead of 4^n.
  The source term is carried by an augmented generator [[L, s], [0, 0]].
- Output: a TransitionList whose lines carry their own decay rates
  (`line_rates`, added to the fitted component rate by both renderers);
  non-oscillating modes are not rendered (listed in the metadata).
- Derivatives: analytic, divided differences in the eigenbasis of the
  augmented generator (X = V^-1 dL V; pairs inside a cluster of equal
  eigenvalues give t exp(lambda t) terms), for couplings and log k.
- Solver: parameter kind `log_exchange` (`Parameterization.add_exchange`,
  tied across isotopologues), MixtureForward switches the component to the
  exchange model when it has such a parameter; `scripts/fit_staged.py
  --exchange LABEL`.
- Tests (tests/test_exchange.py): brute-force propagation in the full
  Liouville space (explicit partial traces, affine expm) at k = 0 ... 400
  1/s and for an unpolarized solvent; k = 0 equals compute_transitions;
  k = 1e5 decouples the spin; line derivatives vs central differences on a
  system with degenerate modes; solver Jacobian vs residual differences;
  refinement recovers k.
- Optional GPU: the dense steps (Casimir basis, eigen-decomposition,
  derivative products) run in NumPy or PyTorch (`exchange.set_backend`,
  `ZULF_LINALG_DEVICE=cuda`, `--device cuda`); torch on CPU agrees with
  NumPy (test). Cost on CPU for 7 spins: about 3 s per evaluation, 65 % in
  the 1001-dimensional eigen-decomposition.
- Not covered: pulses and static fields with exchange, exchange between two
  sites of one molecule, 14N quadrupolar relaxation.

## D46. Known-structure complex fits: rate families, component search, multi-species, model-line passes (2026-10-05)

The amine fits (docs/analysis/2026-10-04_amines_complex-fit.md, 2026-10-05_amines_overnight.md) settled these
choices in `scripts/fit_joint_series.py`.

- **Decay-rate families.** One decay rate per line cluster (`--family-edges`), not one per isotopologue.
  - A single rate cannot carry the narrow and the broad lines of one component: stage 1 objectives were 0.10-0.31,
    against 0.038-0.072 with families.
  - Edges are chosen by hand from the line clusters and are reported with every run as nuisance parameters.
  - A sharp line that shares a family with broad lines gets its own narrow family:
    - triethylamine 0.049 -> 0.038;
    - ethylenediamine 0.0066 -> 0.0051.
    In both, the 1J values moved by < 0.07 Hz.
  - A rate pushed to its bound or a 15x rate spread inside one species is a warning sign. The family may be hiding
    a predicted line the data do not have.
- **Component search on by default** (`--component-search both`).
  - Every isotopologue is searched on the window it dominates, from random starts of its own couplings.
  - Candidates are accepted by the whole objective, never by the window cost: window fits re-solve gains and
    phase, so they overstate the gain.
  - The start pass guides the multi-start; the end pass is followed by a global refit that is kept only if it
    wins.
  - Reason: the N-ethylmethylamine 13CH2 basin. The multi-start never reached it; it was -33 % on the whole
    objective.
- **Several molecules** (`--structure` as a JSON list, `combine_models`): ratios free between parts, fixed within
  each.
  - Seeds must place a new species away from the existing one, because a copy has zero gain and no gradient.
  - A second species is a numerical result until it is identified chemically.
- **Model-line passes** (`--model-line-passes N`, `MixtureForward.add_signal_cores`) are **off by default**.
  - Signal weighting is built from data peaks, so a model line where the data show none is weighted only 0.2.
  - The passes add the model's own line envelope to the cores and refit, as `refine`'s `signal_model_passes`
    does.
  - They are off by default so that older runs stay comparable. The record keeps the objective under the
    original weights for comparison.
- **N-H protons decoupled (fast exchange) by default for amines.** The evidence:
  - no natural-abundance 15N-H lines;
  - static N-H models are worse and fit their N-H couplings to about 0;
  - a fitted exchange rate goes fast.
  The exact exchange model (D45) does not scale to ethylenediamine with N-H kept (9 spins, no reducible group,
  35 GiB), so intermediate exchange there needs an approximation.

## D47. Fitted static field, default zero field (2026-10-06)

Some measurements are made in a known nonzero field (the acetonitrile run of 2026-09-28 is a fringe-field
measurement: its 2J band is split into 271.3 / 273.8 Hz, which a longitudinal field of about 0.08 uT reproduces,
while the zero-field 13CH3 model leaves a relative residual of 0.73). The field is therefore an optional fitted
parameter (`ParameterPolicy.fit_field`, default False: zero field, unchanged results).

- With preparation and detection along z the signal depends only on the transverse size B_transverse and on |Bz|
  (rotation about z, and a pi rotation about x, leave it unchanged). Two parameters `field_transverse_ut` (placed in
  Bx, By = 0, so the real eigensolver stays in use) and `field_z_ut`, both >= 0 (`field_bounds_ut`).
- The signal is stationary at zero field (zero gradient), so starts are nonzero (`initial_field_ut`, default
  0.02 uT each); a fit that should be able to reach zero field still can (lower bound 0).
- `MixtureForward.protocol_for(values)` puts the field into the protocol of every evaluation (transitions,
  coupling derivatives, exchange); the field columns of the Jacobian are central differences (the analytic
  derivative needs the Zeeman directions in physics.derivatives; possible later).
- The global pattern search (solver.search) still uses the fixed protocol: with a fitted field, rely on the
  multistart from the field starts.
- Test: tests/test_solver.py FieldFitTests (data from brute-force propagation in a field; recovery of
  B_transverse and Bz from a nonzero start; the zero-field model is at least 10x worse on such data).
- D19 still holds for reporting: a fitted field is reported with every result that uses it.
- Names: the transverse component was first called `field_perp_ut` (axis `perp`); renamed the same day to
  `field_transverse_ut` (axis `transverse`) on Xuehan's request. Figures write it as B with a perpendicular sign.

## D48. ZULF Studio: one session for the window, the fitter and AI agents (2026-10-06)

Xuehan asked for a standalone real-time simulator with a Python desktop UI (PySide), fitting with the
fit-progress slider, log and terminal inside, and an API for AI. Design:

- `zulf_studio.session.StudioSession` holds the whole state and every operation, without a GUI. The model is
  rebuilt from the structure specification on every change (zulf_hypothesis.structure_spec; build_model takes
  well under 1 ms) and the transitions are exact in the session's static field (D47). Triethylamine updates in
  about 50 ms, so sliders are live.
- The displayed simulation is a sum of complex Lorentzians with one decay rate: a quick look. The fit is
  `scripts/fit_joint_series.py` run as a subprocess from the current couplings and field, so the fitter, its
  processing operator and its objective are the same code as in the analyses (no second implementation); its
  trace frames (`--trace`) feed the progress slider and are drawn as computed.
- `zulf_studio.api` exposes the same operations as JSON tools over local HTTP (127.0.0.1, /api/tools with
  schemas, Anthropic and OpenAI formats). The window listens to the session, so a person and an agent work on
  one state and the window follows the agent.
- Log (session, api, fit, terminal sources; also runs/studio/studio.log), a shell terminal (repository
  directory, the studio's Python first on PATH) and a Python console with the session are tabs of the window.
- `structure_for` / `override_couplings` moved from scripts into `zulf_hypothesis.structure_spec` (the scripts
  re-export them), so a package does not import scripts.
- Tests: tests/test_studio.py (closed forms: XA3 lines at J and 2J, a 13C-1H line split by (gH + gC) B / 2 in a
  transverse field, Lorentzian FWHM = rate / pi; fit command, fit application and trace frames; HTTP API;
  offscreen window and session following each other across threads).

## D49. Studio display scale and figures (2026-10-06)

Two requests of Xuehan (handed over from a read-only session):

- The simulated curve jumped while 1J was dragged (acetonitrile, 136.8-136.9 Hz): the display scale was a
  least-squares fit on the shown real part, which faded to 3e-6 as the line moved through a dispersive data
  feature and then fell back to max|data| / max|sim| (a 300x jump). The scale is now the least-squares fit on the
  magnitudes (positive, continuous in the parameters, the same for every shown part); `lock_scale` freezes it.
  The status bar shows it. Test: DisplayScaleTests (dispersive synthetic data, scan of J).
- Figures: Studio runs scripts/paper_figure.py as a subprocess (PNG at 300 dpi, PDF, SVG, caption.txt; preview
  in the window; export to a chosen directory; API tools figure_command, make_figure, figure_status,
  export_figure). If the parameters are exactly an applied fit, its run is drawn (its options and fit.json);
  after the sliders moved, a fit.json-compatible parameter file of the current state is written and the figure
  and caption say "manual parameters (not a fit)".
- paper_figure.py no longer assumes 4 kHz: the sampling rate comes from the series entry's record (or
  --sampling-rate), the display zero-order phase from the series entry's phasing (or --phase0-deg; the 4 kHz
  calibration only as the fallback); the caption names the fitted field and the sampling rate; --formats, --dpi,
  --manual. scripts/regression_confirmed.py is unchanged.
- make_series_entry.py records the FID as "source_fid" (not "fid": fit_joint_series reads an entry with "fid"
  as a raw FID to process with the 4 kHz recipe). Studio asks for the FID when a series has none.
- Tests: FigureTests (a 136 Hz line of a 2 kHz FID appears at 136 Hz in the display spectrum; a full Studio
  figure of a synthetic 2 kHz series, files and caption).
- Auto phase of the shown data (`auto_phase`, button in the data panel): "model" finds the phase0 and delay
  that best match the current simulation (delay grid +-3 ms, phase0 in closed form); "data" is model-free
  (zulf_core.render.phasing.estimate_phase), delay within +-0.5 ms by default because for lines at J and 2J
  delays about 1 / (2 J) apart are equivalent for it (a +-3 ms search landed 3.5 ms off on synthetic data).
  Acetonitrile: model -17.7 deg, 0.365 ms on top of the series phasing (the 4 kHz calibration). Test:
  AutoPhaseTests (known phase and delay).
- Fixed on the way: the quick-look Lorentzian had the conjugate sign of the processed spectra (numpy FFT
  convention: a gamma / (gamma + i (f - f_k))); the real part was unaffected, the imaginary part and any complex
  match were wrong.

## D50. AI assistant in ZULF Studio: Anthropic and OpenAI models drive the session (2026-10-06)

Xuehan asked for Anthropic and Codex (OpenAI) API access in Studio. `zulf_studio.assistant.StudioAssistant` runs
a tool loop: the model gets the Studio tools (the same definitions as /api/tools, executed in-process on the
session the window shows), every call and result goes to the transcript and the log (source "ai"), arrays of
`simulate` are summarised before they reach the model, and a system prompt says that fitted values are
conditional numerical results.

- Anthropic: Messages API through the `anthropic` SDK, manual loop with client tools, default model
  claude-opus-5-5, server-side refusal fallback on (`fallbacks="default"`); thinking blocks are passed back
  unchanged (append-only history).
- OpenAI: Responses API through the `openai` SDK with function calling; no default model (names of Codex models
  change): the model field or OPENAI_MODEL.
- Credentials from the environment (ANTHROPIC_API_KEY / ANTHROPIC_AUTH_TOKEN / `ant auth login`,
  OPENAI_API_KEY), or entered in the AI tab: kept in memory, or with "Remember" stored encrypted in the macOS
  Keychain (service zulf-studio, zulf_studio/credentials.py; Xuehan: enter once) and loaded at start, an
  environment value winning. Never in a file, config, the repository or the log. The terminal and the Python console are not tools of the model.
- Environment: anthropic 1.11 and openai 3.24 added to the `zulf` env (15 packages, all new; no version
  changes); pyproject extra "ai". User guide docs/STUDIO.md, environment docs/ENVIRONMENT.md.
- Tests: AssistantTests (fake provider clients, no network: tool calls reach the session, unknown tools come
  back as errors, simulate results are summarised).

## D51. Reduced coupling constants K next to J (2026-10-07)

Xuehan asked for the coupling with the gyromagnetic ratios removed. K = 4 pi^2 J / (h gamma_A gamma_B)
(docs/CONVENTIONS.md) is now:
- written by fit_joint_series for every coupling: fit.json couplings gain "nuclei" and "K_at_x" (unit in
  "reduced_coupling_unit": 1e19 N A^-2 m^-3); J_table.csv gains the columns "nuclei" and "K at x" at the end
  (earlier columns unchanged);
- available as `zulf_core.nuclei.reduced_coupling`, `coupling_from_reduced`, `convert_coupling` for starts and
  priors of other isotopes: 15N literature couplings -> 14N (the 14N model, PLAN Phase 3e), 1H -> 2H (labelled
  samples). The primary isotope effect on K (about 1 %) is neglected; converted values are not results.
- Using K in the J-to-structure likelihood (Phase 7) is open for discussion.
Tests: tests/test_spinsystem.py ReducedCouplingTests (CODATA gammas, independent formula); tests/test_j_tuner.py
RunToolsTests (K and nuclei in fit.json for C-H and H-H couplings).

## D52. Route B input mode: J (default) or K (2026-10-07)

`zulf_model.structure.edge_model` gains an input mode (`EdgeModel.mode`, saved with the model; files without it
load as "J"). Mode "K" feeds the equivalent 13C-1H coupling at equal K (13C-1H rows unchanged, other pairs
rescaled by their gyromagnetic ratios); per-row isotopes can be given (`features(..., nuclei=...)`,
`accuracy(..., nuclei=...)`). Measured (docs/analysis/2026-10-07_j-to-structure_j-vs-k.md): equal performance on
13C / 1H data, and K alone transfers to deuterated couplings (0.843 against 0.298 for J). J stays the default;
K is to be used once labelled-isotope or N-coupling observations enter (decision with Xuehan). Test:
tests/test_j_structure.py ReducedCouplingModeTests.

## D53. Fitted gyromagnetic ratio in a known field: naming the heteronucleus from the data (2026-10-07)

In a blind test the isotopes may be unknown. Zero-field line positions depend on J only, so a 13C line and a
15N line with another J can sit at the same place; a field splits lines by amounts proportional to the
gyromagnetic ratios of the coupled nuclei (a 13C-1H line by (g_H + g_C) B / 2 in a transverse field), so in a
known field the data fix gamma.
- `Protocol.gamma_overrides` (default empty): gamma per nucleus for the field term and the gamma weights;
  `zeeman()` takes the protocol's gamma; transitions, derivatives and the brute-force reference use it.
- `ParameterPolicy.fit_gamma` / `initial_gamma` / `gamma_bounds_hz_per_ut` (parameters `gamma_<symbol>`, central-
  difference Jacobian) and `field_fixed` (hold the field at `initial_field_ut`: the known field). gamma and a free
  field trade against each other, so fit gamma only with a known field (from a reference sample, or from lines of
  known nuclei).
- fit_joint_series `--field transverse,z` (known, fixed) and `--fit-gamma 13C [--gamma-start ...]`; fit.json
  `gamma_identification` lists the nearest registered nuclei of the same spin.
- Spin-1/2 heteronuclei (13C, 15N, 19F, 31P) are told apart this way; a change of spin (1H vs 2H, 15N vs 14N) is a
  different spin system and is decided by comparing hypotheses (labelling option of the search).
- Measured: synthetic 13C-1H and 15N-1H pairs in 0.1 uT (brute-force data), a 13C-labelled model with gamma free
  ends at the true gamma (13C from a start of 8.0; 15N, -4.316, from -3.0). Acetonitrile with the field held at
  its fitted value: gamma 11.04 from a start of 8.0, nearest 13C (3 %; 31P 56 %, 15N 139 %).
Tests: tests/test_solver.py GammaFitTests, tests/test_j_tuner.py GammaOptionTests; physics, search, derivative,
exchange and solver suites unchanged (70 tests).

## D54. Labelling hypotheses for samples of unknown labelling (2026-10-07)

Blind samples may be enriched (15N) or in D2O (exchangeable protons as 2H). Labelling is part of the hypothesis
and the data rank it:
- `ProtonGroup.isotope` ("1H" default, "2H"); the builder builds proton groups with their isotope. Fragments
  without it serialise as before.
- `labeling.deuterate_exchangeable(fragment)`: exchangeable groups (flagged, or on N / O / S, the exchange
  variants' rule) become 2H, their couplings converted at equal K (D51).
- `labeling.labeling_variants(fragment, mode)`, modes natural, 15N (N sites at 0.98, up to two labels),
  2H-exchange, unknown (every one that applies). `fit.fit_structure_labelings` fits each with fit_structure and
  ranks the best fits on one scale (chi2 over the common overdispersion, plus the criterion penalty);
  analyze_sample.py `--labeling`, labelings.json, one report per labelling.
- Limit: in fast exchange the exchangeable protons are decoupled whether 1H or 2H, so those two hypotheses tie.
- Measured on synthetic methylamine (slow exchange, three truths): natural, 15N and 2H-exchange each recovered,
  the runner-up worse by at least 8.5e4 in BIC. Test: tests/test_hypothesis.py LabelingHypothesisTests (15N truth,
  about 3.5 min).

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
