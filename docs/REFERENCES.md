# References and prior art

Survey date: 2026-09-25. Status tags: [V] title/authors/venue confirmed by
search, [V-code] repository README read, [U] unverified detail. Full texts were
mostly not accessible from the survey environment; verify [U] items before
relying on them. No published machine-learning model for ZULF spectrum to spin
system inversion was found (confirm with a Scholar search).

## ZULF J-spectroscopy

| Item | Link | Use here |
| --- | --- | --- |
| Barskiy et al. 2025, Prog. NMR Spectrosc. 101558, review [V] | https://doi.org/10.1016/j.pnmrs.2025.101558 | Main physics reference. |
| Butler et al. 2013, J. Chem. Phys. 138, 184202, zero-field multiplets [V] | https://doi.org/10.1063/1.4803144 | Strong couplings truncate weak ones; motivates collective sectors and coarse-to-fine decoding (1J first). |
| Theis et al. 2013, Chem. Phys. Lett. 580, 160, multiplet rules [V] | https://www.sciencedirect.com/science/article/abs/pii/S0009261413008191 | XA_n line rules (J, 3J/2, J and 2J) used as simulator unit tests. |
| Blanchard et al. 2013, JACS 135, 3607, aromatics [V] | https://doi.org/10.1021/ja312239v | Aromatic benchmark spectra. |
| Wilzewski et al. 2017, JMR 284, 66, sub-mHz J fitting [V] | https://doi.org/10.1016/j.jmr.2017.08.016 | Closest prior art for the solver; precision target. |
| Stern and Sheberstov 2023, Magn. Reson. 4, 87, simulation A to Z [V] | https://mr.copernicus.org/articles/4/87/2023/ | Independent reference for detection and polarization conventions. |
| Tayler et al. 2017, Rev. Sci. Instrum. 88, 091101 [V] | https://pubs.aip.org/aip/rsi/article/88/9/091101/950953 | Noise and linewidth priors. |
| Barskiy et al. 2019, Nat. Commun., exchange [V]; urea 2021 JPCL [V] | https://www.nature.com/articles/s41467-019-10787-9 | Justifies dropping exchangeable NH/OH protons; 14N must be exchange-decoupled for spin-1/2 scope. |
| Andrews, Liu, Zumbrunn et al. 2026, arXiv 2604.26071, natural-abundance 13C ZULF with DFT [V abstract] | https://arxiv.org/abs/2604.26071 | First candidate real test set (13 molecules, isotopologue mixtures). |
| Put et al. 2021 Anal. Chem.; Put et al. 2023 Commun. Chem. [V] | https://doi.org/10.1021/acs.analchem.0c04738 | Labelled 13C and natural-abundance 15N cases. |

## Machine learning for NMR inverse problems

| Item | Link | Use here |
| --- | --- | --- |
| Alberts, Zipoli, Vaucher 2023, spectra to SMILES transformer [V] | https://chemrxiv.org/doi/10.26434/chemrxiv-2023-8wxcz | Seq2seq, beam search and top-k evaluation template. |
| Hu et al. 2024, ACS Cent. Sci. 10, 2162 [V] | https://doi.org/10.1021/acscentsci.4c01132 | Multitask heads (composition plus structure). |
| Jonas 2019 NeurIPS, deep imitation for inverse problems [V] | https://proceedings.neurips.cc/paper/2019/hash/b0bef4c9a6e50d43880191492d4fc827-Abstract.html | Alternative graph-building decoder. |
| Sridharan et al. 2022 JPCL, DeepSPInN (MCTS) [V] | https://pubs.acs.org/doi/10.1021/acs.jpclett.2c00624 | Forward-model-guided search alternative to beam search. |
| NMR-Solver 2026 Nat. Commun. [V] | https://www.nature.com/articles/s41467-026-71315-0 | Retrieve-then-refine design, as here. |
| UltraNMR arXiv 2606.20756; NMRGym arXiv 2601.15763 [V] | - | Simulation pretraining then adaptation; scaffold splits. |
| Lemm, von Rudorff, von Lilienfeld 2024, Digital Discovery 3, 136 [V] | https://pubs.rsc.org/en/content/articlelanding/2024/dd/d3dd00132f | Forward-model error limits discrimination; informs J bin width and tolerances. |
| IMPRESSION / IMPRESSION-G2 (Chem. Sci. 2020, 2025) [V] | https://pubs.rsc.org/en/content/articlelanding/2020/sc/c9sc03854j | Near-DFT J for a realistic-J test distribution. |
| Kaggle CHAMPS scalar couplings; Bratholm et al. 2021 PLOS ONE [V] | https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0253612 | DFT J distributions to calibrate generator rules. |
| ANATOLIA total lineshape fitting [V] | https://github.com/dcheshkov/ANATOLIA | Solver benchmark. |
| Physics-informed networks for ABC/ABCD couplings 2026 [V, authors U] | https://www.sciencedirect.com/science/article/pii/S2949747726000217 | Network guess plus QM refinement for fixed topology; we generalize to unknown topology. |
| MolDeTr 2026 Anal. Chem. [V] | https://github.com/smidooo/MolDeTr | DETR-style set prediction on 1D spectra; precedent for the set baseline. |
| DEEP Picker 2021 Nat. Commun. [V] | https://github.com/lidawei1975/deep | Synthetic-only training transferring to real spectra. |

## Inference, sets, generation, canonical forms

| Item | Link | Use here |
| --- | --- | --- |
| Cranmer, Brehmer, Louppe 2020 PNAS; sbi toolkit [V] | https://github.com/sbi-dev/sbi | Possible NPE head for continuous J given topology; SBC calibration checks. |
| DETR (Carion 2020), Deep Sets, Set Transformer [V] | - | Hungarian set loss over components. |
| Surge structure generator [V] | https://github.com/StructureGenerator/surge | Exhaustive small-molecule enumeration for coverage and hard negatives. |
| GuacaMol, RDKit [V] | https://github.com/BenevolentAI/guacamol | Distribution checks for procedural graphs. |
| nauty / pynauty [V] | https://github.com/pdobsan/pynauty | Exact canonical certificates for spin graphs (optional dependency). |

## Software

| Tool | Link | Note |
| --- | --- | --- |
| ZULFPy (Ajoy lab, GPL-2.0) [V-code] | https://github.com/Ajoy-Lab/ZULFPy_release | Simulator cross-check, NMRduino I/O, work-in-progress gradient J fitting. GPL: compare, do not copy code. |
| Spinach [V] | https://spindynamics.org | Gold-standard check (MATLAB). |
| SpinDynamica [V] | https://www.spindynamica.soton.ac.uk | Symbolic checks. |
| nmrsim [V] | https://github.com/sametz/nmrsim | High-field only; not used. |
| MRSimulator [V] | https://github.com/deepanshs/mrsimulator | Weak coupling only; unsuitable for ZULF. |
| QuTiP [V] | https://qutip.org | Relaxation and exchange extensions. |
| LSD / pyLSD, Sherlock [V] | https://nuzillard.github.io/LSD/ | Later J-to-structure stage. |
| NMRShiftDB2, NP-MRD [V] | https://nmrshiftdb.nmr.uni-koeln.de | Experimental J distributions. |

## Coupling values (generator priors; signs as observed J)

| Coupling | Values | Status |
| --- | --- | --- |
| 1J(13C,1H) | about 500 s(C): sp3 about 125, sp2 160-165, sp about 250 Hz; plus 20-30 Hz with electronegative substituents; benzene 158.3 Hz; positive | [V snippet] |
| 2J(H,H) | sp3 CH2 -10 to -15 Hz; sp2 =CH2 0 to +3 Hz | [V snippet] |
| 3J(H,H) | Karplus 0 to about 13 Hz; free rotation averages about 7 Hz; Haasnoot-de Leeuw-Altona 1980 | [V], parameters [U] |
| 2J(13C,1H) | about -5 to +5 Hz; sp3 usually negative; benzene about +1.1 Hz | [U] |
| 3J(13C,1H) | 0 to about 10 Hz, Karplus-like; aromatic 7-8 Hz | [U] |
| 1J(13C,13C) | single 35-45, aromatic about 56, double about 68, triple about 172 Hz | ranges [V], exact [U] |
| 1J(15N,1H) | negative; NH4+ -73.5; sp3 amines about -61 to -67 [U]; anilines -78 to -90; amides and pyrroles -89 to -97 Hz | partly [V] |
| 2J/3J(15N,1H) | usually magnitude up to 5 Hz, variable sign; azines 2J up to about -16 Hz | [V title] |
| 1J(13C,15N) | small, usually negative; amines -4 to -6.5, amides -13 to -15, nitriles about -17 Hz | partly [U] |

## Design consequences adopted

1. Global J sign is unobservable for the real protocol (H to -H leaves
   frequencies and amplitudes unchanged). Matching and evaluation compare up to
   one global sign flip; canonical forms fix the largest-magnitude observable
   heteronuclear coupling positive. Relative signs remain observable.
   (`docs/DECISIONS.md` D12.)
2. Components with only one nucleus type give no zero-field signal; the
   generator rejects components without an observable heteronuclear coupling.
3. Theis and Butler line rules are simulator unit tests (`tests/test_physics.py`).
4. The set baseline uses Hungarian matching over components (DETR style).
5. Splits are by spin-graph topology (scaffold-style), never by render.
6. The rule-based J table (`configs/couplings_v1.json`) is a first prior;
   calibration against DFT J data (CHAMPS, IMPRESSION-G2) is a Phase 1 task.
7. First real benchmark candidates: Andrews et al. 2026 molecules and the
   local isopropylamine data; ZULFPy is a simulator and solver comparison.

## Second survey: optimizations for the current bottlenecks (2026-09-25)

Bottlenecks from the CPU verification runs: (a) J-token accuracy within one
bin about 12 percent; (b) weak second components missed; (c) 13C versus 15N
and proton-count confusion; (d) top-1 far below top-10; (e) rendering
throughput. Tags as above; arxiv, PubMed and ScienceDirect pages could not be
opened from the environment, so [V] here means confirmed from search snippets
or opened GitHub pages.

| Item | Link | How to apply |
| --- | --- | --- |
| Farebrother et al. 2024 ICML, "Stop Regressing" (HL-Gauss) [V] | https://arxiv.org/abs/2403.03950 | Our soft neighbour-bin targets are close already; decode J by the local expectation over bin probabilities instead of the argmax, sigma relative to the nonuniform bin width (a). |
| Imani, White et al. 2026 JMLR, histogram loss [V] | https://github.com/marthawhite/Histogram_loss | Gains come from optimization; keep classification for J (a). |
| Song and Bahri 2025 TMLR, decoding-based regression [V] | https://arxiv.org/abs/2501.19383 | Coarse-to-fine J: a coarse token (about 2 Hz) then a fine token (about 0.1 Hz), smaller better-populated softmaxes (a). |
| Zausinger et al. 2025 ICML, Number Token Loss [V] | https://github.com/tum-ai/number-token-loss | Add a Wasserstein term over J-bin centres to the token loss (a). |
| Golkar et al. 2023, xVal [V] | https://arxiv.org/abs/2310.02989 | Continuous embedding of already decoded J as decoder input; weaker as an output head. |
| Li, He et al. 2024 NeurIPS, per-token diffusion loss [V] | https://arxiv.org/abs/2406.11838 | Heavier continuous J head; low priority. |
| PANPE 2024, reflectorch (MIT, GPU simulator) [V] | https://github.com/schreiber-lab/reflectorch | Network proposals plus likelihood refinement of many candidates; template for on-GPU simulation in training (d, e). |
| SpecCal 2026 [V] | https://arxiv.org/abs/2607.27788 | Rerank beam candidates by solver residual after a short fit; add mutation neighbours (13C/15N swap, group size +/-1) (c, d). |
| Casanovo 2024 Nat. Commun. (Apache-2.0) [V] | https://github.com/Noble-Lab/casanovo | Cheap physical consistency filters inside beam search (d). |
| BayesFlow [V] | https://github.com/bayesflow-org | Calibrated structure classifier as reranker (d). |
| DreaMS 2025 Nat. Biotechnol. (MIT) [V] | https://github.com/pluskal-lab/DreaMS | Peak tokens with high-resolution Fourier frequency features and pairwise-difference attention; differences encode J in ZULF (a). |
| NMRTrans 2026 [V] | https://arxiv.org/abs/2602.10158 | Peak-set input for NMR performs well on experimental spectra (a). |
| Pix2Seq 2022 ICLR [V] | https://arxiv.org/abs/2109.10852 | Fake-component sequence augmentation and EOS down-weighting to raise recall of weak components (b). |
| NMRQNet 2023; Raman iterative unmixing 2025; DeepRaman 2023 [V] | https://github.com/XiaqiongFan/DeepRaman | Explain-away: fit the dominant component, subtract, re-run for the next (b). |
| Asymmetric Loss, Ridnik et al. 2021 ICCV [V] | - | Existence and weight tokens; oversample 1-10 percent weights (b). |
| Diverse Beam Search, Vijayakumar et al. AAAI 2018 [V] | https://arxiv.org/abs/1610.02424 | Group beams by structure prefix so top-10 differ in structure, not only in J (c, d). |
| Scheduled sampling for transformers 2019 [V] | - | Two-pass decoding to reduce J-token error compounding (a). |
| JMR 368:107792 2024, neural nets for strongly coupled spectra [V] | - | Simulated-only training; real-data obstacles reported. |
| PINN for ABC/ABCD J extraction 2026 [V] | https://www.sciencedirect.com/science/article/pii/S2949747726000217 | Auxiliary losses from multiplet rules; peak-list input; physics refinement (a). |
| 2D-JCOG 2026 Anal. Chem. [V] | - | Auxiliary head for which peak pairs share a J (a, c). |
| Zaleski and Prozument 2018 J. Chem. Phys. 149, 104106 [V] | - | Rotational spectra: classify Hamiltonian type, then regress constants; supports structure-then-J. |
| Stern and Sheberstov 2023 [V] | https://mr.copernicus.org/articles/4/87/2023/ | Block-diagonalize by total F_z (adopted: `physics.transitions` M blocks). |
| PyTorch issue 175585 [V] | - | Batched CUDA eigh is fast up to n = 32 and slows sharply above; keep blocks at or below 32 for a GPU renderer (e). |

No dedicated ML-for-ZULF paper from 2022-2026 was found; a DDNP-ZULF
perspective (arXiv 2307.06973) proposes ML for the inverse problem [V].

Ranked plan (gain per effort): (1) rerank top-k by short solver fits plus
mutation neighbours; (2) GPU batched renderer on F_z blocks; (3) J decoding:
local expectation, Wasserstein term, coarse-to-fine tokens; (4) explain-away
second pass and Pix2Seq-style augmentation for weak components; (5) peak-token
stream with pairwise frequency-difference attention and diverse beam search.
