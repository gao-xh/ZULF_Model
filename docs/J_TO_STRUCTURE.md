# J to structure (Phase 7)

Question (Ashok, 2026-10-05): the J couplings carry far more information than chemical shifts (n^2 values
instead of n), so in principle they should say which spins are connected, and the structure could be read
from the J couplings alone. This document records how the toolkit does that and where the limits are.

Three routes are kept side by side. None replaces another.

| Route | Input | Output | Code |
|---|---|---|---|
| Hypothesis pipeline (kept) | the spectrum | ranked fragments fitted to the data | `zulf_hypothesis.search`, `fit_structure` |
| A: rule-based direct readout | a fitted J network | ranked heavy-atom graphs | `zulf_hypothesis.j_structure` |
| B: learned likelihood | a fitted J network | ranked heavy-atom graphs | `zulf_model.structure.edge_model` |

The hypothesis pipeline starts from the spectrum and fits structures. Routes A and B start from the couplings
a fit (or a person) has already read off, and ask which bonding graph explains them.

## Observation

`configs/j_networks/*.json` (format of `JObservation.from_dict`):

- `units`: one entry per 13C site, with the number of protons on the carbon (from the 1J pattern) and the
  number of symmetry-equivalent copies (e.g. the two CH3 of isopropylamine are one unit with two copies).
- `protons`: proton groups, each attached to one copy of a unit.
- `couplings`: J(unit, group) from the 13C of copy 0 of the unit, and J(group, group).
- `sigma`: optional per-coupling uncertainty; the default is 0.1 Hz.
- `truth`: optional, a bond list for checking. Bonds are written `-`, `=` or `#`; unseen atoms are `X1`, `X2`.

13C-13C couplings are not used, since they need doubly labelled molecules.

## Candidates

Candidates are connected graphs over the carbon copies plus up to `max_unseen` (default 2) unseen heavy atoms X.

- **Unseen atoms.** An X is N, O, or a carbon without protons. None of these has a 13C line of its own carrying
  a 1J, so the fit never sees them directly.
- **Bonds.** Carbon degree is at most 4 minus its protons. X degree is at most 3, its valence at most 4, and two
  X atoms are never bonded. At most one ring is allowed.
- **Bond orders.** Every assignment of up to two double or triple bonds within the valences.
- **Symmetry.** The copies of a unit must stay equivalent (equal Weisfeiler-Lehman colours, with bond orders
  included). Distinct units must stay distinct.
- **Duplicates.** Isomorphic labellings are merged, keeping the one that explains the data best.

## Score

score = sum over couplings of log p(J | nuclei, bond count) + log prior.

- **Bond count.** For J(C, H) it is the C-C distance to the proton's carbon plus 1. For J(H, H) it is the
  distance between the two carbons plus 2.
- **Bond orders.** They do not change bond counts. They enter through the 1J of each carbon, whose size depends
  on hybridization:
  - sp3: about 118-150 Hz (including the shift from electronegative neighbours);
  - sp2: about 145-180 Hz;
  - sp: about 250 Hz.
- **Prior.** Each unseen atom costs -1, each ring -2, each double bond -1.5, each triple bond -3, and each
  unfilled carbon valence -3.

The two routes differ only in p(J | ...):

- **Route A.** Literature-style ranges (`couplings_likelihood`, `one_bond_likelihood`), each broadened by the
  fit sigma and 0.3 Hz, plus a 2% wide tail.
  - 1J(CH): one range per hybridization (as above).
  - 2J(CH): -7.5 to 4 Hz. 3J(CH): 0 to 9 Hz.
  - 2J(HH): -17 to -6 Hz. 3J(HH): 4 to 9.5 Hz, broadening to 0 to 18 Hz.
  - Long range: a spike at 0 plus a Laplace tail with a 0.4 Hz scale.
- **Route B.** Two MLPs trained on generator pairs (`observation_from_graph`, rules
  `configs/couplings_v1.json`):
  - p(bond count | J, kind, proton counts of the two carbons);
  - p(hybridization | 1J, protons).

  Ranking uses p(class | J) / p(class), the ratio to the class frequency in training. By Bayes this is
  p(J | class) up to a factor that does not depend on the candidate.

## Results (2026-10-05)

`python scripts/j_structure_benchmark.py --samples 100 --model runs/models/j_edges_v1.json`

Route B was trained on 4000 graphs (48436 couplings) in about 20 s on one CPU thread. Its held-out bond-count
accuracy is 0.85.

| Set | Route A | Route B |
|---|---|---|
| synthetic, 2-5 heavy atoms, at most 5 carbon copies (100) | top-1 0.91, top-3 1.00 | top-1 0.91, top-3 1.00 |
| N-ethylmethylamine (fast fit) | rank 1, margin 13.4 | rank 1, margin 14.7 |
| ethylenediamine | rank 1, margin 1.0 (prior) | rank 1, margin 1.0 |
| isopropylamine | rank 1, margin 2.0 (prior) | rank 1, margin 2.0 |
| triethylamine | rank 1 (only candidate) | rank 1 |

### What limits the ranking

**Single bonds only.** With single bonds only, synthetic top-1 was 0.54. All 46 failures were cases where the best
candidate had exactly the same bond count as the truth for every observed coupling. The generator's truths had
double bonds; the candidates filled the missing valence with unseen atoms instead.

**Bond orders and the hybridization 1J.** Adding them raised top-1 to 0.91. The remaining 9% are again
equal-bond-count ties, which no coupling of this observation can separate. Examples:

- a ring through X versus two separate X atoms;
- which X carries the double bond.

**Real amines.** For ethylenediamine and isopropylamine, the runner-up differs from the truth only in where the
unseen N atoms sit. Both have the same likelihood, so the margin is the prior alone. The J couplings show where the
carbons are and that something is bonded at the free valence. They do not show what that something is.

### Caveats

- **Biased networks.** The real J networks come from fits of the true structure, so they are biased toward it. A
  fair test fits a structure-free J network (group-level couplings without a motif) and then ranks.
- **Synthetic favours route B.** The synthetic test uses the same coupling rules route B was trained on. The real
  networks are the fair comparison, and there the two routes agree.
- **Ring size.** Ring strain and three-membered rings are not penalized beyond the ring cost.

## Next

1. Route B as a GNN on the whole network (all couplings at once, edge classification), and compare it with A.
2. Feed in the coupling uncertainty from the fit (sigma per coupling) instead of the default 0.1 Hz.
3. Add 15N (and 13C-13C, for labelled samples) observations. They resolve the X ambiguity, for example by
   telling which X carries the protons.
4. Run the readout on a structure-free J network from the hypothesis pipeline. This closes the loop: spectrum,
   then J network, then structure, then a structure fit.
