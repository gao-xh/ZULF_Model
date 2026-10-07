# J to structure, route B: learning on J versus on the reduced coupling K (2026-10-07)

- Data: synthetic generator pairs (configs/couplings_v1.json, 13C / 1H); the four real amine J networks
  (configs/j_networks). No experimental spectra.
- Question: does route B learn better from K (D51) than from J? (Xuehan: "K和J学习对比".)
- Code: 444f93c plus the K input mode of zulf_model.structure.edge_model (EdgeModel.mode, D52) and
  scripts/j_k_learning_comparison.py.

## Setup

K mode feeds the model the equivalent 13C-1H coupling J (g_13C g_1H) / (g_A g_B), i.e. K up to a constant, so
13C-1H rows have identical features in both modes and H-H rows are scaled by g_13C / g_1H = 0.2515. Same MLP,
same rows, same seeds for both modes; only the input changes.

    python scripts/j_k_learning_comparison.py --seeds 1,2,3 --out runs/processed/j_vs_k

4000 training graphs per seed (about 36000 coupling rows), 20 epochs; held-out rows from 500 graphs
(seed + 1000; 5414 rows for seed 1); 100 synthetic ranking cases (seed 7, 2-5 heavy atoms); the deuterated test
replaces every 1H by 2H in the held-out rows (J(C,D) = J(C,H) x 0.1535, J(D,D) = J(H,H) x 0.0236), untrained.

## Results

Mean +- sample standard deviation over seeds 1, 2, 3 (runs/processed/j_vs_k/summary.json):

| Input | Bond-count accuracy, 1H | Bond-count accuracy, 2H (transfer) | Synthetic top-1 | Synthetic top-3 | Real networks (rank of truth) |
|---|---|---|---|---|---|
| J | 0.844 +- 0.007 | 0.298 +- 0.013 | 0.91 | 1.00 | 1, 1, 1, 1 (all seeds) |
| K | 0.843 +- 0.006 | 0.843 +- 0.006 | 0.91 | 1.00 | 1, 1, 1, 1 (all seeds) |

Baseline: always predicting the most frequent bond count (3) scores 0.322 on the held-out rows (CH rows 0.349,
HH rows 0.500).

## Conclusion

- On data of one isotope J and K learn equally well (difference 0.001, within the seed spread); ranking
  synthetic and real structures is identical. This is expected: the model already knows the element pair, and
  K only rescales each pair by a constant.
- K carries over to other isotopes without retraining: on deuterated couplings the K model keeps its accuracy
  (0.843), the J model falls to 0.298, below the majority-class baseline (0.322). For labelled samples (2H, 15N)
  or mixed isotopes the K input is required; the primary isotope effect on K (about 1 %) is neglected in this
  test, so the real transfer will be slightly worse than 0.843.
- Decision proposed (to discuss with Xuehan): keep J as the default input while all observations are 13C / 1H;
  switch route B to K before any labelled-isotope or N-coupling observations are used.

## Open points

- N couplings (15N, 14N) are not in the generator observations; a K-based model with an N-H / C-N kind needs them.
- The primary isotope effect on K is not modelled.
- Route A (rule likelihood) is unchanged (J ranges per element pair).

Commit: (this file's commit)
