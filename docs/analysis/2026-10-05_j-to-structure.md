# 2026-10-05 J to structure: direct readout (route A) and learned likelihood (route B)

- Samples: isopropylamine, ethylenediamine, triethylamine, N-ethylmethylamine (J networks from the fast-exchange
  fits, configs/j_networks/*.json); synthetic generator pairs.
- Question: Ashok, can the structure be read from the J couplings alone? The user asked for both routes, keeping
  the direct readout.
- Code: zulf_hypothesis/j_structure.py, zulf_model/structure/{observations,edge_model}.py,
  scripts/j_structure.py, scripts/j_structure_benchmark.py; design in docs/J_TO_STRUCTURE.md.

## Commands

```
python scripts/j_structure.py configs/j_networks/*.json --top 3
python scripts/j_structure_benchmark.py --samples 100 --model runs/models/j_edges_v1.json \
    --out runs/processed/j_structure_benchmark/report.json
python scripts/j_structure_benchmark.py --samples 40 --heavy-atoms 4,6 --max-atoms 6 \
    --model runs/models/j_edges_v1.json --out runs/processed/j_structure_benchmark/report_4-6.json
```

## Results

| Version | Synthetic top-1 / top-3 (2-5 heavy atoms, 100) | Amines rank 1 |
|---|---|---|
| A, single bonds only | 0.54 / 0.79 | 4 / 4 |
| B, single bonds only | 0.54 / 0.79 | 4 / 4 |
| A, bond orders + hybridization 1J | 0.91 / 1.00 | 4 / 4 |
| B, learned bond count + hybridization | 0.91 / 1.00 | 4 / 4 |

Route B was trained on 4000 graphs (48436 couplings) in 22 s. Its held-out bond-count accuracy is 0.85.

### Failure analysis

- **Single bonds only.** All 46 misses were equal-bond-count ties: the best candidate predicts the same bond count
  as the truth for every observed coupling. The truths had double bonds, and the candidates filled the free
  valence with unseen atoms instead.
- **With bond orders.** The 9 remaining misses are ties of the same kind: a ring through X versus two X atoms,
  and which X carries the double bond.
- **Real amines.** The runner-up for ethylenediamine (margin 1.0) and isopropylamine (margin 2.0) differs from the
  truth only in the X placement, and both have equal log L. The J values do not identify the unseen atom.

## Conclusion

On these networks, the J couplings determine the carbon skeleton and the bond orders through 1J. They do not
determine the unseen heteroatoms. Routes A and B agree on all four amines.

Caveat: the amine networks come from fits of the true structure. A structure-free fit is the next test.
