# Lactic acid: slow O-H exchange test (2026-10-11)

- Data: raw/2026-03-16-l-lactic-acid/L_+_lactic_acid/data (series runs/series/batch_lactic-acid).
- Question (Xuehan): the fast-exchange fit misses the CH3 J-band line at 127.33 Hz and places the CH line at
  149.58 Hz off; is the hydroxyl proton coupled (slow exchange)?
- Prior check: grid scan of 2J(C2,H1) x 3J(H1,H2) with the OH decoupled: the current values are the only deep
  minimum, so no coupling combination of the fast-exchange model produces the two lines.
- Command: `python scripts/run_batch.py configs/batch_2026-10-10_v8.json` (exchange "slow", OH proton kept);
  results table docs/analysis/2026-10-10_batch_results_v8.md; fits runs/processed/batch_lactic-acid_v8_s1, _s2.

## Results

| model | stage 1 | stage 2 | 1J(C1,H) | 2J(C1,H2) | 3J(H1,H2) | J to OH |
|---|---|---|---|---|---|---|
| fast exchange (batch_lactic-acid_s2) | 0.346 | **0.223** | 147.07 | -4.26 | 5.98 | - |
| slow exchange (batch_lactic-acid_v8_s2) | 0.382 | 0.263 | 146.49 | -6.94 | 7.69 | 2J(C1,HO) -6.9 +- 0.8 |

Slow exchange is worse in both stages (stage 2 took 35 min against 8 min). The OH coupling is poorly determined
(std 0.8 Hz) and pulls 2J(C1,H2) away from its prior.

## Conclusion

A coupled hydroxyl proton does not explain the 127.33 / 149.58 Hz lines; the fast-exchange fit stays the result
(conditional numerical result of the given structure). Open: sample composition (lactide / oligomers of lactic
acid form in concentrated solution and give CH3 / CH bands close to those of the monomer).

Commit: (this commit)
