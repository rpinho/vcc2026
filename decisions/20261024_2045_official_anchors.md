# Decision 2026-10-24 20:45 ET — official anchors RECOVERED from the public leaderboard

## What was blocking
`harness_audit_20261023/REPORT.md` and both overnight jobs concluded: "no trustworthy local ranking
check is established", "no calibrated local-to-official noise bound", "no D/E/F setting was fitted",
"the retained H1 incumbent replay is bitwise identical; numerical repeatability is not transfer
uncertainty". Every promotion decision was blocked on this. Both jobs burned the night on sweeps that
could not be promoted even if they had won.

## What we found
`final_push/fable_tests_1024/leaderboard.json` (fetched 10-24 11:00, 1635 val entries) contains, for
every published entry, the **raw value and the scaled score of all six metrics**
(`pds_cosine`+`score_pds`, `expr_mse_unbiased_capped_norm`+`score_mse`, …).
Since `scaled = (baseline − raw)/(baseline − replicate)`, each entry with an interior score is one
equation in the two unknown anchors. 626–1440 usable entries per metric pin them down.

Recovered official val anchors (robust IRLS, median |residual| 4e-5 … 1e-3):

| metric | baseline b | replicate anchor | interior n |
|---|---:|---:|---:|
| pds_cosine | 0.50070 | 0.95112 | 1440 |
| expr_mse_unbiased_capped_norm | 0.98963 | 0.03661 | 626 |
| de_wilcoxon_lfc_nmae | 1.00041 | 0.39456 | 1113 |
| de_wilcoxon_direction_fidelity_yield_raw | 0.51261 | 0.80878 | 302 |
| de_wilcoxon_direction_reach_raw | 0.07893 | 0.96555 | 1421 |
| de_wilcoxon_sig_jaccard | 0.03032 | 0.40432 | 369 |

Consistency checks passed: the MSE fit reproduces the `s==0` clipping constraint (b=0.9896 <= 0.9935 =
min raw among score-zero entries); PDS baseline ≈ 0.5 and reach baseline ≈ 0.079 match the metric's
documented no-skill construction.

## Why this matters
Our local harness anchors are **b=0.95870, replicate=0.00488** versus official **b=0.98963, rep=0.03661**.
The local harness treats the replicate anchor as nearly perfect (0.0049) when officially it is 0.0366,
and local "predict the control" maps to 0.9587 instead of 0.9896. That is exactly the mis-ordering
`harness_audit_20261023` could not explain — the local anchor pair, not the scorer code, was wrong.

## The real deficit, now measured on the official scale
| metric | ours | top-10 median | rank 10–50 median |
|---|---:|---:|---:|
| pds_cosine | 0.706 | 0.831 | 0.813 |
| **MSE** | **0.061** | **0.513** | **0.445** |
| NMAE | 0.116 | 0.282 | 0.244 |
| fidelity | −0.066 | 0.080 | 0.045 |
| reach | 0.186 | 0.302 | 0.272 |
| jaccard | −0.016 | 0.020 | 0.012 |
Overall 0.16438 = rank 451/1635. **MSE alone is 0.384 of the 1.0-scale deficit — 45% of the entire
gap**, and we are *below the baseline* on fidelity and jaccard.
