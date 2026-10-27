**Progress remains plausible.** The sign audit supports testing selective calls; it does not establish that local precision gains transfer officially.

Run these three arms against paired B0, retaining its bulk model, amplitude and seeds:

1. **`B0_K250`** — reduce the forced-call budget to 250, preserving confidence order.
2. **`B0_K100`** — the same at 100, testing whether precision and Jaccard gains outweigh lost coverage.
3. **`B0_R20_K885`** — retain 885 calls; allocate the strongest p-value positions to 20 genes ranked by source-derived sign certainty × predicted reference detectability. Preserve the remaining order. Use no held-out labels for ranking.

These isolate budget and ordering before introducing adaptive budgets or Skadi selection. Measure **actual significant calls**, including incidental calls. Under the [official definition](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md), fidelity equals precision only when n_pred ≥ n_real; median n_real=46 does not guarantee that for every target.

**Decision instrument:** accept the 53 weak+mid targets as a diagnostic, not a submission gate. Report paired changes separately for weak, mid and strong targets. Reconcile signaudit’s 9/44/26 strata with calib3’s differing counts before weighting.

Use two separate projections:

- Stratum-weighted raw component deltas, converted through official scaling slopes.
- Calib3’s global component-ratio projection, as a sensitivity check.

Do **not** apply global ratios after reweighting: that risks correcting the same mismatch twice. The 31.5/65/3.5 mix is inferred, not observed; sweep plausible mixtures. Recompute MSE at panel level, including B and claims.

**Submission recommendation within ten days: yes, conditionally.** Require an independently held-out recipient/context, paired seed replication, gains distributed across targets, and robustness to stratum/calibration uncertainty. My proposed bar: conservative projected gain ≥0.005, exceeding the documented 0.002–0.003 seed floor, with no material PDS/MSE/NMAE regression. Otherwise retain B0. Ricardo authorizes any submission.