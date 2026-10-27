# Amendment N1-CD4 — confirmation rule, registered 2026-10-26 ~14:45 UTC before any CD4 build or score

Authority: decisions/consult_1026_amp round 3 (unanimous on route (b) CD4 with fallback (a)).
Composition rule used for the pass criteria: a criterion is adopted at the strength required by >= 2 of 3 advisors
(a stronger requirement implies the weaker one). Stricter single-advisor criteria are reported, not binding.

When: only after N1 passes its full 79-gate (AMENDMENT_N1.md). If N1 fails on the 79, CD4 is not run.
Transform: the locked N1 transform with alpha 1.25, no retuning: A = mu_CD4 + 1.25(B0_CD4 - mu_CD4), rescaled per
target to ||B0_CD4_i||, self-gene from B0_CD4. B0_CD4 = shipped-model CD4 predictions (mse_1025/inputs/inference_CD4.npz).
mu_CD4 legality: the mean of CD4 truth responses over CD4 targets NOT in the CD4 search fold (same CD4 condition and
response coordinates as the search fold). No other-cell-type substitute. If this cannot be formed legally ->
route (a), and the route taken is recorded.

CD4 PASS (all binding):
1. paired composite point estimate N1 - B0 on CD4 > 0            (Opus; implied by Astra's LB > 0; Fable >= 0)
2. official raw-MSE ratio N1/B0 on CD4 <= 1.03                  (all three)
3. raw PDS(N1) >= raw PDS(B0) - 0.005 on CD4                     (Opus; implied by Astra's >= B0)
Reported, not binding: composite paired 95% LB > 0 (Astra), PDS >= B0 (Astra).
A completed CD4 failure blocks promotion; route (a) may not override it (Astra; consistent with Fable/Opus).

Route (a), only if no legal mu_CD4: full 79-gate + composite >= B0 + 0.008 + raw-MSE ratio <= 1.03 (all three).

Recommendation to Ricardo (only if a route passes): last entry counts (CONSTRAINTS.md), so N1 would REPLACE B0
blind. One page: transform, sha256s, 79 results with paired CIs, raw-MSE ratio, CD4 result or why skipped, expected
gain and worst case, "improvement on the hidden final panel remains unknown; a submission gives no feedback before
the deadline". One N1 submission, by 11-04, B0 file kept as rollback. P1a is a separate decision. Orchestrator never submits.
