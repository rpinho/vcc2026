# Amendment RENDER_D — 2x layout call budget (1770) on B0 (registered before scoring)

Authority: decisions/consult_1026_amp round 6 (2026-10-26 18:12 UTC), unanimous.
- U885 accepted as valid under the corrected check below. This is a post-result correction of a mistaken
  invariant (Astra): the original clause "PDS and raw MSE bit-identical" (AMENDMENT_RENDER_U.md) stays in the audit
  trail unchanged. The preregistered branch d < -0.002 fired: confidence allocation contributes at 885 calls
  (Opus: real effect, but CI upper end -0.00056 sits inside +/-0.002, so its size is uncertain).
- One diagnostic arm, then no further arms (Opus, Astra: no budget sweep). Nothing is submitted without Ricardo.

Arm D1770: B0 response verbatim (results/predictions_D1770.npz), B0 confidence ranking verbatim, call budget 1770
instead of base['renderer']['budget'] = 885. Seed base 20260910, anchors and depth identical. Targets with fewer
eligible genes get all of them (selection takes order[:nconf]). Runner score_render_d.py = score_run.py
confidence_all/run with textual substitutions only: confidence cache path; the layout budget argument (and the
render-json budget record) -> BUDGET.get(arm, base budget).

Validity (corrected clause; per-item median of R6): the rendered pseudobulk (count_sums) and PDS are bit-identical
to B0; the scaled-MSE contribution to the composite is unchanged (0 for both); raw expr_mse_unbiased_capped_norm is
reported and no more than 0.1% worse than B0 (Opus). Raw and scaled MSE deltas are reported separately (Astra).
Failing any of these voids the run.

Stop rule (Opus R5, binding, all three R6): paired composite d = D1770 - B0 (20k boot, seed 20261026).
Point estimate < +0.002 OR 95% lower bound <= 0 -> CLOSE renderer research; B0 is final.
Otherwise PASS (diagnostic only) -> back to the advisors for Astra's untouched-target confirmation design and
Opus's R5 Q2 list. No promotion from H1 SEARCH.
