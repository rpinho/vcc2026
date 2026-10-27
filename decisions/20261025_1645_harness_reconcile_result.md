# 2026-10-25 16:45Z — reconciliation finished: my "2x gap" was wrong; the rho-identity needs correction

The bounded diagnostic that all three advisors authorized ran and **stopped in Phase 1**, because its
own stop rule fired. It refuted the claim it was sent to explain. Record:
`final_push/harness_reconcile_1025/REPORT.md`; verification PASS in `verification.json`;
package SHA unchanged (`ce78c81c…9dc4`, 2,595,584,000 B) before and after.

## What was refuted (all mine)

1. **The 2x gap was an anchor-provenance error.** `0.961342 / 0.003359` are the *local H1* endpoints,
   not official billmA anchors. Applying them to an official score is invalid. Substituting the
   retained approximate official fit (`0.989634 / 0.036609`) moves the back-solved u by 0.028594,
   which is **35.46%** of the claimed gap on its own — enough to trigger the stop rule. Remaining
   cross-panel difference: unresolved and not testable with retained artifacts.

2. **`u_min = 1 - rho^2` is not an identity for this scorer.** The exact decomposition is
   `u = (1 + a^2 - 2*a*rho_global - e - q) / (1 - e - c)`, where the prediction-budget term
   `q = K/T` and reference-sampling terms `e`, `c` survive. The across-target budget binds hard here:
   **B/claim = 0.448973** (the per-target `jk_pred <= jk_real` cap binds on **0/79** targets, so the
   cap was never the binding constraint — the budget was). rho alone does not determine u.
   The conditional rho needed for 0.061 is **0.186-0.214**, not 0.31.

3. **`0.1535` was three different estimands conflated.** Mean per-target cosine 0.153497, globally
   pooled cosine **0.121114**, amplitude-oracle projection RMS 0.152113. Only the middle one enters a
   global aggregate.

4. **The local replay itself is arithmetically sound.** Independent reconstruction from raw frozen
   H5AD cells reproduces u = 0.983537088913 to 4.77e-15, with reference means and jackknifes exact.

## What this corrects in the record

`decisions/20261025_0045_mse_ceiling_is_rho.md` — the identity is wrong as stated. Its *conclusion*
(amplitude/layout knobs cannot cross the ceiling) is still supported empirically by the closed sweeps
and by item 2 above, but the clean `1 - rho^2` argument must not be cited again.

`decisions/20261025_0115_mse_headroom_is_unlocked.md` — the qualitative finding survives and is now
better supported: 289 entries at our PDS carry far more direction signal, so the gap is a *modelling*
gap, not a data gap. But its `rho_implied = sqrt(1 - u)` column used the refuted inversion, so those
rho values are illustrative only. Its concrete proposal (stop masking the expression leg, keep the
call set) is **not** refuted and remains the single highest-value untested mechanism.

## Why this does not reopen the submission

Fable's rule is unchanged and decisive: the server is closed, so nothing produced now can be scored or
selected. Astra and Kimi's caps were 4-16 h of learning-only work; this consumed ~1 h and produced a
result. The verdict stays **preserve `5NeeoWn96lfm23hgExcm`**.

## What the next cycle should build (from the report, priority order)

1. Preserve per-context raw numerator/denominator rows, bulk moments, both JK terms, both correction
   bounds and exact anchor identifiers with official results. Report mean cosine, global cosine and RMS
   projection under separate names. Never back-solve official correlation from one aggregate loss.
2. Preserve exact recipe/config/emitter identities and reference-cell selections with native emitted
   sufficient statistics, so controlled pairs exist.
3. Build representative native SEARCH folds (target coverage, truth depth, control construction
   matched), with a predeclared ordering test. Do not fit weights to known official ranks.

No follow-up job launched. `learnednorm` remains void. No submission.
