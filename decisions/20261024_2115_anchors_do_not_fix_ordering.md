# 2026-10-24 21:15 ET — recovered anchors fix SCALING, but do NOT fix historical mis-ordering

Correction to the 20:45 "breakthrough" framing. The anchors are real and well-fit (holdout median
error 1e-4…1e-3, structural checks pass), but they do **not** rescue the harness's rank calibration.

## Test
Re-scored every retained `final_push/artifacts/de_check/H1_rich_search/*` arm with the recovered
anchors. The audit's known reversal is v5b-native vs v9-layout-at-885 (officially billmA 0.16438 >
v5b 0.1528):

| arm | local PDS raw | local MSE norm | official-equivalent overall |
|---|---:|---:|---:|
| incumbent (v5b native emission) | 0.7537 | 0.9835 | 0.1471 |
| levers_base_885 (v9 layout, 885) | 0.7235 | 0.9811 | 0.0776 |

**Still reversed.** Officially v9 layout *raises* PDS (v5b 0.680 → billmA 0.706, and v8a's
layout-only-on-v5b-totals kept it at exactly 0.680); locally the v9 arm *lowers* PDS by 0.030.
So the reversal lives in the RAW metric, not in the anchor pair.

## What this changes
1. Do **not** use the local harness to rank emitted-layout variants. That is the artifact the
   `harness_audit_20261023` report was describing, and the anchors do not remove it. Local PDS also
   sits ~0.05 above official (0.7537 vs 0.680 for the same recipe), consistent with the local replay's
   deeper truth (400 real cells/target) making the ranking problem easier.
2. **Do** trust direct scorer-space direction measurements (`rho_bar`, `rho_optimal`, `u = 1 - rho_opt^2`).
   Those do not depend on emitted layout at all, and the model reproduces the harness's capped norm to
   5 decimals. That is why the amplitude-closure result stands.
3. Net: the anchors are still worth having (they give correct *absolute* scale, which is what made the
   MSE deficit legible and the amplitude-closure provable), but the "we can now rank candidates locally"
   claim in the 20:45 entry is **withdrawn**. A local win is necessary evidence, never sufficient.

## Consequence for promotion
The gate `historical_order = false` on every arm stays correct and must not be waived. No package.
