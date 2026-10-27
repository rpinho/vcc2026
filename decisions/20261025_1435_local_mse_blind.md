# 2026-10-25 14:35Z — the local harness is structurally blind to scaled MSE

## Finding

The frozen H1_rich_search harness derives scaled MSE from a measured baseline/replicate anchor
pair in `rsi/artifacts/scores/H1_rich_search/anchors.json`:

- baseline `expr_mse_unbiased_capped_norm` u = **0.961342**
- replicate u = **0.003359**

Scaling is linear between those two anchors, and the bounded-error class clamps at 0. So any
arm whose u is at or above **0.961342** scores exactly `0.000000`, with no resolution at all.

Every arm we have ever scored locally sits above that line:

| recipe | local raw u | local scaled MSE | **official MSE** |
|---|---:|---:|---:|
| v1 | 3.054069 | 0.000000 | 0.000 |
| v2b | 2.976428 | 0.000000 | 0.000 |
| v3b | 1.045012 | 0.000000 | 0.000 |
| v4a | 1.789320 | 0.000000 | 0.000 |
| v4b | 0.998805 | 0.000000 | 0.000 |
| **v5a** | **0.985618** | **0.000000** | **0.078** |
| **v5b** | **0.985006** | **0.000000** | **0.082** |

v5a and v5b are the two arms that officially *earn* MSE credit — 0.078 and 0.082, i.e.
**8.7% and 8.9% of their official overall**. Locally both are pinned to exactly zero.

All 24 arms in the two Fable-jobs' local artifacts have scaled MSE = 0.000000, range exactly
0.0. This is not because the arms are all equal; it is because the anchor pair on this task is
calibrated to a scale where our whole arm family, incumbent included, is *worse than the
baseline*, and the clamp erases the difference.

## Why this matters

1. **The incumbent is locally worse than the trivial baseline on expression.** u = 0.9835
   vs the anchor baseline's 0.9613. That was never surfaced because the clamp shows 0.
2. **The single largest unclaimed component is invisible to our optimizer.** Official MSE is
   worth ~0.013 of overall for v5b and ~0.012 for billmA (0.061/6). Everything we have tuned
   locally optimized five components with resolution and one component with none.
3. **My own pre-registered rule in tonight's job is unsatisfiable.** I wrote "scaled MSE > 0"
   as a pass condition for `final_push/dense_bulk_1025`. Given this anchor, *no* arm can
   satisfy it, including a hypothetical perfect one that still sits above 0.961342. That
   clause must be recorded as void-by-construction and the arm judged on the 5 live
   components plus raw u, with raw u reported. This is a correction to my own design, made
   before seeing the arms' results.

Both consultants independently warned about this. Astra: "Diagnose the all-zero MSE component
before interpreting it as evidence that expression accuracy cannot improve." Kimi: "a
capped/floored component can remain zero despite genuine raw-MSE improvement… record raw/
uncapped MSE." They were right, and the mechanism is now identified rather than suspected.

## What is actually measurable, and what a fix would require

Raw u *is* retained per arm and does differ (0.9988 to 1.7893 across the historical seven;
0.9711 to 0.9850 across tonight's local arms). So a *relative* expression comparison is
available; the *scaled* credit is not.

To restore sensitivity, the anchor pair would have to be re-measured on a scale where our arms
can fall below the baseline — that is, a baseline anchor built the way the official one is,
from real replicate splits of the same panel, rather than one where our predictions are
uniformly worse. That is a real piece of work, and it must not be done by fiddling constants
until our arms score positive: the anchor is a measurement, not a knob. Any re-derivation needs
its own receipt (what was held out, what was scored, how the split was chosen) and an
independent check that it does not reverse the historical pairs it is supposed to preserve.

This is *not* a path to a better package tonight — a re-anchored local MSE still would not
transfer (see `decisions/20261025_1425_transfer_error_scale.md`: the H1 panel inverts 6/15
pairs outside v3b, and local PDS is anti-correlated with official PDS for the v4b/v5a/v5b
family). It is a diagnosis of why the harness has been pointing us at the wrong levers.

## Action taken

- No change to the frozen anchors. Nothing is re-scored or re-anchored tonight.
- `final_push/dense_bulk_1025` continues; its MSE clause will be reported as void rather than
  failed, and its verdict will rest on the five live components and raw u.
- Recorded so the same "MSE is exactly zero" observation is not again read as "expression
  cannot improve".
