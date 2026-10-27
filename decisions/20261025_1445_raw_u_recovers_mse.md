# 2026-10-25 14:45Z — raw local u recovers the official MSE ordering, scaled MSE does not

## Correction to the previous record

`decisions/20261025_1435_local_mse_blind.md` is right about the clamp but wrong about the
consequence. I wrote that "MSE is invisible to our optimizer". The clamp does erase the
*scaled* signal, but the *raw* u underneath it is retained per arm, and it turns out to carry
the ordering. My first pass at this also had the sign backwards (I compared a
lower-is-better error against a higher-is-better score); the corrected version is below.

## Result

Official MSE on the validation panel vs. raw local `expr_mse_unbiased_capped_norm` (u, lower
is better) on the H1 rich panel, seven recipes:

| recipe | local raw u | local scaled MSE | official MSE |
|---|---:|---:|---:|
| **v5b** | **0.985006** | 0.000000 | **0.082** |
| **v5a** | **0.985618** | 0.000000 | **0.078** |
| v4b | 0.998805 | 0.000000 | 0.000 |
| v3b | 1.045012 | 0.000000 | 0.000 |
| v4a | 1.789320 | 0.000000 | 0.000 |
| v2b | 2.976428 | 0.000000 | 0.000 |
| v1 | 3.054069 | 0.000000 | 0.000 |

**11 of 11 pairs with a non-tied official MSE are ordered correctly by raw local u** — including
the hardest one, v5a vs v5b (local u 0.985618 → 0.985006, official 0.078 → 0.082). The scaled
column orders 0 of them, because every value is identically zero.

Seven of the 21 pairs are officially tied at MSE 0 and are uninformative in either direction.

## Why this is significant

The component that the local harness appeared to be blind to is in fact the component whose
*ordering* it recovers best — better than any of the five live scaled components
(`decisions/20261025_1425_transfer_error_scale.md`: PDS 8/15, NMAE 12/15, fidelity 10/15, reach
9/15, Jaccard 7/15, all excluding v3b).

Caveats, stated plainly:

- **11/11 is 11 pairs, not 11 independent experiments.** The pairs are not independent: they
  share recipes, and eight of the eleven involve v5a or v5b versus an arm with officially tied
  zero. The genuinely informative comparisons are v5a vs v5b (0.004 apart locally, 0.004
  officially — almost suspiciously close) and the four "does this arm earn *any* MSE credit"
  comparisons. The correct reading is: **raw u separates "earns MSE credit" from "does not",
  and it correctly ranks the two arms that both do.** That is a real, useful signal, but it is
  not a calibrated predictor of MSE magnitude.
- The near-exact agreement on v5a vs v5b could be coincidence. Do not treat 0.9856 vs 0.9850
  as implying 0.078 vs 0.082 by proportion.

## What it means for tonight's job

The pass rule I pre-registered for `final_push/dense_bulk_1025` said "scaled MSE > 0". That
clause is void: the anchor makes it unreachable for every arm in the family. The replacement,
recorded *before* any arm has been scored, is:

- report **raw u** for each arm and compare it to the same-run incumbent;
- a *lower* raw u than the incumbent is evidence the arm genuinely improves expression, which
  is the mechanism H4 predicts;
- this is a **diagnostic**, not a promotion criterion. Nothing here overrides the transfer
  finding: the local panel still gets the *overall* ordering wrong 6/15 outside v3b, and PDS
  is anti-correlated with official PDS for exactly the v4b/v5a/v5b family we are working in.

So: a `dense` arm that lowers raw u meaningfully is a mechanistically interesting result worth
Ricardo's attention. It still is not a package.

## Preserved

No anchors changed, nothing re-scored, no package touched. Safety package
`ce78c81c…9dc4` unchanged.
