# DECISION 2026-10-24 22:55 ET — local MSE is BLIND to layout (CORRECTED 2026-10-25 00:45: the ratio_of_sums numerator is invariant under layout; NOT the trace cap)

## The measurement
`mse_direction_1024/artifacts/comparison.json`, the paired H1 run (same head/amplitude/totals/B/scorer):

| arm | layout | median calls | raw MSE | raw fidelity | raw reach |
|---|---|---:|---:|---:|---:|
| incumbent | v9 confidence-ordered | 1,037 | **0.9835370889129174** | 0.45333 | 0.11782 |
| layout_native | native dense | 4,475 | **0.9835370889129174** | 0.54916 | 0.06638 |

**The two MSE values are identical to all 16 published digits.** PDS is identical too (0.75365141187926).

## Why (mechanism, from the vendor source)
`metrics/delta.py:690` applies `pred_term = np.minimum(pred_term, PRED_TRACE_CAP_K * real_tn[rows])`
with `PRED_TRACE_CAP_K = 1.0`, frozen in this job's contract. The v9 layout inflates the delete-1
cell jackknife relative to the original integer emission:

| arm | median `jk` | median `original_jk` | **`jk_ratio`** |
|---|---:|---:|---:|
| incumbent (v9) | 4.5946 | 2.2597 | **2.0364** (range 1.789–2.191) |
| native dense | 2.2597 | 2.2597 | **1.0000** (exactly) |

**CORRECTION 2026-10-25 00:45 ET (final) — my first two explanations were both wrong. The correct
mechanism is that the metric is a ratio_of_sums whose numerator sum is invariant under the layout.**

Direct measurement from `layout_cd4_1025/artifacts/per_target.csv`:

| recipient | per-target numerator max abs diff | **sum of numerators** | denominator sum | ratio_of_sums |
|---|---:|---|---:|---:|
| H1 | 1.34e-05 | 0.15641351 vs 0.15641351 (identical) | 0.15903163 | **0.983537 both** |
| CD4 | 8.44e-05 | 0.33867613 vs 0.33867613 (identical) | 0.30258111 | **1.119290 both** |

`expr_mse_unbiased_capped_norm` aggregates as `ratio_of_sums`, and `expr_distance_unbiased` (the
denominator) **reads only the real side, so it is submission-independent by construction**
(`delta.py:1210`, "SUBMISSION-INDEPENDENT ... identical for every submission"). The layout moves
mass between rows — per-target numerators differ by up to 8e-5 — but the weighted sum is preserved,
so the ratio cannot move.

**What is NOT the reason:** the trace cap. Measured directly, `jk_pred` median 4.5946 < `jk_real`
median 5.9781, so `jk_pred > jk_real` on **0 of 79 targets**: the cap never binds on H1. (My earlier
note that the `jk_ratio` 2.036-vs-1.000 difference was erased by the cap is withdrawn.) Also *not* the
reason: a bit-identical `P` — `np.array_equal` did return True, but the *emitted cell matrices* differ
and so do the per-target numerators; the recomputed bulk means agree only to float rounding.

Everything in the "binding consequences" section below still holds; only the explanation changes.

## Why this matters more than any single arm
1. `OFFICIAL_SCORES.md` records v5b→v8a (the same "keep v5b totals, apply the v9 layout" change) as
   **scaled MSE 0.082 -> 0.016, a -0.066 loss**, which alone cost -0.011 overall and is the main
   reason v8a (0.1470) lost to v5b (0.1528).
2. **Locally that loss is exactly 0.000000.** The harness is not mis-scaling MSE — it is structurally
   insensitive to it, because the frozen cap saturates on both arms.
3. Therefore **no layout arm can be ranked locally on the component that decides it.** Combined with
   `20261024_2240_local_to_official_transfer.md` (fidelity attenuates 10x, reach amplifies 1.5x, and
   MSE is blind), the local harness has *no* calibrated channel for a layout decision.

## Binding consequences
- **Stop proposing layout arms as score improvements.** They cannot be evaluated before submission,
  and every visible slot spent learning this is a slot not spent on a measurable lever.
- `G8_jk_dominance` PASSing on both arms is not evidence of safety: it passes *because* the cap binds.
  Report `jk_ratio` alongside every future layout receipt so the inflation is never invisible again.
- Any future arm that changes per-cell dispersion or within-group allocation must state, up front,
  whether it moves `jk_ratio`; if it does, the local MSE number carries no information about it.
- This is the concrete mechanism behind `OFFICIAL_SCORES.md`'s 10-19 note that "the within-group
  layout drives MSE far more than min(B, claim) says" — it is the cap's saturation boundary, and the
  official panel sits on the other side of it from our SEARCH replay.
