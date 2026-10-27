# DECISION 2026-10-24 22:40 ET — local layout deltas transfer with 6-10x ATTENUATION to the official scorer; price local arms by the attenuated number

## The comparison that makes this measurable
We now have the same layout change measured twice: once on the local H1 SEARCH replay, once on the
official `vcc2026-val-1` panel. `v5b` (~3.5k calls, dense) -> `v8a` (v9 885-call layout) is the
official instance; `layout_native` (~4.5k) -> `incumbent` (v9, 885) is the local instance. Both hold
head/amplitude/totals/B fixed, so both isolate layout.

Anchors ($b$ = no-skill baseline, $r$ = replicate anchor) recovered in `final_push/official_calibration.py`
convert scaled component back to raw: `raw = b + scaled*(r-b)`.

| component | official Δraw (v5b→v8a) | local Δraw (native→v9) | **transfer ratio** |
|---|---:|---:|---:|
| fidelity | **-0.0098** | **-0.0956** | **0.103 (10x attenuation)** |
| reach | +0.0754 | +0.0514 | 1.465 (1.5x amplification) |
| Jaccard | -0.0056 | -0.0385 | 0.157 (6x attenuation) |

## What this means
1. **A local fidelity delta of 0.10 is worth about 0.010 officially** — `0.010/0.29617 = +0.034`
   scaled, `/6 = +0.006` overall. Any fidelity arm whose local gain is < ~0.05 raw is
   **inside the transfer haze and cannot be promoted**; it is not a result.
2. **Reach is the opposite**: local *understates* it by 1.5x. A local reach gain of +0.05 is worth
   ~+0.075 officially. That is why every density arm "looks great" locally and why v8a still lost.
3. **The official fidelity level calibrates reasonably but its SLOPE does not.** Official v9a implies
   raw fid ≈ 0.4931 (`b + (-0.066)(r-b)`); local H1 v9 raw fid is 0.4537. The *level* is within 0.04;
   the *response to density* is 10x smaller. Mechanism: fidelity is `k/max(n_pred, n_conf)`. The
   denominator only punishes under-calling, so density moves fidelity only on targets where
   `n_conf > n_pred`. On the official panel that regime is evidently much rarer than on H1
   (local: 26/79 targets under-called; total `n_conf` 103,308 vs v9 `n_pred` 80,480).

## Consequence for the fidelity lever (Astra's rank-1, and ours)
The lever is real in *direction* but the density route to it is 10x weaker officially than locally.
**Do not launch a density arm as a fidelity play and expect the local number to survive.** A fidelity
candidate must clear ~0.05 local raw to be even visible above the transfer haze, and even then the
official test is the only decision.

## Consequence for the reach lever
Because reach amplifies 1.5x locally, a local arm that gains reach and loses nothing else is
*under*-priced by our own harness — but `v8a` shows the MSE cost that comes with it is not. Pair
every reach arm with a frozen MSE non-regression gate on the official six-component aggregate.

## Standing rule
Never promote a component on its local delta alone. Convert to expected official delta with the
ratios above, then apply the `/6` composite weight, then compare against the arm's MSE/fidelity cost.
