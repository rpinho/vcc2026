# DECISION 2026-10-25 04:25 ET — CD4 completes; the expression leg is closed

## Result
Official-anchored deltas from `price_official.py CD4_depth192_search` (positive = better):
- `no_mask`: PDS +0.0144, MSE +0.0024, NMAE +0.0015, fidelity +0.0180, reach −0.0234, Jaccard −0.0009;
  **sum +0.0121, Δofficial overall +0.0020**.
- `no_pc`: PDS −0.0040, MSE −0.0004, NMAE −0.0017, fidelity −0.0009, reach −0.0357, Jaccard −0.0002;
  **sum −0.0428, Δofficial overall −0.0071**.
- `unmasked`: PDS −0.0090, MSE +0.0001, NMAE +0.0053, fidelity +0.0159, reach −0.0529, Jaccard −0.0011;
  **sum −0.0417, Δofficial overall −0.0069**.

These reproduce the values in the 00:10 decision summary exactly. The job's own `score.json` `scaled` column uses
**local** anchors that clip scaled MSE to 0 at this operating point, so it is not a contradiction — it is a
different anchor set. Never mix the two conventions in one table.

## Verdict
- All four arms pass every required gate on both tasks; independent verification and archive replay PASS.
- **Nothing is promoted.** Every CD4 MSE bootstrap interval includes zero, and raw reach falls in every arm
  (`no_mask` −0.0207, `no_pc` −0.0316, `unmasked` −0.0469), exceeding the pre-registered 0.02 loss allowance.
- H1's `no_pc` gain (+0.0016 official overall) does **not** replicate in CD4 (−0.0071). Under the
  both-recipient rule, `no_pc` is closed as a package candidate.
- **The expression-leg programme is CLOSED** (decision `20261024_2225` lineage). No further mask/PC arms.

## Mechanism, retained for future reading
Removing the mask collapsed `alpha` from 0.444 to 0.261 on H1 and from 1.546 to 0.795 on CD4, at constant
B/claim. The `B=claim` solve converts a shape change into an amplitude change, so any "more signal" arm
measures the solver unless it reports `alpha` and B/claim alongside the six components.

## What this does NOT close
The `family_quality` **source-direction** lever — a different mechanism, because it replaces the direction
*before* calibration instead of editing the incumbent's mask. That job launched 04:21 ET as
`vcc_family_quality_matched_1025` with three arms (`incumbent`, `family_quality`, `family_equal`) under the
gate in `decisions/20261025_0410`.

No submission, no package change. Safety entry `5NeeoWn96lfm23hgExcm` and package
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4` unchanged.
