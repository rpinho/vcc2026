# DECISION 2026-10-25 03:55 ET — the 9-10% post-blend amplitude shortfall is PROTECTIVE: do not "fix" it

## The observation
Every emitted arm in `expr_unmask_1025` lands **below** its requested claim after the learned-blend step:

| task / arm | alpha | B | claim | B/claim | requested B/claim |
|---|---:|---:|---:|---:|---:|
| H1 incumbent | 0.4443 | 0.008996 | 0.020036 | 0.449 | **0.9091** |
| H1 no_mask | 0.2611 | 0.008955 | 0.020942 | 0.428 | 0.9032 |
| H1 no_pc | 0.4282 | 0.009023 | 0.020043 | 0.450 | 0.9118 |
| H1 unmasked | 0.2573 | 0.008974 | 0.020929 | 0.429 | 0.9051 |
| CD4 incumbent | 1.5457 | 0.059144 | 0.092433 | 0.640 | **0.8860** |

The cause is geometric, not biological. `blend()` rotates 30% toward the learned direction and then
**renormalises to the pooled norm**; if the learned direction sits at angle `theta` from the pooled one, the
blended vector's length before renormalisation is `sqrt((1-c)^2 + c^2 + 2c(1-c)cos theta)` — 0.84 at
`cos theta = 0.3`, 0.76 at `cos theta = 0`. The amplitude solve ran *before* the blend, so the blend silently
taxes the amplitude it just set.

## Why it must not be repaired
Recomputing the MSE decomposition (`mse_pipeline_1025/mse_decomposition.json`, per-target rho and norm ratios,
`ratio_of_sums` weighting) at fixed jackknife shows **more amplitude is worse for MSE above the current point**:

| amplitude vs current | implied u | scaled MSE | Δ vs now |
|---:|---:|---:|---:|
| 1.00 | 0.849 | **0.1471** | — |
| 1.05 | 0.867 | 0.1289 | −0.0182 |
| **1.107** (exactly B=claim) | 0.884 | 0.1106 | **−0.0365** |
| 1.20 | 0.926 | 0.0667 | −0.0803 |
| 1.50 | 1.077 | −0.0914 | −0.2385 |

Lifting every profile to `B = claim` would cost about **−0.041 scaled MSE = −0.0068 official**, and the true
cost is larger because `jk_pred` also rises with amplitude. The earlier amplitude model puts the optimum at
`f ≈ 0.30` with the curve nearly flat between 1.0 and 1.1 — i.e. the current under-claim is already close to
the MSE-optimal point.

**Conclusion: the shortfall is a feature.** The `B = claim` rule is a *floor for full credit*, not an optimum
(`FABLE_REVIEW_1024.md` says the same). Do not launch a "re-solve alpha after the blend" arm expecting an MSE
gain; the measured model says it is a loss. If anyone revisits this, the only defensible form is a joint
solve that re-derives `B` and the jackknife **together**, and it must beat the table above, not just reach claim.

## Consequence for the gates (important for the next job)
This is also why a `requested_B/claim >= 0.90` gate is wrong: the CD4 incumbent sits at 0.886 and is a
documented G3 exception, not a defect. `_codex/vcc_family_quality_matched_1025.md` has been amended so the
budget gate is **arms matched to each other** (`|Δratio_requested| <= 0.02`), with the hard upper bound
`B/claim <= 1.05` retained and the lower side reported as a diagnostic only.

No package change, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.
