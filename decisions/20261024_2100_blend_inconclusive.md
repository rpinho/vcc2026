# 2026-10-24 21:00 ET — blend-coefficient result is INCONCLUSIVE (my method, not the model)

## What I ran
`ceiling_normmatched.py` and `blend_sweep.py`, both trying to compare the pooled-source direction, the
learned MLP direction, and the shipped 0.3 blend with amplitude held fixed.

## The problem
They disagree on the same quantity:

| run | pooled-direction rho_opt |
|---|---:|
| `ceiling_normmatched.py` | 0.1410 |
| `blend_sweep.py` (coef=0.0) | 0.1004 |

Cause: the two scripts match amplitude to **different** targets. `ceiling_normmatched` matches to the
budget of the actually-emitted counts (`B = 0.1428`); `blend_sweep` matches to the budget of
`ctrl + ensemble_response` (`B` computed from a direction, not from emitted counts). Cosine here is
amplitude-sensitive because every candidate is mapped to counts as `ctrl + f·direction` before the
log transform, so a different `f` moves the cosine materially. My 400-point linear grid over
f ∈ [0.02, 20] also cannot resolve the optimum: it picked f = 0.120 and f = 0.170, both at the extreme
low end, meaning the true optimum is below the grid's useful resolution.

## What this does and does not license
- **Does not** license "the blend is harmful". The measured gap between coef 0.0 and 0.3 is
  rho_opt 0.1004 vs 0.0974 (+0.003) — the same order as the discrepancy between my own two scripts.
- **Does not** license "raw sources beat the shipped direction".
- **Does** establish that this specific comparison needs the emitter's own `solve_alpha` (which bisects
  to an exact budget) instead of an ad-hoc scale grid, before any conclusion is drawn.

## decision summary
Recorded as a precise open question for `vcc_mse_direction_1024`: *does blend coefficient 0 beat 0.3
when both are amplitude-solved by the emitter itself, on the same fold, with rho_opt as the statistic?*
Method requirement: use `v5.solve_alpha`, not a scale grid; report the achieved B for each arm so the
match is auditable. If the answer is within ±0.01 rho_opt, close the question as null.

This is the third thing tonight I have had to correct by checking rather than trusting. The pattern is
consistent: derived quantities (cosine at unverified amplitude, metric ratios over tiny denominators,
u = 1−2rho) have all misled; only the exactly-reproduced identities held. Keep the standard.
