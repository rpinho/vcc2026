# DECISION 2026-10-25 00:45 ET — MSE is a function of rho alone, and rho is the wall

## The algebra
`expr_mse_unbiased_capped_norm` is `ratio_of_sums(expr_mse_unbiased_capped / expr_distance_unbiased)`
(`catalog.py:345`), where per target

```
numerator_p   = (||mu_pred - mu_real||^2_G - min(jk_pred, 1.0*jk_real) - jk_real) / G
denominator_p = (||mu_real - mu_ctrl||^2_G - jk_real,p - jk_ctrl/N) / G
```

Both sides are *unbiased* squared distances, so the sampling corrections cancel in expectation and
the score reduces to the **relative direction error**. Writing `a = ||pred||`, `b = ||real||`,
`rho = cos(pred, real)`, the numerator is `a² + b² - 2ab·rho` and the denominator is `b²`. At the
optimal scalar amplitude the numerator minimises at `b²(1 - rho²)`:

```
u_min = 1 - rho^2      =>      scaled_MSE = (0.989634 - (1 - rho^2)) / 0.953025
```

## The conversion table (this is the useful artifact)
| rho | u_min | scaled MSE |
|---:|---:|---:|
| 0.15 | 0.9775 | **0.0127** |
| 0.262 | 0.9315 | 0.061 (what we emit) |
| 0.393 | 0.8457 | 0.151 (Sidechain) |
| 0.618 | 0.6180 | 0.390 (best entry inside our PDS band) |
| **0.707** | 0.5007 | **0.513 (top-10 median)** |

**Correction (2026-10-25 00:40) — our local measurement is not what I first wrote here.** The local
run at `bulk_target_sum=50000` gives `jk_pred` median 4.5946 < `jk_real` median 5.9781, so the cap
**never binds** on H1. Local mean per-target rho on the SEARCH replay is **0.1535**, which under the
formula implies scaled MSE ≈ **0.013–0.016**, yet the harness emits 0.061. That gap means the
amplitude is not exactly optimal and/or the aggregate weighting differs from a per-target mean
(`ratio_of_sums`, not mean-of-ratios). Treat 0.1535 as an upper bound on usable signal and 0.061 as
what our current (slightly non-optimal) amplitude actually converts it into.

## What it means
1. **MSE headroom is rho headroom.** Every point of scaled MSE above us corresponds to a direction
   correlation we do not have. Top-10 median 0.513 needs rho ≈ 0.71. We have 0.15.
2. **No amplitude, layout, masking or allocation knob can cross this.** Those knobs move `a` (and the
   small correction terms); the ceiling they approach is `1 - rho²`. Amplitude is already within a
   few percent of optimal, which is consistent with the closed amplitude sweep.
3. **The 0.061 → 0.39 headroom seen inside our own PDS band (260 entries) is real but is bought with
   rho**, i.e. with response-direction information, not with emitter tuning.
4. **This is consistent with, and explains, the cross-line transfer ceiling**: rho_opt ~0.12–0.14
   cross-line, 0.1435 same-line, 0.732 same-experiment. The final panel is cross-line, so
   rho ≈ 0.13–0.15 is the honest expectation and MSE ≈ 0.01–0.06 is the honest ceiling for our
   current *information*. The reason we are at 0.061 and not 0.015 is that our amplitude conversion
   is imperfect, not that we are leaving a large MSE lever unused.

## Binding consequences
- **Do not launch further MSE-tuning arms** (amplitude, bulk multiplier, deflation, dispersion,
  learned-norm): each is bounded above by the same `1 - rho²` and the empirical sweeps all returned
  f=1/optimal already. `mse_shared_5e4_1025` confirms this at the correct coordinate.
- **The only lever that moves MSE materially is acquiring/reconstructing response-direction signal**
  that raises rho on the 300 final targets. That is a data/biology problem, and the measured same-line
  premium (1.19x) makes it a poor bet within 10 days.
- When any future arm claims an MSE gain, require the implied rho change. If it does not raise rho, it
  cannot be real.
