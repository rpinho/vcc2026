# DECISION 2026-10-25 00:10 ET — layout effect REPLICATES on a second recipient, but the lever is CLOSED for the final package

## Verdict (frozen rule, applied exactly)
`layout_cd4_1025` COMPLETE, exit 0 both arms, independent full-table + archived-state verification PASS.
**REPLICATES.** All six frozen checks true: fid Δ>+0.020, jac Δ>+0.010, reach Δ<−0.005, three signs match
H1, all three paired-target intervals exclude zero, frozen input hashes unchanged.

| recipient | Δfid | Δjac | Δreach | ΔNMAE | fid 95% interval |
|---|---:|---:|---:|---:|---|
| H1 (79 tgt) | +0.095827 | +0.038508 | −0.051443 | −0.000122 | [+0.0640, +0.1296] |
| CD4 (79 tgt) | +0.058834 | +0.022778 | −0.051629 | +0.000148 | [+0.0354, +0.0854] |

Dense (native) beats sparse (v9) on fidelity and Jaccard, loses reach, in both recipients, with non-zero
intervals. The effect is a genuine property of the layout trade, not H1 noise. **That question is answered.**

## Why it still does not license a package change
The official instance of this same lever already exists: `v5b` (dense) **0.1528** vs `v8a` (v9 sparse)
**0.1470** — dense wins officially too, direction agreeing with local on all three components:

| component | official v5b→v8a (dense→sparse) | local native→v9 (dense→sparse) | transfer |
|---|---|---:|---:|
| fidelity | −0.0098 raw (scaled −0.033) | −0.0956 raw | 10x attenuated |
| Jaccard | −0.0056 raw (scaled −0.015) | −0.0385 raw | 6x attenuated |
| reach | +0.0754 raw (scaled +0.085) | +0.0514 raw | 1.5x amplified |
| MSE | **−0.066 scaled** | **0.000000** (blind) | not measurable |

Applying the measured transfer ratios to the CD4/H1 deltas, dense would buy fidelity+jaccard
≈ (+0.034+0.016)/6 = +0.008, cost reach ≈ −0.084/6 = −0.014, and carry v8a's measured MSE risk −0.011.
**Net ≈ −0.017 overall.** Dense layout is a loss for the final package. Keep the v9 layout.

## The two findings that come out of this and generalize
1. **MSE is structurally invisible locally** (`20261024_2255`). Both arms returned
   `0.9835370889129174` — identical to 16 digits — because `PRED_TRACE_CAP_K=1.0` saturates on both.
   `jk_ratio` is 2.036 for v9 and exactly 1.000 for native; the cap pins both to `1.0 * real_tn`.
   Any future arm that moves `jk_ratio` cannot be judged locally on MSE at all.
2. **The effect is conditional on the hidden reference having many calls.** Deltas are concentrated
   almost entirely in the `n_conf ≥ 1000` stratum:

| stratum | H1 Δfid (n) | CD4 Δfid (n) | H1 Δreach | CD4 Δreach |
|---|---:|---:|---:|---:|
| `<10` | +0.0103 (12) | +0.0264 (23) | **−0.1567** | +0.0311 |
| `10–999` | −0.0003 (41) | +0.0334 (50) | −0.0505 | −0.0925 |
| `≥1000` | **+0.2869 (26)** | **+0.3949 (6)** | −0.0044 | −0.0008 |

The whole gain lives where the reference calls thousands of genes *and* where the reach cost is
negligible. Everywhere else dense either does nothing (10–999) or is actively harmful (the `<10`
stratum costs −0.157 reach on H1 for +0.010 fidelity). A blanket dense policy is therefore a bad
approximation of a good rule — but the good rule needs `n_conf`, which is hidden at scoring time.

## What is left open, honestly
A *conditional* dense policy (dense only where the reference is dense) would be worth roughly
+0.095 H1 / +0.022 CD4 fidelity for ~0.0014 reach on H1 — 98.5% of the fidelity gain for 2.8% of the
reach cost. It is **not implementable today**: `n_conf` is a hidden-scoring quantity, and the
label-free proxies are weak (predicted response magnitude Spearman 0.347 with `n_conf`; native call
count is nearly constant per target — CD4 2920–3256 — so it carries almost no per-target signal;
v9's own call count correlates 0.039, essentially nothing).

Recorded as a real but blocked lever, not a candidate. Do not spend visible slots on it.

## Closed
- Global dense layout for the final package: officially and now doubly locally measured as a loss.
- Layout arms as a fidelity play: transfer attenuates 10x and MSE is unmeasurable.
Keep the safety entry `5NeeoWn96lfm23hgExcm`; no package change.
