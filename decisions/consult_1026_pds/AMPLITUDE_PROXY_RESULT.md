# The amplitude experiment, previewed at zero cost — and the advisors' premise is half wrong again

Computed on the saved h2h2 matrices (79 targets, response space, self-gene excluded). Proxy only; the
official renderer applies token/layout/count effects this cannot model. But it answers the structural
question the advisors asked for Monday.

## Result 1: Astra was right that PDS must be MEASURED, and it turns out it IS invariant
Pure positive rescaling of B0 leaves the retrieval rank **exactly** unchanged:

| alpha on `alpha*B0` | mean rank of own truth among 79 |
|---|---:|
| 0.5 | 16.747 |
| 1.0 | 16.747 |
| 1.5 | 16.747 |
| 2.0 | 16.747 |

Identical to 3 decimals, exactly as the cosine algebra predicts. **This confirms the mechanism is
sound**: scaling cannot hurt PDS through the retrieval geometry. (It could still move DE significance
via the renderer — that remains unmeasured and is the one thing worth an official run.)

## Result 2: the `mu + alpha*(B0-mu)` form does NOT behave like scaling
The advisors' formula is not pure scaling — it interpolates toward the training mean. Measured:

| alpha | rank | est. PDS | raw MSE |
|---|---:|---:|---:|
| 0.000 (mu only) | 40.00 | 0.506 | 0.00321 |
| 0.343 (mse_1025's scalar) | 31.19 | 0.636 | 0.00331 |
| 0.500 | 27.20 | 0.695 | 0.00348 |
| 0.750 | 20.99 | 0.787 | 0.00394 |
| **1.000 (= B0)** | **16.75** | **0.850** | **0.00461** |
| 1.250 | 15.48 | 0.869 | 0.00549 |
| **1.500** | **15.32** | **0.871** | **0.00658** |
| 2.000 | 15.57 | 0.867 | 0.00941 |

Two things fall out:

1. **alpha < 1 is catastrophic for PDS.** At the CD4-derived scalar 0.343, estimated PDS is 0.636 —
   that alone is **−0.0745 composite points**, versus an MSE gain of only +0.0050. Shrinking is a
   terrible trade, and this independently explains `mse_1025`'s fidelity drop from 0.1767 to 0.1251.
   **The advisors' "fit alpha by LOO" instruction would very likely land on alpha < 1 and destroy the
   score.** That is the single most actionable thing in this document.

2. **alpha slightly > 1 looks mildly positive here (1.25-1.5), but the MSE cost is severe** (0.00461 →
   0.00658 at 1.5, a 43% worsening) and the whole apparent gain comes from the *oracle* `mu`.

## The fatal caveat — and it invalidates the apparent win
`mu` in this sweep is the **mean of the 79 held-out truths**, which is not available at prediction time.
The real Monday arm must use the 168-target training mean, which is a much noisier and smaller-sample
estimate (and B4, the actual training-mean arm, scored rank 40.0 — as bad as the oracle-mu arm here,
so this proxy likely overstates what a legitimate `mu` can do).

**So the honest reading is:**
- Pure scaling: PDS-safe (proven), but MSE-neutral — it cannot buy points, because B0's amplitude is
  already near its own MSE optimum at alpha≈0.10 in *response* space yet the official metric's own
  optimum is elsewhere. Scaling alone is a no-op or a loss.
- Mean-interpolation: PDS-destroying below alpha=1, and the apparent gain above alpha=1 is an oracle
  artifact that will not survive using the real training mean.

## What this means for Monday
The advisors converged on amplitude calibration as "the one mechanism B1 did not isolate." The proxy
says it is **probably a dead end for the reasons they each half-anticipated**: Fable thought alpha>1
might fix under-calling (it mildly helps rank but costs MSE), Astra and Kimi both warned PDS
preservation must be measured (it is preserved for pure scaling and destroyed for mean-interpolation).

**Recommendation to put back to them:** the only version worth an official run is **pure scaling**
(`alpha*B0`, no mean term), because it is the only variant that is PDS-invariant by construction and
therefore the only one that cannot spend our crown jewel. Test alpha in {0.9, 1.1} — a narrow band
around 1, since MSE is steeply worse outside it — and measure the full ledger once. Expected value is
low but the cost is ~30 min and one arm.

Everything above is a proxy on 79 targets with an oracle mean where noted. It does not substitute for
the official run, and it does not authorize any package.

---

## Result 3: the proxy passes a hard sanity test, and it says NO arm beats B0

B4 (training-mean) is by construction **the same vector for all 79 targets**. A constant prediction has
identical cosine to every truth, so its PDS must sit at chance — i.e. a rank near the middle of 79
(≈40) and PDS ≈ 0.50. The proxy returns **rank 40.01, PDS 0.5054**. That is an independent confirmation
that the rank→PDS calibration is behaving correctly and is not an artefact.

| arm | mean rank | est. PDS | raw MSE | PDS pts | ≈net vs B0 |
|---|---:|---:|---:|---:|---:|
| **B0 shipped** | 16.75 | **0.8500** | 0.00461 | 0.0000 | **0.0000** |
| B1 F1 | 36.68 | 0.5547 | 0.00329 | −0.1028 | −0.1028 |
| B2 F3 | 29.05 | 0.6678 | 0.00357 | −0.0634 | −0.0634 |
| B3 pooled | 16.80 | 0.8493 | 0.00459 | −0.0003 | −0.0003 |
| B4 training-mean | 40.01 | 0.5054 | 0.00337 | −0.1200 | −0.1200 |
| B5 shrink 0.5 | 30.63 | 0.6443 | 0.00355 | −0.0716 | −0.0716 |

(B0 and B1 reproduce by construction — the calibration is fitted on them. B2-B5 are genuine predictions.)

**No arm beats B0.** B3 (the pooled donor branch) is essentially tied and marginally worse. B2 gains
MSE but gives back 0.063 points of PDS. B4 and B5 — the shrinkage family the advisors hoped would be
free — are the two worst arms on discrimination, exactly as the "constant vector has no
discrimination" argument requires.

**Caveat that keeps this from being final:** ≈net models PDS + MSE only. Fidelity, jaccard, nmae and
reach are not modelled, and B1 gained **+0.0243** on those four combined. If B2/B3/B5 gain a similar
+0.02 on DE terms, B3 could plausibly edge B0. That is the one live possibility, and it is exactly what
the official B2-B5 run will settle. It is also the *only* reason to still want those numbers.
