# DECISION 2026-10-25 02:35 ET — what an unmasking gain is worth, and the NMAE break-even

## Why this note exists
`expr_unmask_1025` is running. Its headline is raw MSE, but the official composite is the unweighted
mean of six scaled components, so **one scaled unit of MSE buys 1/6 = 0.167 points of official
score and nothing more**. This note fixes the decision rule before the result lands, so a
large-looking MSE improvement cannot be mistaken for a win.

## The arithmetic to apply on receipt
Let `ΔS_mse` be the change in the **scaled** MSE component (higher is better), and similarly for
`ΔS_nae`, `ΔS_fao`, `ΔS_reach`, `ΔS_pds`, `ΔS_jac`. The run reports these directly, using the official
calibration anchors and the same score function as the incumbent. Then:

    Δofficial = (ΔS_pds + ΔS_mse + ΔS_nae + ΔS_fac + ΔS_reach + ΔS_jaccard) / 6

Do not convert raw MSE to a decision by eye. Use the harness's own scaled columns. If the local
scaled MSE is clipped at 0.0 (as it was for the incumbent in both earlier paired runs), treat the
MSE delta as **not independently measurable** and do not award it; record only the raw direction as
evidence.

## The cost that decides this test
The `v8a -> v8b` global amplitude experiment bought **+0.073 net scaled MSE-equivalent**
(+0.079 MSE, +0.028 PDS, +0.004 fidelity, +0.016 reach) and still lost **−0.006 overall** because
NMAE fell **−0.161**, three times the MSE gain (`decisions/20261025_0105`). Unmasking is directionally
the same kind of change: it restores response structure outside the DE call set, which is exactly
the region where NMAE's reference-gated denominator lives. So NMAE is the binding risk, not MSE.

**Pre-registered break-even for this arm:** if the other five components move by their measured
amounts, the arm is overall-positive only when

    |ΔS_nae| < ΔS_mse + ΔS_pds + ΔS_fac + ΔS_reach + ΔS_jac

With the other components expected roughly flat by construction (same call set, same layout), the
practical rule is **`|ΔS_nae| < ΔS_mse`**. If NMAE falls by more than MSE rises, the arm loses even
if raw MSE improves. Do not promote it on raw MSE alone.

## What this arm can and cannot establish
- **Can:** whether the expression-leg mask and the 50-control-PC projection suppress usable
  response direction at all. That is a mechanism question, and a clean paired answer is useful even
  if the answer is no.
- **Cannot:** establish that an official-score gain will follow. Local fidelity attenuates ~10x,
  Jaccard ~6x, reach amplifies ~1.5x, and MSE is the least transferable component. A local MSE gain
  is a reason to run a controlled official comparison, not a reason to ship.
- **Cannot:** reopen global amplitude, rho, pool_k, B/claim, source breadth, per-target shrinkage,
  global dense layout, or the learned-blend coefficient. Those are closed with measured evidence.

## Decision on receipt
1. If all three variants fail the frozen both-recipient rule, record a clean negative, close the
   expression-mask/deflation line, and spend the remaining night on the next untested mechanism or
   a written decision summary — not on another transform of the same banks.
2. If a variant passes on both recipients **and** its priced NMAE cost is below its MSE+PDS gain,
   write it up as a candidate for one controlled official validation comparison (Ricardo decides
   whether to spend a slot). Do not build or submit a package.
3. If it passes the raw-MSE rule but fails the NMAE break-even, record it as a **measured trade, not
   a win**, and put the gene-dependent amplitude question (per-expression-decile amplitude at fixed
   call set) at the top of the next queue.
4. Any outcome: no package change, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.

## Status
Written while `expr_unmask_1025` runs. Preregistration in the job directory;
no results inspected at the time of writing.
