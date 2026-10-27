# DECISION — 2026-10-26T02:00Z (H1 result round)

Advice on record: `kimi.md`, `fable.md`, `astra.md`. Two of three returned before this was written;
Kimi is still running and will be appended when it lands. Measurement overrides advice where they
conflict, and each override is stated.

## The measured fact that drove this round

Official H1 fold, four arms, same frozen protocol and same run:

| arm | score_avg |
|---|---:|
| incumbent | **0.021555** |
| A-shared | 0.008394 |
| B-real | 0.011010 |
| B-shuffled | 0.009601 |

The incumbent wins every pairwise comparison. **Recipient-control conditioning is closed** as a
negative result. Independently verified; see `dirmap2_1025/REPORT.md`.

## Options considered, with arithmetic

**A. Replace the expression head with the control mean.** Rejected, unanimously and correctly.
Official raw MSE 0.973, scaled +0.061; control mean is raw 1.0 → scaled 0. Direct loss 0.061/6 =
**−0.010 composite**. Worse, `pds_cosine` (0.706 scaled) is computed from the same expression output,
so giving every target the identical control profile destroys perturbation discrimination too.
Astra: optimistically 0.16438 − (0.061+0.706)/6 = **0.03655**. Fable: same −0.010 conclusion.
*My earlier note calling this a "free fix" was wrong and is withdrawn in VCC_OVERNIGHT_LOG.md.*

**B. Manifold experiment.** Already running (`vcc_manifold_1026`). Fable computed EV ≈ +0.008
assuming ~30% pass odds; Astra correctly refused to assign a probability, noting it is not
identifiable without scored outputs. Both agree it has the better upside and should be allowed to
finish. It is the only candidate mechanism that could produce non-orthogonal expression direction
rather than rescaling a vector of length zero.

**C. DE-side effect-size gating** (Fable's proposal). Emit only genes with predicted |log2FC| > 0.8,
with a per-target budget drawn from the skewed reference count prior (median 31, p10 3, p90 565).
Rationale from our own measurements: sign accuracy rises to 0.636 above |log2FC| > 0.8, top-10
precision is 0.752, and realization delivers ~89% of chosen genes. This attacks fidelity, the
component where we are at chance (raw 0.4929).

**D. Gene-wise sign calibration** (Astra's fallback). Regularized logistic sign predictors from
incumbent signed LFC and gene identity, cross-fit by target, magnitudes and call sets frozen.

## DECISION

1. **Let the manifold experiment finish.** It is the only expression-direction mechanism left and is
   already running. Its own preregistered gate governs: pass = cosine ≥ 0.20 and MSE ratio vs
   do-nothing ≤ 0.90 in ≥4/7 families with the permutation control at or below do-nothing.
2. **Launch Fable's effect-size-gated calling experiment now**, on CPU, in parallel. C and B are
   independent (DE call sets vs expression basis) and the machine has headroom.
3. **Do not build any submission.** The provisional fold 5NeeoWn96lfm23hgExcm stays last and frozen.
4. **Adopt Fable's kill test as written:** the gated arm must beat the incumbent's raw
   direction_fidelity by ≥0.05 AND not give up more scaled reach than it gains in scaled fidelity,
   with a random-sign control sitting near 0.5. If it fails, gating is dead.
5. **Adopt Astra's caution:** no claim of composite improvement until a paired held-out comparison
   survives uncertainty and beats a shuffled-sign control.

## Ceilings, as stated by the advisors (planning judgments, not measurements)

- Fable: 0.24 realistic, 0.26 absolute if submission reopened; frozen official stays 0.16438.
- Astra: 0.20 prospective; reaching it needs 6 × (0.20 − 0.16438) = 0.214 additional summed
  component points.

Both are far below the leader's 0.47203. On current evidence the leaderboard top is not reachable
this cycle; the realistic objective is to move up from rank ~463.
