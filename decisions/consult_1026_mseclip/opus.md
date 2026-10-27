# Verdict: pick (d). Model the shared perturbation response explicitly

## VERDICT

Your model's perturbation-specific delta has roughly zero correlation with the truth in expression space. The evidence is that the optimally amplitude-calibrated arm landed exactly on the control reference. If the scale is fit optimally, MSE ≈ (1 − r²) × control MSE. Landing on control therefore means r ≈ 0. So the problem is **shape, not amplitude**.

My hypothesis is that the leader's 0.235 is not mostly perturbation-specific structure. It is most likely the **shared perturbation response**: the offset common to essentially all perturbed cells versus control (targeting/stress/library/batch effects). "Control + mean perturbation delta" is the known baseline that beats most models in Perturb-seq benchmarks. Your pipeline predicts around the control mean, so it pays the full shared offset on every target. That one vector would explain three things:

- why "predict control" scores about 1.0;
- why you sit near 1.0;
- why 61.7% of the board is parked at zero.

It also feeds the DE side. True DE-vs-control sets contain the shared genes for every target. That should lift REACH, NMAE and Jaccard together, through calls that are correct by construction rather than gate-widening calls of low precision.

## Answers to your four questions

**Q2. Is "beat control" the right target?**
Yes as a target, no as a method. Shrinkage cannot get past control, as your item 2 shows. A correct mean shift can.

Measure local skill against the local control (1.1087), not against 1.0. The local/panel disagreement is plausibly a fold-to-fold difference in how large the shared component is relative to your noise.

**Q3. Wrong amplitude vs wrong shape?**
Run this per target on pseudobulk deltas (no rendering needed):
- Compute r = cosine(predicted delta, true delta).
- Compute slope = ⟨predicted, true⟩ / ‖predicted‖².

Interpret it as follows:
- **Amplitude problem:** r is high and slope ≠ 1.
- **Shape problem:** r is near 0.
- **Expression-weighting problem:** recompute r with gene-wise z-scoring. If r rises, the error is concentrated in highly expressed genes.

Your data already implies r ≈ 0, so fix shape first. The cheapest large shape fix is the shared component.

**Q4. Ranking:** (d) shared response > (b) > (a) > (c).
- (a) is capped at control once r ≈ 0.
- (c) trades REACH for NMAE.

## EXPERIMENT

Predict, in log1p space:

pred_p = ctrl_target + α · s + β · specific_p

**Step 1: estimate the shared vector s.**
- If any target-cell-line perturbed cells exist in training, set s = mean over training perturbations of (pert mean − control mean) in that line. Use this whenever possible.
- Otherwise, average this vector across source lines. Apply it to the target control in log space, which makes it multiplicative in counts.

**Step 2: define the other terms.**
- specific_p = incumbent delta_p minus its own mean across targets.
- Grid: α ∈ {0.5, 0.75, 1.0, 1.25} and β ∈ {0, 0.5, 1}.
- Render through the existing v9 integer layout unchanged.

## PROTOCOL

**Phase 0: offline pseudobulk check (CPU minutes, zero renders, about 30 min).**
On both local folds, compute normalized MSE skill (1 − MSE / MSE_control) for:
1. control;
2. incumbent;
3. oracle leave-one-target-out shared delta from fold truth (the ceiling);
4. transferred s, across the α × β grid;
5. the shape diagnostic from Q3.

Also compute the fraction of true-delta variance explained by the shared component.

**Phase 1 (renders 1–2, fold A).**
- Render the best (α, β) from Phase 0, and the same α with β = 0.
- Score both with the official scorer.

**Phase 2 (renders 3–4, fold B).**
- Confirm the winner from Phase 1.
- Also render α × 1.25. Extra amplitude helps Wilcoxon calls on the integer counts, so check whether it trades MSE for REACH favourably.

**Phase 3 (renders 5–6, spare).**
Re-run if PDS drops. Shared offsets can compress cosine-based discrimination, so use β = 1 to preserve the specific term.

**Expected effect (conditional on the shared component being large).**
- If the oracle shared delta explains 40–70% of delta variance, raw MSE moves from about 0.97 to 0.45–0.65.
- That is roughly +0.35 to +0.5 on the scaled MSE component, or about +0.06 to +0.08 composite from MSE alone.
- Add maybe +0.02 to +0.04 from REACH and NMAE.
- Plausible landing: about 0.25–0.28 versus 0.164 today.

## KILL TEST

| Stage | Condition | Action |
|---|---|---|
| Phase 0 | Oracle LOO shared delta gives skill < 15% vs local control on both folds | Kill. The leader's gain is perturbation-specific; pivot to (b). |
| Phase 0 | Oracle is strong but transferred s gives skill < 5% on both folds | Kill the transfer route. Effort goes to estimating s, not to rendering. |
| Phases 1–2 | Official composite fails to beat incumbent by ≥ 0.02 on both folds | Do not adopt. |
| Phases 1–2 | PDS falls by > 0.03 and β = 1 does not recover it | Do not adopt. |
| Phases 1–2 | Fidelity falls below 0.47 | Do not adopt. |

Adopt only if both folds agree in sign. One fold winning is the coin-flip pattern you have already been burned by.