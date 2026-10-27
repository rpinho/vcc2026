**VERDICT — Pick (a): sign-preserving magnitude-profile recalibration.** Rank: **a > d (control/context-offset correction) > c > b.** Test whether redistributing response magnitude across genes adds transferable signal.

Three corrections materially affect that decision:

- **0.4929 is below 0.5023.** Fidelity/PDS ratios to the leader do not establish that your response vectors have the right direction: PDS measures retrieval rank.
- **1.0453 is not an exact leaderboard cutoff.** Your cached leaderboard (source reference: `mse_1025/inputs/leaderboard.json`; not included) contains **36 entries above it with positive MSE scores and 79 below it with zero scores**. Use each fold’s actual anchors, not the pooled regression.
- REACH measures the deepest confidence-ranked prefix achieving **90% directional purity** among reference-significant genes. NMAE also uses **reference-significant genes**, not just your calls. Neither is simply a call-budget objective. [Official definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md)

**EXPERIMENT — An eight-parameter shape correction at fixed response norm.**

Divide genes into **eight bins by input-control expression**. Learn one positive multiplier per bin from leakage-free training predictions; constrain multipliers to **[0.25, 4]**. Apply them to the predicted response, then restore each perturbation’s original response norm:

\[
\hat d'_p=\|\hat d_p\|\,\frac{D\hat d_p}{\|D\hat d_p\|}.
\]

Fit in the scorer’s normalized/log pseudobulk coordinates. Keep zero responses zero.

The mechanism is specific: **move response energy away from consistently overestimated genes and toward underestimated genes without globally shrinking the response.** Signs are preserved at this stage; rendered significance and PDS still require verification.

**Effect-size hypothesis:** a **0.10 reduction in raw MSE**. Using your fitted span only as a rough sensitivity, that contributes approximately **+0.014 composite when unclipped**, before changes in other metrics. If both values remain clipped, its MSE contribution is zero. That is a testable hypothesis—not an evidence-backed forecast.

**PROTOCOL**

1. **First 30 minutes: amplitude-versus-shape diagnostic.** For each perturbation, calculate actual response cosine, predicted/true norm ratio, and the oracle scalar
   \[
   \alpha_p^*=\max\!\left(0,\frac{\hat d_p^\top d_p}{\|\hat d_p\|^2}\right).
   \]
   Use a common input-control origin; separately report its offset from the reference control. Compare residual error before and after oracle scaling, by expression bin. Oracle scalars are diagnostic only.

   High cosine plus wrong norm indicates amplitude error. Large residual after optimal scaling indicates shape error or missing signal. **A transferable improvement from the eight-bin correction distinguishes fixable shape error from merely having poor predictions.**

2. **Fit and screen using training families only.** Cross-fit by entire cell family, including generation of the predictions used for calibration. Require the correction to improve mean response cosine by **at least 0.02**, with improvement in a majority of held-out training families. Freeze one recipe; no parameter sweep overnight.

3. **Spend four full evaluations:** incumbent and corrected candidate on **CD4 and one preselected second cell-line fold**, approximately **4–8 hours**. Reuse valid incumbent results where available. Match renderer, cells, seeds, and anchors. Record all six official components, raw MSE, and post-render response norms/cosines.

**KILL TEST**

Stop before rendering if the training-family shape screen fails. After rendering, reject unless:

- Raw MSE improves on **both** folds.
- Official composite improves on **both**, averaging **at least +0.01**.
- Post-render alignment improves, confirming that the proposed shape correction survived rendering.

**Beating the matched control predictor is the right scientific benchmark; it is not a sufficient scoring objective.** A clipped improvement earns nothing directly. Also, hidden raw MSE 0.973 establishes improvement over the nominal no-skill reference—not a measured 2.7% advantage over your exact control-emission pipeline. This experiment tests structural signal; it offers no basis for expecting leader-level performance overnight.