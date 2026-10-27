**VERDICT: Choose the decisive go/no-go experiment. Stop new modelling tonight. Build a verified six-metric comparison against the actual shipped predictor, and stop this DE-rendering programme if it cannot clear +0.02 composite on untouched validation. Do not build an 08:00 submission package.**

The biggest remaining risk is drawing another confident conclusion from incomparable measurements.

Your local diagnosis (source reference: `diagnostics/THE_REAL_DIAGNOSIS_2026-10-26.md`; not included) already records that the near-zero cosine came from CD4, while official scores came from a different panel. That limits the scope of the negative result. But its subsequent claim that an H1 within-line fold therefore matches the competition is also unsupported: Arc describes this challenge as prediction in previously unseen cellular contexts. Neither conclusion licenses another modelling campaign. [Official task description](https://arcinstitute.org/news/virtual-cell-challenge-2026)

**1. Adaptive budgets: theoretically useful, currently unsupported.**

Skewed reference counts establish that different targets have different optimal budgets. They **do not establish that your outputs identify those targets**.

Your local budget diagnostic (source reference: `diagnostics/NCONF_IS_UNPREDICTABLE_2026-10-26.md`; not included) reports essentially no predictive signal: best feature \(R^2=0.019\), best Spearman \(0.143\). The [decision record](DECISION.md) reports family-held-out \(R^2=-0.108\). On that tested distribution, **reject the adaptive-budget proposal**.

The correct objective would be to choose \(B\) maximizing predicted
\[
E\!\left[\frac{K_B}{\max(B,N_{\rm conf})}\mid\text{model outputs}\right],
\]
rather than making budget proportional to response magnitude. You have not demonstrated that this expectation is predictable.

“Already calls enough” kills the undercalling explanation. It does **not** kill replacing weak calls with better-ranked calls. That is the surviving hypothesis.

Two qualifications matter:

- Your 89.17% realization rate measures recovery of selected calls. Without counting significant **unselected** genes, it does not establish control of the complete emitted set.
- Fidelity measures directional correctness among adjudicable predicted calls; those genes need not themselves be reference-significant. Therefore your sign-precision evidence does not establish DE-set precision or Jaccard improvement. [Official metric definitions](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md)

**2. My actionable-headroom ranking: fidelity → Jaccard → NMAE → reach.**

This is an experiment-priority ranking, not four identified performance ceilings.

| Rank | Metric | What the numbers actually support |
|---|---|---|
| **1** | **Fidelity** | Raw incumbent approximately **0.492**, versus **0.705** top-50 directional precision. If that precision survived rendering and coverage penalties, the raw improvement would be **0.213**. With the published fidelity anchor spans, that is approximately **+0.108–0.127 composite** from fidelity alone. This is substantial—but conditional, not a forecast. |
| **2** | **Jaccard** | Current scaled **−0.016** leaves a plausible opportunity if the shipped package overcalls. But sign accuracy **0.752** does not tell you how often those genes belong to the reference DE set. You lack the overlap evidence needed to quantify this opportunity. |
| **3** | **NMAE** | **0.9171 versus null 1.0000** is an **8.29% raw error reduction versus null**. It does not demonstrate improvement over the incumbent. The failed inflation experiment supplies no evidence for another magnitude adjustment. |
| **4** | **Reach** | **0.186 versus 0.892** looks inviting, but reach demands a prefix with **90% correct directions**. Your strongest reported aggregate precision is **75.2%**, and accuracy on reference DE genes is **56.44%**. There is no demonstrated broad, sufficiently reliable prefix. |

Fidelity’s published anchor span is approximately 0.28–0.33; reach ranks reference-significant genes by prediction confidence and requires 90% directional purity. [Official metric reference](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md)

Coverage still matters after pruning. At \(B=50\), a target with \(N_{\rm conf}=153\) has fidelity at most \(50/153=0.327\), even with perfect signs. Applying the pooled 0.705 precision gives approximately **0.230**, not 0.705.

Also, **0.71 is an observed operating point, not a proven ceiling**. Conversely, the leader’s reach gap is not recoverable merely by issuing more calls.

**3. Replicates: a legitimate conditional route, not an explanation of the leader.**

First retract the numerical inference. For ordinary uncorrected squared-error geometry,
\[
R=1+a^2-2a\rho=(a-\rho)^2+1-\rho^2.
\]
Thus \(R=0.2325\) would imply a **lower bound** of approximately 0.876 on positive cosine under those assumptions—not identify the cosine. The official noise-corrected panel statistic does not directly satisfy that simple identity. Neither calculation establishes that the leader used measured replicates.

Arc permits external public and private training datasets, while withholding challenge perturbation measurements. Independent, legitimately obtained same-line measurements could therefore be useful; provenance and any applicable competition restrictions still need checking. [Official data policy description](https://arcinstitute.org/news/virtual-cell-challenge-2026)

The minimum viable test is:

1. Establish that an independent dataset actually covers the relevant **cell context and perturbations**.
2. Estimate responses using its own matched controls.
3. Evaluate transfer into an entirely separate experimental batch or study, with no shared cells or fitted validation labels.
4. Compare against the incumbent through the same renderer and scorer.

A random cell split from one experiment measures repeatability. It does not validate transfer to a new experiment. Unperturbed controls alone do not reveal perturbation-specific responses.

**I would not call this the only path to the top. There is insufficient evidence for that claim.** If qualifying data are absent, this route is unavailable tonight.

**4. Build one falsification report, with a binding stopping rule.**

Use the existing magnitude-ranked renderer as the single candidate. Compare it with the shipped predictor and the existing fixed-budget comparator. Do not launch another architecture or broad budget sweep.

The report must establish:

- The incumbent prediction path matches the shipped code, configuration, data access and rendering.
- Validation preserves the challenge’s held-out-context restriction. An H1 within-line result is supplementary evidence.
- All six metrics are scored from emitted integer counts, including unintended significant calls.
- Candidate selection is frozen before the final evaluation; uncertainty is estimated over perturbation groups, not individual genes.

**Pass:** candidate improves the validated six-metric composite by **at least +0.02**, improves mean raw fidelity by **at least +0.05**, and the paired 95% confidence interval for composite improvement excludes zero. An adaptive rule earns credit only if it also beats the fixed-budget comparator.

**Fail:** those conditions are not met. Stop this model-and-renderer programme and write up the result. If baseline reproduction or validation integrity fails, stop tuning and document that measurement failure; do not label it biological impossibility.

The leaderboard gap is **0.30765 composite**, requiring **1.8459 summed component points**. Even the optimistic fidelity calculation closes only part of it. The decision tonight is whether your existing signal supports a reproducible improvement—not whether another plausible story can fill the gap.