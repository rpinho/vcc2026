**Choose A, with a bounded calibration check inside the experiment.** It tests a documented defect in the fidelity denominator without requiring better response directions. B has greater theoretical upside, but no measured mechanism yet establishes that proxy directions help; C alone produces no candidate improvement. The earlier global budget sweep did not test target-specific allocation. Expect a modest gain at best—this will not close the MSE gap.

Use this **single, three-hour experiment**, with the following rules frozen before scoring.

1. **First 20 minutes: check score arithmetic.** Feed the archived official v5b/v9a component vectors through the aggregation code and require reproduction of their official overall scores within reported rounding precision. Check normalization, weights, clipping, and aggregation order. Do **not** fit a correction to those two points or require H1 to reproduce the official ranking: different evaluation populations can reverse the trade-off. If the arithmetic remains unexplained, continue measuring components, but declare any official-score implication unresolved.

2. **Freeze one public-source budget estimator.** For each target, estimate responder-gene counts using the evaluator’s DE definition, gene universe, and matched cell counts. Exclude the evaluation cell line and its entire experiment from source evidence. Take the median count across source contexts within each source family, then the median across families, rounded and bounded by the eligible gene count. Use **885 when evidence is missing**. No fitting against H1 labels, no alternative estimators after seeing results, and no 885 ceiling for supported targets—the upper tail matters.

3. **Run these arms on the identical 79 H1 targets:**

   | Arm | Layout and budget |
   |---|---|
   | Baseline | v9 confidence order, fixed 885 |
   | Candidate | Same order, predicted per-target budget |
   | Reference | Native v5b layout |
   | Negative control | Same predicted budgets randomly reassigned among source-covered targets, using a frozen seed |

   Hold the response head, amplitude, per-gene integer totals, B, controls, cell counts, seeds, and scorer constant. Verify those invariants directly. Record **realized scorer calls**, since requesting a budget does not establish that the renderer achieved it. MSE/PDS must remain identical within numerical tolerance.

4. **Predeclare the decision rule.** Candidate versus fixed 885 must achieve all three raw-component thresholds:
   - Fidelity **≥ +0.030**
   - Jaccard **≥ +0.005**
   - Reach **≥ −0.005**

   Additionally, require a positive lower bound for the paired target-bootstrap 95% interval on fidelity improvement. Against shuffled budgets, require fidelity **≥ +0.010**, with neither jaccard nor reach worse by more than **0.005**. These are experiment acceptance hurdles, not predicted gains. Native is a reference, not an arm selected by local overall.

5. **Use the remaining time for CD4 with the frozen rule**, if its existing pipeline fits the three-hour cap. Require positive fidelity and jaccard deltas and reach loss no worse than 0.005. A CD4 failure rejects a transfer claim; an unavailable CD4 run leaves an H1-only result. Report source-covered and uncovered targets separately: improvement on covered targets cannot be extrapolated across all 300 when only 139 have evidence.

**Abandon this version** if it misses any H1 hurdle, cannot preserve the invariants, or fails CD4 replication. Do not spend the remaining night rescuing it through estimator or threshold sweeps.

**The biggest way it could fool us is that generic budget shrinkage looks like successful prediction of target-specific response size.** The shuffled-budget arm catches that: it preserves the complete budget distribution, including its heavy tail, while destroying target assignment. If the candidate cannot beat it, the proposed mechanism has failed.

I would spend the hours on this experiment, but I would **not promise a material official-score gain tonight**. Even a pass establishes a plausible component improvement; the known population-dependent ranking reversal still limits what H1 can tell us about D/E/F.