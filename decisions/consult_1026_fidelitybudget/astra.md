**VERDICT: The budget approach is alive, but +0.053 is a conditional scenario—not a forecast. “Power, not sign” overstates the evidence. Start with B=100; choose the final budget using all six metrics after rendering.**

I found the completed **B=885** rendering result locally: **89.2% realisation**, fidelity **0.1767 → 0.5706**, unchanged PDS, essentially unchanged NMAE/Jaccard, but reach **0.2796 → 0.2397**. This establishes feasibility for that renderer and fold, not B=150 or leaderboard transfer. Its fidelity remains below the **local baseline 0.5866**. Result (source reference: `gate_1025/results/official_CD4.json`; not included)

**1. Arithmetic and conventions**

Yes, each metric has weight **1/6**; each context is scaled against its own anchors before aggregation. Using an averaged anchor span is only an approximation. [Official scoring guide](https://vcc-cli-wiki.virtualcellchallenge.org/)

Your displayed calculation gives the **new contribution**:

\[
(0.5892-0.512)/(0.301\times6)=0.04275.
\]

The **gain over current fidelity** is:

\[
(0.5892-0.4929)/(0.301\times6)=\mathbf{0.05332}.
\]

Thus **0.2177** is arithmetically consistent if the new raw fidelity transfers and everything else stays constant.

The fidelity formula is correct, but its definitions matter:

- `n_pred`: predicted-significant genes whose reference direction is defined, excluding the perturbed gene.
- `k`: matching signs among those genes. **Reference significance is not required.**
- `n_conf`: reference-significant count after target exclusion.

Consequently, 0.73 “precision” does **not** mean 73% correct DE discoveries. Nor does 0.5702 establish useful prediction against an imbalanced-sign baseline. [Direction implementation](https://github.com/ArcInstitute/cell-eval2/blob/main/src/cell_eval2/metrics/direction.py)

There is also an audit defect: fixed_budget.py (source reference: `diag_fid_1025/fixed_budget.py`; not included) labels a truth-budget oracle-prefix calculation `actual`, and filters reference-zero directions **before** selecting top B. Select using predictions alone, then let the scorer adjudicate. Until repaired, I would not certify that script’s curve as leak-free.

**2. Can other metrics erase the gain? Which B?**

**Jaccard:** \(J=I/(B+R-I)\), where \(I\) is significant-set overlap. At \(R=30\), even perfect reference recovery caps Jaccard at **0.60 for B=50**, **0.30 for B=100**, and **0.20 for B=150**. Your actual overlap is unknown; directional fidelity cannot supply it.

However, your current official Jaccard is already near its floor. Using the documented validation floors and your recorded scaled Jaccard −0.016, its remaining downside is roughly **0.011 composite**. Jaccard alone is unlikely to erase +0.053 on those anchors. [Metric policies](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/metrics.md)

**NMAE:** changing significance alone has **zero direct effect**. It evaluates predicted LFCs on a **reference-defined significant set**, regardless of predicted significance. Altering counts can change those LFCs. A **0.320 scaled NMAE loss** would erase the entire +0.0533 fidelity gain. [DE implementation](https://github.com/ArcInstitute/cell-eval2/blob/main/src/cell_eval2/metrics/de.py)

**Reach:** it is **not recall**. It measures how far the prediction-ranked reference-significant genes sustain **90% directional purity**. Rendering can reorder that list and hurt reach—as your completed run demonstrates.

**MSE:** preserving group totals does not guarantee preservation of the scored, noise-corrected MSE. Cell layout changes its sampling correction.

My provisional choice is **B=100**. From your curve, B=50→100 buys **0.0223 composite** through fidelity; B=100→150 buys only **0.0106**. The latter is easily consumed by collateral losses. There is insufficient evidence to declare any B composite-optimal.

**3. Realisation per unit compute**

My ranking is conditional on the current renderer:

| Rank | Lever | Main failure mode |
|---|---|---|
| 1 | Selective count redistribution to reduce spread | Wilcoxon has no adjustable dispersion parameter; you must change actual cells. Rank shifts, extra calls and lost MSE noise credit can offset gains. |
| 2 | Inflate selected LFCs | Cheap power, expensive errors: wrong signs and magnitudes damage NMAE and expression metrics. |
| 3 | Generate genuinely deeper counts | Helps sampling/zero noise only; costs density and memory. Simply multiplying every cell’s counts leaves normalized values—and Wilcoxon ranks—unchanged. |
| Unavailable | Increase cells per target | Official submissions require **400 cells per perturbation**. Extra generated cells are not an admissible power knob. [Submission requirements](https://vcc-cli-wiki.virtualcellchallenge.org/) |

Use the already-demonstrated layout mechanism first, with explicit preservation checks. Do not infer that halving LFC spread proves insufficient statistical power; effect-size spread and sampling uncertainty are different quantities.

**4. Reachable signal without strong expression cosine**

**NMAE:** yes, potentially. Absolute error on reference DE genes is a different objective from genome-wide squared expression error. But budget expansion itself supplies no NMAE improvement, and matching the reference LFC standard deviation is not a valid objective.

**PDS:** yes. It measures **relative retrieval rank**, not absolute matched cosine. A small correct-target advantage over competing targets can score well despite poor overall alignment. Your existing scaled PDS around **0.706** already demonstrates this. Preserve it; count redistribution with unchanged group totals preserves its input. The direct target-knockdown shortcut is excluded. [PDS definition](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/metrics.md)

Also, the leader’s **scored MSE** does not by itself imply cosine 0.875: that inference needs an appropriate raw error ratio and assumptions about amplitude and noise correction. This does not reopen your empirically dead modeling routes.

**5. One next experiment**

Run **one count-rendered budget sweep: incumbent versus B=50,100,150**, keeping the predictor, cell count, row depths and group gene totals fixed. Repair selection first. Use CD4 to choose B, then freeze it and evaluate on one untouched held-out family with the exact six-metric scorer.

**Pass:** on that confirmation family, **≥80% intended-call realisation**, **≥+0.02 composite**, and a paired target-bootstrap **95% lower bound above zero**, using matched anchors. Report unintended calls and all six component changes.

**Fail:** otherwise, reject deployment of that budget change. High realisation or higher fidelity alone is not a pass.