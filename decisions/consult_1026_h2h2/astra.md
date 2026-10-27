**My decision: if manifold fails, run Sidechain staged source progression. Optimize the composite to move up from rank 463. No GPU.**

The files now contain results newer than your brief. B2 has completed scoring; B5 and the final verification were incomplete when I checked:

| Local official metric | B0 shipped | B1 STRING | B2 F3 |
|---|---:|---:|---:|
| Expression MSE ↓ | 1.41072 | 0.95190 | 0.96653 |
| PDS ↑ | 0.85005 | 0.55469 | 0.65596 |
| Composite ↑ | **0.02155** | −0.05569 | −0.01302 |

These are **local SEARCH scores, not leaderboard scores**. Both manifold arms currently lose the composite. Finish the locked analysis and verification before recording the formal verdict. B0 (source reference: `h2h2_1026/results/official_B0.json`; not included), B1 (source reference: `h2h2_1026/results/official_B1.json`; not included), B2 (source reference: `h2h2_1026/results/official_B2.json`; not included)

1. **After a manifold loss: Sidechain staged source progression. Reliability weighting stays closed.**

   There is already a substantial negative experiment on reliability-weighted pooling: 24,593 rows, eight held-out families, **0.206% worse** family-macro expression MSE, and a confidence interval crossing zero. Recommending it again without identifying a materially different mechanism ignores the experiment ledger. Completed report (source reference: `source_weight_1025/REPORT.md`; not included)

   **Protocol:** one bounded source-addition ablation. Pin the documented Sidechain source progression, identify its first genuinely additional eligible source tier relative to our inventory, and add that tier to the existing predictor. Hold architecture, pooling rule, amplitude calibration, renderer, and scoring fixed. No simultaneous weighting or scalar sweep. If the tier is already represented or cannot be reproduced, stop before training.

   Select any necessary settings inside training-only folds. Evaluate with whole recipient families excluded from donors and fitting. The repeatedly inspected H1 SEARCH panel is development evidence now; reserve an uninspected confirmation panel before running.

   **Kill test:** require improved official expression MSE **and** composite, with a paired target-bootstrap composite lower bound above zero on confirmation. My minimum worthwhile composite gain is **+0.005 locally**, declared before running—not a predicted leaderboard gain. Failure ends this source-addition experiment; no expanding the source shopping list afterward.

   The reason to choose this is narrow: it tests additional information. Reweighting the existing information has already failed.

2. **A manifold win must survive the full score, trivial controls, and independent confirmation before packaging.**

   The 168/79 disjointness is correct. But the 79 must be absent from **all** fitted transforms and selection—not merely ridge fitting. Same-target donors from other families are legitimate under this protocol.

   Require the existing preregistered gates against B0, B4, and B5; exact B0 reproduction; all 79 targets retained, including missing-STRING fallbacks; and full verification. Then freeze one candidate and confirm it on an uninspected panel relevant to deployment. H1 target generalization alone does not establish transfer to unseen cellular contexts.

   **First invalidator: the rendered composite loses despite better response cosine or MSE.** That is already the point-estimate result for B1 and B2. B1 earns approximately **+0.00164** composite points from expression MSE while surrendering **−0.10282** through PDS. That trade is fatal. Do not move the gates to rescue it.

3. **The “below-chance direction” diagnosis is wrong. Do the bounded integrity check, then stop chasing an inversion.**

   Official fidelity is
   \[
   F=\frac{\text{correct-sign calls}}{\max(n_{\mathrm{pred}},n_{\mathrm{real}})}.
   \]
   Its fair-coin expectation, conditional on valid calls, is therefore \(0.5\,n_{\mathrm{pred}}/\max(n_{\mathrm{pred}},n_{\mathrm{real}})\), not universally 0.5. [Official definition](https://github.com/ArcInstitute/cell-eval2/blob/main/docs/vcc2026_metrics/vcc2026-metrics.md)

   The corrected local report puts B0’s coverage-matched expectation around **0.40**; **0.4929 was another panel’s result**, not chance. Thus 0.4487 does not establish below-chance performance. Nor does poor fidelity prove that a model is out of sample—provenance establishes that.

   Check gene/target alignment, units, control identity, and shipped-output reproduction. The report records B0 reproduction passing. There is no evidence here warranting an open-ended bug hunt.

   **What carries B0 locally is perturbation discrimination:** PDS contributes approximately **+0.12186** composite points; expression MSE contributes **zero**. Fidelity and Jaccard together subtract approximately **0.10998**. B0 can distinguish perturbations while predicting many individual gene directions poorly. Those are different capabilities. This decomposition does not establish what carries the hidden leaderboard score.

4. **Objective: (a), move up from rank 463.**

   Stop optimizing expression MSE in isolation. A large share of the leaderboard gap identifies potential headroom; it does not establish an accessible intervention. The observed manifold trade demonstrates that distinction.

   Stop reopening closed DE, conditioning, amplitude, and weighting lines without a new falsifiable mechanism. Stop treating cosine improvements as package evidence. Stop using the same inspected panel as fresh confirmation. There is no defensible forecast of 0.25 yet.

5. **GPU rental: no for the next 48 hours.**

   Finish frozen scoring, verification, and the bounded source experiment on CPU. The recorded operational failure concerns disk headroom, not GPU throughput. A GPU purchase does not resolve it.

**The frozen safety entry stays last. Nothing is submitted before the review deadline and authorizes it.**