**None; finish verification and hand off.** On the supplied evidence, there is no defensible new experiment with a demonstrated path from a local improvement to a better hidden final score.

The central problem is the missing **local-to-official transfer relationship**. The large discrepancy between the local incumbent and its official analogue does not establish that candidate rankings transfer—or that they fail to transfer. It leaves that relationship unknown. Consequently, even a statistically convincing local gain could optimize the wrong target. Trying additional variants on the same ablation risks selecting evaluation noise.

There is also a distinction between improving a candidate and improving the submitted entry: under **NEVER SUBMIT**, a local experiment cannot change the hidden score of the existing entry. Its value would be evidence for a later decision.

My ranking is:

1. **Finish the two pending audits and hand off.** Highest defensible operational value; expected direct score improvement is zero. Scientific risk is negligible. The production PASS and strict reseal already support integrity; the remaining value is resolving any independent audit findings and documenting exactly what is verified.

   **Minimal run:** zero new model fits, zero candidate searches, zero new datasets. Let the existing read-only audits finish; inspect their receipts and any discrepancies. Use existing artifacts and audit outputs, local CPU only, with no new large materializations while maintaining the storage reserve.

   **Exact gate:** both audit revisions must report PASS for their stated scope; package size must remain **2,595,584,000 bytes**, SHA-256 must remain **ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4**, and production coverage must remain **900/900 groups**, with no unresolved contradiction between receipts. Any failure means hand off the discrepancy as unresolved; it does not authorize changing the published package.

2. **Conditionally reopen a specific, newly evidenced implementation defect.** Potentially meaningful candidate upside, but no such defect is established in this brief. Risk is moderate: a superficially plausible “fix” can change intended behavior or exploit the local evaluator.

   **Trigger:** an audit identifies a concrete mismatch between intended scientific computation and actual execution, with an independently checkable expected result. Another tuning hypothesis is insufficient.

   **Minimal experiment if triggered:** reproduce the defect on the smallest affected case, then compare the incumbent and a corrected scratch implementation on a fixed, previously unused labeled holdout—**only if one already exists**. Reuse cached inputs; avoid a full package rebuild.

   **Gate:** the correction must fix the independently specified error, preserve applicable constraints, and improve the prespecified primary metric on that untouched holdout. That establishes a locally useful correction. It still does **not** establish official-score improvement without transfer evidence.

I would not rank another parameter, blend, source-breadth, or expression arm as a credible third option.

Reopening a closed lever requires evidence addressing **why it was closed**:

- A budget failure requires a demonstrated feasible configuration under the frozen budget—not a relaxed tolerance.
- A negligible or metric-negative result requires an independently supported mechanism and confirmation on unused evaluation data.
- A transfer-limited result requires historical paired local/official candidate comparisons under matched conditions, sufficient to estimate whether improvements transfer and with what uncertainty.

There is **no supportable numerical local-gain threshold** in the brief. Inventing one would manufacture confidence. If suitable paired historical results already exist, they could justify a prespecified transfer gate; collecting new official results conflicts with the current constraints.

The decision summary should therefore state: package integrity verified to the completed audit scopes, no new candidate justified, and hidden-score upside remains unestablished.