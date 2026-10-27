 **VERDICT:** (a) Per-gene magnitude/shape recalibration.

**EXPERIMENT:** Gene-wise Responsiveness Calibration (GRC). Learn a per-gene scalar $c_g$ that corrects systematic shrinkage/expansion of predicted log-fold-changes relative to the empirical distribution in training data.

**PROTOCOL:**
1. **Compute empirical responsiveness:** On training data (source cell lines), for each gene $g$ calculate $s_g = \text{median}_{p}|\text{LFC}_{g,p}|$ (robust measure of how much the gene actually varies across all perturbations).
2. **Compute model responsiveness:** On the same training split, calculate $\hat{s}_g = \text{median}_{p}|\hat{\text{LFC}}_{g,p}|$ from your current model’s continuous predictions.
3. **Derive factor:** $c_g = s_g / \hat{s}_g$ (winsorize at 0.2 and 5.0 to avoid volatile low-count genes). If $\hat{s}_g \approx 0$ but $s_g > 0$, set $c_g = 3$ (arbitrary boost for silent-but-responsive genes).
4. **Apply:** For the target cell line, multiply the continuous predicted expression of gene $g$ by $c_g$ (additive in log-space: $\hat{\mu}'_g = \hat{\mu}_g \cdot c_g$).
5. **Evaluate:** Run one full held-out fold with the v9 integer renderer.

**MECHANISM:** Your fidelity (0.49) and PDS (0.82) prove you have the **right direction** but wrong magnitude structure. Global amplitude scalars fail because they uniformly push genes across the DE significance threshold, trading MSE for REACH (what you saw). Per-gene factors fix the “shape” by amplifying the specific genes that truly vary (improving REACH) while dampening house-keeping noise (improving MSE), without destroying the direction signal.

**EXPECTED EFFECT SIZE:** Raw MSE drops from ~1.0 to ~0.65–0.75 (scaled gain +0.25 to +0.30). REACH improves from 0.24 to ~0.35–0.40 as responsive genes cross significance. Composite gain estimate +0.06 to +0.10 (200–300 rank jump).

**KILL TEST:** If on the held-out fold MSE does not fall below 0.90 (i.e., fails to solidly beat the control reference of 1.0) **OR** if REACH drops (indicating you amplified noise), kill immediately.

---

**ANSWERS TO SUB-QUESTIONS:**

1. **Highest-value experiment:** GRC as above. It attacks the 45% MSE gap directly while collateralizing the 13.5% REACH gap, unlike REACH-only hacks that poison NMAE.
2. **Is “beat control” the right target?** No—it is a trap. The metric is log-linear in the regime $0.2 < \text{raw MSE} < 1.0$ with a hard floor at 1.045. Beating control by 2% (0.973 vs 1.0) yields only +0.023 scaled points, essentially zero. The leader at 0.235 is in a different statistical regime; we need **structure**, not shrinkage. GRC provides structure.
3. **Diagnostic for shape vs amplitude:** Stratify genes by their mean expression in control samples. Compute **mean signed error** $(\hat{\text{LFC}} - \text{LFC})$ per stratum. If error is uniformly positive/negative across all strata → global amplitude error. If error flips sign (e.g., overpredict for low-expr, underpredict for high-expr) → shape error. **Fix shape first;** amplitude is a scalar you can tune in post-processing, but shape determines which genes enter the DE set.
4. **Ranking:** (a) > (b) > (d) > (c).
   - (c) Shrinking LFCs is suicide: it collapses REACH (already your second-worst metric) and gives diminishing returns on NMAE.
   - (b) Calling more genes without fixing magnitude profile adds low-precision noise (your partial run already showed this).
   - (a) Fixes the root cause evidenced by high fidelity + high MSE.