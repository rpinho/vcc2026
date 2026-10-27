**Decision: B0 remains shipped unless something beats its full official score. Finish B2–B5, then close the direction-search program. If B2 fails, Monday’s sole new experiment is amplitude-only calibration of B0—not reliability weighting, not Sidechain. No GPU.**

### 1. “Contest expression MSE” is not the objective

The advisors were **half right**: stop the failing DE direction work, but stop calling MSE the objective.

The objective is **official composite improvement with B0’s PDS protected**. DE fidelity and jaccard headroom are not expected points; they are merely points you would recover *if* you reached baseline. Three failed DE programs and B1 provide no evidence that you can cash that headroom without destroying PDS.

B1 earned **+0.02558 outside PDS** and paid **−0.10282 in PDS**. That is the coupling problem. Until an experiment decouples DE improvement from PDS rotation, the large DE headroom is not a plan.

MSE is worth pursuing only when it is nearly free.

### 2. PDS ≥0.80 plus better MSE is **not** a B0 replacement

Using the observed B0→B1 exchange rate:

- Dropping PDS from 0.8500 to 0.8000 costs approximately
  **0.10282 × 0.05 / 0.2953 = 0.01741 points**.
- Maximum MSE recovery is **+0.00164**.
- Therefore B2 at 0.80 starts roughly **0.01577 points behind**, before counting any DE changes.

With DE otherwise unchanged, B2 needs approximately **PDS ≥0.845 even if it captures all the MSE headroom** just to break even. At 0.80, it needs substantial additional DE gains.

**Packaging gate:** the frozen candidate must beat **0.021555 on the full official score**, with B0 winning ties. No special gate based on “good MSE” or “PDS roughly held.” B4’s superior response-space MSE is irrelevant until its official composite survives.

If B2 clears that gate, package it. Otherwise keep B0 and stop extending the manifold program.

### 3. Yes—run amplitude-only calibration Monday, but timebox it brutally

This is the one remaining cheap hypothesis that directly tests decoupling.

For a prediction \(y=\mu+\alpha d\), any **positive** \(\alpha\) preserves the response-space PDS cosine:

\[
\cos(\alpha d,t)=\cos(d,t).
\]

Meanwhile, MSE is quadratic in amplitude, so the LOO estimate \(\alpha^\*\approx0.10\) is a credible reason to expect substantially lower raw MSE without rotating B0’s direction. B5 at 0.5 is only one coarse point on that curve; it does not test the proposed optimum.

But **do not ship 0.10 directly**. The renderer and amplitude-sensitive DE calls—especially jaccard/reach—can break the apparent invariance. Fit the fixed per-family scalar using **LOO on the 168 only**, check the full component suite rather than minimizing MSE alone, and give the selected candidate **one official evaluation**. It ships only if it beats B0’s composite.

It jumps ahead of reliability weighting and Sidechain because it is cheap and isolates amplitude from direction—not because +0.0016 is strategically large. One CPU-scale calibration, one official score, then either ship the improvement or kill the idea.

### 4. Objective: **rank**

Not 0.25. Not a component.

Even repairing fidelity and jaccard to baseline and collecting the MSE baseline headroom, while retaining everything currently above baseline, gives only approximately:

\[
0.021555+0.0799+0.0301+0.0016 \approx 0.1332.
\]

So 0.25 requires substantial above-baseline performance for which you currently have no evidence.

**STOP:**
- MSE-only selection and celebration;
- new manifold/direction arms after B2–B5;
- effect-size gating, recipient conditioning, and further DE work without fixed-PDS evidence;
- reliability weighting and Sidechain for this cycle;
- tuning anything on the 79.

### 5. GPU: **no**

You are not capacity-bound. You are objective- and coupling-bound. More GPU would mostly generate more candidates that improve secondary metrics while taxing PDS.

**Ship rule: highest full official composite, B0 on ties. Lower MSE with a lower composite is a worse entry. Full stop.**