**Keep B0 deployed. Drop “contest expression MSE” as the objective. Monday’s experiment is one source-fitted, amplitude-only calibration, judged on the full composite. No GPU.**

There is newer evidence on disk: B2–B4 have finished. I verified matching targets, controls, reference, anchors, and scoring configuration.

| Arm | Composite | PDS | Raw MSE |
|---|---:|---:|---:|
| B0 shipped | +0.021555 | 0.850049 | 1.410715 |
| B2 manifold F3 | **−0.013018** | **0.655956** | 0.966530 |
| B3 pooled | **+0.025615** | **0.858650** | 1.409414 |
| B4 training mean | −0.099178 | 0.500000 | 1.022118 |

B2 result (source reference: `h2h2_1026/results/official_B2.json`; not included), B3 result (source reference: `h2h2_1026/results/official_B3.json`; not included), B4 result (source reference: `h2h2_1026/results/official_B4.json`; not included).

**1. The advisors chose the wrong optimization objective.**

Replace it with **improving the official composite while preserving demonstrated discrimination**. MSE is a diagnostic and a component, not the mission.

Two arithmetic corrections matter:

- Reaching the MSE anchor baseline earns **zero** points. B1’s +0.00164 comes from beating that baseline. That is the demonstrated gain, **not MSE’s theoretical ceiling**.
- Your component table sums to **+0.13153 − 0.10998 = +0.02155**, not +0.1592 − 0.1376.

The DE deficits are opportunities to investigate, not recoverable cash. Fidelity rewards both correct directions and sufficient coverage; suppressing calls can improve precision while destroying fidelity. Three failed interventions reject those interventions—not every possible DE improvement. The [official definition](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md) makes that coverage penalty explicit.

B2 now strengthens the conclusion: **all four DE components improved**, together earning +0.0330, but PDS cost −0.0676. Its substantial raw MSE improvement earned **nothing** after scaling.

**2. B2 is rejected. B3 advances to confirmation, not deployment.**

Even hypothetically, **PDS ≥0.80 is not a packaging gate**. Dropping from 0.85 to 0.80 costs approximately **0.0174 composite points**. With only B1-sized MSE gains and unchanged DE, break-even PDS would be approximately **0.8453**.

B3’s observed improvement is +0.00406, mostly from PDS. It is the only completed challenger here that deserves promotion testing.

My packaging gate:

- Freeze the candidate before further evaluation.
- Confirm a positive **full-composite** improvement against B0 in a second held-out context, with no tuning on either evaluation panel.
- Require the combined prespecified comparison to have a positive paired 95% lower confidence bound; uncertainty estimation must respect panel-dependent PDS and expression scoring.
- Verify the packaged artifact reproduces the evaluated predictions.

The repeatedly inspected 79-target SEARCH panel cannot supply fresh confirmation. Until that gate passes, B0 stays.

**3. Yes: run amplitude-only calibration ahead of reliability weighting and Sidechain.**

The hypothesis is **not** “PDS and MSE are inherently anti-correlated.” They measure different properties. For a nonzero effect vector and positive scalar,

\[
\cos(a\hat\delta,\delta)=\cos(\hat\delta,\delta).
\]

Thus amplitude can improve squared error without changing cosine retrieval **if the final scored effect remains exactly collinear**. The scorer uses normalized, logged pseudobulk effects relative to reference controls; the invariance applies there. [Official PDS definition](https://raw.githubusercontent.com/ArcInstitute/cell-eval2/main/docs/vcc2026_metrics/vcc2026-metrics-brief.md).

That condition is not guaranteed by multiplying an upstream response. Control mismatch, count normalization, clipping, integer emission, and sampling can change the final direction. DE significance also changes even when PDS holds.

Monday’s run is therefore:

- Scale the **entire shipped perturbation effect relative to control** by one positive scalar per predefined family. Do not interpolate toward the training perturbation mean.
- Fit using leakage-free out-of-fold predictions on the 168; exclude each target from every target-dependent fitting step. Use pooled fallback for poorly supported families.
- Select calibration using the **full official composite**, with scalar 1 available. Freeze one configuration before evaluating the 79.
- Measure all six components and the rendered direction change.

**Do not hardcode 0.10 from the 79.** It is diagnostic evidence, not an admissible fitted parameter.

The earlier CD4 scalar experiment (source reference: `mse_1025/REPORT.md`; not included) already reduced MSE while reducing fidelity. That lowers expectations. It does not settle the missing six-component comparison. One bounded experiment is justified; an amplitude sweep on the 79 is not.

**4. Pick (a): rank.**

Operationally, pursue validated composite gains across contexts. A component win is evidence, not success. **0.25 has no supporting trajectory:** merely repairing fidelity and Jaccard to baseline would bring B0 to roughly **0.1315**, assuming everything else held.

Stop:

- MSE-only selection and promotion from response-space proxies.
- Further manifold variants from this failed direction.
- Repeating effect-size gating or recipient conditioning without a new, testable mechanism.
- Reliability weighting and Sidechain work during the bounded calibration test.
- Treating local SEARCH improvements as demonstrated leaderboard gains.

If calibration fails, close that experiment and retain the validated incumbent. Do not manufacture another shrinkage variant.

**5. GPU: still no.**

The immediate uncertainty is whether a cheap calibration preserves useful discrimination and improves the actual score. More compute does not answer a different question. B3 confirmation and one CPU calibration experiment are the next moves.