# From 0.0877 to 0.16438, and then learning why we could not improve

Our best recorded Virtual Cell Challenge 2026 validation score is **billmA
0.16438**, up from **0.0877** for our first transfer model. Its recorded rank was
about 463 of 1,580. We did not reach the strong teams. After validation closed on
October 22, we found local improvements, large failures and several mistakes in
our reasoning, but no sufficiently credible replacement for billmA. Final scores
are hidden. Keeping the incumbent is a decision under uncertainty, not evidence
that every alternative would lose on the final panel.

This account is for other teams who might recognize the same failure modes. The
[official ledger](OFFICIAL_SCORES.md) records submitted validation scores. Later
numbers come from local public-data experiments or vector diagnostics. Running
the official scorer locally does not make those numbers official results.

## What the submissions taught us

The first model, **v1, 0.0877**, put public same-target signatures on a recipient
control mean. It established some perturbation discrimination, but had zero
scaled MSE and negative fidelity. **v2b, 0.0617**, transported real control cells
with full-strength signatures. More realistic-looking cells did not compensate
for inaccurate transferred responses. **v3b, -0.2482**, added a common
average-perturbation shift and shrank target signatures. Fidelity fell to -1.539.
That was an expensive demonstration that a plausible common shift can be
disastrous when imposed on every target.

We should have reproduced stronger public methods earlier. **AtlasShift** gave
us the useful response-transfer and integer-count foundation; **Sidechain**
informed shrinkage and emission work. **v4a, 0.1438**, using AtlasShift anchors,
was the largest early advance. **v4b, 0.1243**, combined dual summaries, adaptive
shrinkage and a guide, but lost ground. Combining reasonable modifications did
not make them jointly beneficial. Attribution and source-license boundaries are
documented in [THIRD_PARTY.md](THIRD_PARTY.md).

**v5a, 0.1499**, calibrated spread to the scorer's noise-credit budget. Scaled MSE
finally moved above zero, to 0.078, while PDS stayed around 0.678. **v5b,
0.1528**, added cross-source sign agreement and removed the top 50 control axes.
It improved modestly, with MSE 0.082. Those gains taught us that count rendering
and scorer correction terms mattered alongside response estimation. They did
not establish that the underlying biology was accurate.

The next submissions exposed the tradeoffs. **v8a, 0.1470**, changed layout on
the same v5b totals: reach rose from 0.077 to 0.162, while MSE fell to 0.016.
Preserving group totals did not preserve the variance correction. **v8b,
0.1412**, increased amplitude: PDS reached 0.708 and MSE 0.095, but NMAE collapsed
from v8a's 0.102 to -0.059. Improving two members could lose the composite.

**v9a, 0.16135**, used confidence-ordered DE layout on v5b totals and became the
new reference. **v9b, 0.12049**, showed how fragile that success was; its component
scores were not retained, so we cannot give a reliable component-level
explanation. **slim_v5b, 0.08833**, sampled real controls with excessive spread,
roughly 6.5 times its intended noise credit. It got zero MSE and NMAE -0.276.
A diagnostic estimate had been treated as a sufficient gate instead of checking
the emitted counts.

Finally, **billmA, 0.16438**, blended a response MLP with v5b amplitude and v9
layout. PDS reached 0.706; MSE was only 0.061. This was our reconstruction of a
publicly described idea, not the original author's model or score reproduction.
Its gain over v9a was about 0.003. v6 was a generalized builder without a recorded
official score; v7 attempts were refused after the phase changed. Neither belongs
in a table of successful scored submissions.

## Why local progress stopped supporting a new entry

The [final entry history](README.md#final-phase) includes a direction variant
and restoration of the billmA safety recipe after that variant lost about 0.008
locally. The current last entry remains the safety recipe. With the validation
server closed, there was no feedback path for testing a new official improvement.

That alone would not prohibit a well-supported local change. The larger problem
was that our harness could not reliably order close historical recipes. One
seven-recipe audit got 18 of 21 pairs right, yet reversed important v4/v5
comparisons. Four recent layout/amplitude surrogates got three of six pairs
right, with Spearman zero. Missing exact replays, different emission paths,
source-family holdouts, target distributions and anchors all limited the comparison.
We had mathematical scorer parity without empirical candidate-ranking calibration.

The [component audit](calib3_1026/REPORT.md) made that failure concrete. On the
layout contrasts, local fidelity deltas were **9–12 times** official deltas;
Jaccard was **7–9 times**; reach was only **0.6–0.75 times**. DE signs agreed,
but their relative magnitudes reversed the composite. Our initial statement
that the harness got the direction wrong was too coarse: it got the overall
direction wrong while preserving the component directions.

The likely mechanism was panel composition. Strong, under-called targets were
much more prevalent locally. A fitted weak/mid/strong mixture explained some
contrasts, but effectively rested on one official layout comparison. Its
leave-one-contrast-out mapping got only one of four official signs right, or
two without MSE. It was useful for exposing uncertainty, not for predicting points.

MSE was worse. The local clamp often returned zero, hiding changes. Removing
the clamp did not repair transfer: the amplitude contrast gave **-0.001212**
locally versus **+0.0791** officially. The
[instrumentation study](h2h2_1026/scratch/mse_instr/REPORT.md) tried capped,
uncapped, trace-cap, budget and bulk-distance formulations. None passed its
sign/ratio gate. Numerical replay passed; the instrument failed. We should have
made that distinction much earlier.

## Lines we closed, and the evidence

The morning [H1 experiment](RESEARCH_UPDATE_20261026.md) compared STRING and
donor/manifold candidates after fitting 168 targets and evaluating 79 disjoint
targets. B1 improved raw MSE **1.4107→0.9519** but lost the composite
**0.02155→-0.05569**. B2 also lost. B3's +0.00406 composite point gain had a
paired interval crossing zero. Better expression distance was not sufficient.

The [afternoon tests](RESEARCH_UPDATE_20261026b.md) repeated that lesson. N1
improved pre-render direction but lost rendered PDS **0.8500→0.8397** and
composite by 0.00881. P1a's gene-wise magnitude calibration cut raw MSE by about
32%, yet PDS fell to 0.7970 and composite to 0.00579. Scaling by 0.75, 1.25 or
1.50 supplied no positive composite interval; larger amplitudes increased raw
MSE sharply. These failures closed the tested direction and magnitude recipes.

Renderer seed replicates moved the composite by up to 0.00076. Uniform confidence
lost 0.00543. Doubling calls to 1,770 gained **0.03449 locally**, with a positive
paired interval, but its apparent advantage depended heavily on strong targets.
Calibrated projections shrank toward the transformed seed floor. We retained
D1770 as an interesting diagnostic, **not a known official loss**, and did not
promote it.

The [R9 follow-up](h2h2_1026/results/r9_gate.md) found a top-20 reach head inert
(-0.00002), a 250-call budget badly worse (-0.07007), and detectability selection
worse (-0.00280). The latter two also failed validity checks. The
[sign audit](signaudit_1026/REPORT.md) rejected an easy renderer-blame story:
flipped signs agreed with truth more often than retained signs, and the proposed
reference-null down-sign bias was absent. The response head itself was weak.

Selective scatter had already lost roughly 0.10 locally with biology held fixed;
amplitude recalibration did not rescue it. Broader family-balanced billmB training
traded a small H1 improvement for worse CD4 fidelity and reach. The tested
DE-call and retraining combinations were closed, rather than rebranded as
unexplored opportunities.

The [decoupled-amplitude probe](dm_1026/REPORT.md) tried whole-cell depth
multipliers that preserve CPM vectors. Even continuous reweighting barely moved
amplitude along the response. Gene-wise depth-tilted placement remained unbuilt:
it could change Wilcoxon ranks, and its purported MSE upside could not be measured
reliably. An algebraically plausible emitter was not an established improvement.

## What the biology diagnostics actually found

The [training168 decomposition](bio_gap_1026/REPORT.md) was sobering. B0 had mean
cosine **0.084**, despite retrieval PDS around 0.856. Normalized SSE was 1.642;
58.3% of error was attributed to direction and 41.7% to magnitude. Retrieval can
identify the right target among alternatives while predicting its actual response
poorly. Integer rendering changed SSE by only 0.003% in this audit.

Direct source rows covered 161 of 168 targets and carried **99.5% of error**.
Perfectly fixing missing rows would therefore buy at most **0.5%** on this panel.
That ceiling concerns missing-row coverage only. Better donor contexts for
already-covered targets could still help; these experiments do not bound that gain.

Our working conclusion that “every SSE gain was shrinkage” needs a precise scope.
The large distance gains in the initial correction/shared-response campaign did
not establish new aligned signal. Global scaling reduced SSE by about 40% while
leaving just 1.74% of centered spread. Family correction reduced SSE about 10.5%,
but spread fell to about 80% and unrescaled normalized cross fell
**0.1154→0.1086**. Adding a shared mean plus shrunken B0 cut normalized SSE to
0.943 while PDS collapsed to 0.523. Under a spread-limited noise credit, losing
energy can erase an apparent distance improvement.

We then asked about aligned signal at matched centered spread. In the
[in-context experiment](bio_gap_1026/align/REPORT.md), kNN using other targets'
true H1 responses raised cross **0.1154→0.4230**; ridge/RRR reached 0.3984.
Their PDS dropped substantially. These are informative diagnostics, but D/E/F
provide controls, not those recipient perturbation outcomes.

The [leave-one-context-out experiment](bio_gap_1026/loco_map/REPORT.md) excluded
each recipient from training labels and sources. Full maps failed the PDS rule
in all five contexts. Small bounded blends passed a cross-CI/mean-PDS rule in
five of five, preserving a real, modest transferable result. Against full-axis
H1 B0, the best control-conditioned blend raised cross only **+0.0027**;
PDS changed -0.0023, with an interval reaching -0.0052, missing stricter
noninferiority.

That does not license an official score conversion. Cross uses truth-energy
units, and matching centered spread does not fix uncentered energy. The proposed
net +0.0004 overall was an unsupported scenario. Large gains were demonstrated
using in-context truth; we did not prove all larger future gains require it.
These limitations matter as much as the negative result.

## Advisors helped most when evidence corrected them

Fable, Opus, Astra and Kimi participated across the campaign; the afternoon
[rounds 5–12 and decision history](decisions/consult_1026_amp/DECISION.md) preserve
the later record. Kimi participated in earlier consultations, not this series.
The useful process was explicit proposals, fixed tests and written corrections.
Model agreement itself was not validation.

Earlier, Kimi's adaptive-budget suggestion assumed response magnitude predicted
DE counts; measurement gave R² 0.006 and Spearman 0.054. Another suggestion that
global LFC inflation would fix NMAE failed measurement. Astra caught a supposed
family holdout whose grouping variable was unused and an invalid inference of
cosine from MSE. Fable flagged incompatible official and local metrics, prompting
inspection of the compared contexts and renderers.

Later, advisors repeatedly proposed validation feedback that no longer existed.
We corrected that premise. We also withdrew claims that pure scaling was
rendered-PDS-invariant, and that fixed pseudobulk implied fixed raw MSE: the
variance correction depends on cell layout. The R7 rejection of D1770 as an
automatic official loss was revised after component calibration. In round 12,
Astra correctly distinguished matched centered spread from matched total energy
and rejected a confident score forecast. Even the final advice needed checking.

The process also cost time. Several rounds revisited closed amplitude and
renderer mechanisms under new descriptions. A second panel sharing the same
instrument bias was sometimes mistaken for independent confirmation. The
79-target development panel was inspected repeatedly; target bootstrap intervals
could not restore its independence or measure uncertainty over unseen contexts.

## What we would do differently

Start from the strongest reproducible public baseline, with attribution and
licenses settled. Save exact recipes, emitted-count invariants, anchors and all
six components for every official submission. Preserve the evidence needed for
later diagnosis, including failures, without redistributing restricted data.

Before optimizing, require the local harness to reproduce important historical
rankings and component tradeoffs. Keep a genuinely untouched recipient panel.
Separate source expansion, masking, sampling, response estimation and rendering
so each comparison answers one question. Measure render-seed variation early.

For biology, report distance, aligned cross, centered spread, uncentered energy
and retrieval together. Distinguish missing coverage from poor transfer on
covered targets. For emitters, check actual integer counts, moments, correction
budgets and DE outcomes; a proxy budget is not a gate.

Finally, write down what evidence would change the action before another run or
advisor round. Our most useful late result was a clear account of why appealing
local gains did not justify replacing the incumbent. We are publishing those
negative results so another team can test the assumptions earlier than we did.
