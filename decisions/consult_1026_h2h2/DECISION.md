# DECISION — 2026-10-26 07:45 EDT (Fable, Astra, Kimi all consulted)

Per Ricardo's standing rule, the advisors decide; the orchestrator executes. This document records
what each said, where they agree, and what therefore runs. Nothing here is the orchestrator's own
scientific judgment.

## Unanimous agreements (all three)

1. **No GPU rental in the next 48 h.** All three said no, independently, with the same reasoning:
   the live work is 168x79 linear algebra and bootstraps; nothing identifies compute as the
   bottleneck; renting would "buy the feeling of progress" (Fable).
2. **Beating B0 alone is insufficient.** Every arm must also beat **B4 (training-mean)** and
   **B5 (shipped shrunk toward training-mean)**. A cosine win is not a result.
3. **The first invalidator is B5.** If B5 matches or beats a manifold arm, the apparent manifold win
   is shrinkage wearing an embedding, and STRING dies.
4. **Objective = improve the official score / rank, with expression MSE as the route.** Not "reach
   0.25" (no effect size supports it). Fable says "contest expression MSE, nothing else"; Kimi and
   Astra say expression MSE is the route but the *objective* is official score. Compatible.
5. **Stop DE work entirely**, stop cosine-only claims, stop re-litigating the three closed lines.
6. **B0's 0.4487 is not a bug hunt.** Astra and Kimi both corrected the orchestrator's "below chance"
   claim: fidelity is sign-precision x coverage, so its matched null is ~0.40, not 0.4929. Kimi adds
   that with 79 targets the standard error near chance is ~0.056, so the apparent deficit is ~0.8 SE
   — noise. Only cheap alignment/sign assertions, no inversion investigation.

## Where they disagree, and how it resolves

Fable **reversed its own earlier call** and now says: *if neither B2 nor B5 beats B0, run Astra's
response-reliability weighting Monday; Sidechain waits* — because sidechain staged progression is a
selection search over source subsets on a 79-target panel with no replication, and will overfit.

Kimi says: *if the manifold loses, run Sidechain staged source progression Monday; reliability
weighting waits* — because "reliability is not transferability."

Astra says: run **replicate-based reliability weighting of the shipped MLP's training loss** (not
another donor-average weighting, which already failed in `source_weight_1025`).

**Resolution by content, not by vote:** two of three (Fable, Astra) name reliability weighting; only
Kimi names Sidechain, and Fable explicitly reversed *to* reliability for stated overfitting reasons.
Kimi's objection is answered by Astra's refinement — Astra does not propose re-weighting donor
averages (which failed) but re-weighting the *training loss* using independent-replicate agreement.
That is a different mechanism and has not been tried. Sidechain's own audit (`SIDECHAIN_1023.md`)
already contains unsuccessful amplitude and limited source-addition tests.

## Therefore, mechanically (conditional on the h2h2 result)

**If h2h2 shows a manifold arm beats B0, B4 and B5 on official score with CI excluding zero:**
package-gate per §2 of Fable/Astra/Kimi — exact B0 replay, independent metric verification, gains
preserved after integer rendering, fixed fallback for the 2 STRING-missing targets (KDR, SOX2),
leave-one-family-out for >=2 non-H1 families, no calibration touching the 79, deterministic
regeneration. The deliverable is a go/no-go decision.

**If h2h2 shows no manifold arm clears B4/B5 (currently the likely outcome):**
1. Close STRING/manifold permanently in the report headline (all three agree).
2. Monday: **replicate-based reliability weighting of the shipped MLP training loss**, Astra's frozen
   protocol — per-source reliability from independent experimental replicate agreement on
   perturbation-minus-matched-control responses, training targets only; weights frozen before
   evaluation; matched weighted/unweighted retrains with identical architecture, seeds, budget;
   family-held-out evaluation; kill at <3% family-macro MSE improvement, <6/8 families, or any CI
   spanning zero. Astra's precondition: **no credible independent-replicate estimates => no run.**
3. Sidechain staged progression stays queued behind it, not deleted.

## Open item the advisors flagged
Fable's "shipped-shrunk-toward-training-mean wins => Monday is shrinkage calibration, nothing else."
The orchestrator's response-space diagnostics already show B5 (fixed 0.5/0.5) beating B0, and the
training mean alone beating B0. If B5's official score also beats B0, Fable's branch fires and the
Monday job becomes: fit lambda per family by leave-one-out on the 168, score once on the 79. That is
a cheap, high-value check the advisors have pre-authorized.
