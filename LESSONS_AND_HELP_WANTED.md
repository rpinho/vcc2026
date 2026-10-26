# What worked, what failed, and where help would matter

Our best official validation result is **0.16438** (billmA), about rank 463 of
1,580 when recorded. We are behind the strong entries. The point of this release
is to make the evidence and implementation easy to inspect, including failures.
The [README](README.md) separates official validation scores from blind final
submissions. We do not know our final score.

## Where the score comes from

We transfer measured responses from public source screens, condition on the new
context's control cells, and generate integer counts. AtlasShift's source transfer
and count generation gave us the first useful baseline. Spread calibration got
MSE above zero. Cross-source agreement and control-axis removal helped a little.
Confidence-ordered DE calls increased reach, but also made fidelity worse.

billmA adds our reconstruction of a publicly described response MLP: train only
on other source families, ensemble three seeds, and mix its direction into the
pooled response before count emission. It improves the official PDS component
from about 0.680 to 0.706. It is our implementation of a disclosed idea, not the
author's unavailable model or a reproduction of the author's score.

## Where we lose

**MSE is the largest gap.** Our useful later validation entries get roughly
0.06–0.095 scaled MSE, versus about **0.44** for the median of ranks 10–50 in
the earlier leaderboard snapshot. billmA gets 0.061. This means that recognizing
which response belongs to which target is not enough: the predicted response
still has the wrong size or geometry. The component gap alone is about 0.06
overall under equal averaging. Closing it is not a promised gain, and moving one
component can damage the others.

**More DE calls are not a free improvement.** v8/v9 roughly doubled reach relative
to v5b while fidelity and Jaccard suffered. Lower call budgets, an independent
mean-expression head, and PIE ordering did not establish a useful improvement.
The PIE H1 ordering check even had pretraining exposure and was only an upper
bound. We treat those tested DE levers as closed unless new evidence changes the
case. Some checks were uninformative because the harness missed the official
NMAE collapse; an inconclusive check is not a successful experiment.

**Our proxy ranked candidates in the wrong order.** The later audit of seven
historical recipes got 18 of 21 pairs right (Spearman 0.786), but reversed the
close v4/v5 comparisons. Four H1 layout/amplitude surrogates got only three of
six pairs right (Spearman zero). Exact local replays of v9b and slim_v5b were
missing, and the billmA sanity run used a different emission path. Running the
official scorer locally does not make the public task distribution match the
hidden contexts. We have no calibrated conversion from local gain to official gain.

**Selective scatter failed a real-cell DE check.** It looked attractive because
its estimated dispersion credit nearly doubled. On a real public H1 SEARCH task,
with billmA biology held fixed, changing scatter alone reduced local overall by
about 0.10 relative to the incumbent and sharply damaged fidelity and Jaccard.
Recalibrating amplitude still lost. A mechanically correct budget does not ensure
the resulting cells make useful DE calls. We kept the incumbent emitter.

**The slim experiment violated its intended spread budget.** It emitted about
6.5 times the noise credit, got zero official MSE and negative NMAE. A diagnostic
had been allowed to stand in for a required gate. New emission experiments must
fail closed on the actual emitted counts and check all six metrics.

## What we would try with more time or GPU

First rebuild faithful historical recipes on representative, untouched public
contexts and test whether their ordering transfers. Then expand source coverage
with matched controls, explicit missing-gene masks and honest source-family
holdouts. Compare more data at fixed model/emitter settings before changing the
architecture. A broader, family-balanced billmB has already failed the fixed local acceptance
check: a small H1 gain came with worse CD4 fidelity and reach. Source expansion,
training sampling and measured-gene masks changed together. We would isolate
those factors before another GPU run; repeating that combined recipe is not an
untested opportunity.

We would also separate response estimation from rendering: hold the mean/bulk
profiles fixed while checking actual integer-cell moments, DE signs and response
amplitude. Broad architecture searches and further DE sweeps are hard to justify
while the ranking harness fails. The current last safety entry stays the reference.

## Five questions where outside help would change the result

1. **Can a public evaluation panel rank our recent recipes faithfully?** Help build
   exact billmA/v9a/v9b/slim replays, independently hold out source families and
   controls, and explain which distribution or anchor differences reverse the order.
2. **How should response magnitude transfer to a context with only controls?**
   We need a measured comparison against the pooled source norm, using complete
   recipient-family exclusion and all six metrics rather than just a direction score.
3. **Which additional public screens add independent signal?** A reproducible
   adapter with batch-matched controls, native normalization, measured-gene masks
   and coverage accounting would be more useful than another model name.
4. **Can an emitter preserve useful dispersion without creating wrong DE calls?**
   Show actual count-level results with profiles held fixed, all six metrics and
   independent mean/bulk/depth checks. Selective scatter is the failed comparator.
5. **Can someone reproduce a stronger released 2026 recipe on the same public
   tasks?** Distinguish the score of the released code from the author's later
   leaderboard score, and state the license, omitted inputs and compute needed.

Please open an [issue](https://github.com/rpinho/vcc2026/issues) or send a PR.
Small reproducible comparisons, negative results and corrections are welcome.
Share only material you have rights to publish; public contributions do not
change registered team membership or submission rules.

## October 26 follow-up

The [new local diagnostics](RESEARCH_UPDATE_20261026.md) show an expression-MSE
improvement that loses the six-component composite. The fitted STRING/manifold
candidates failed the registered adoption rule on a disjoint-target H1 test.
The final safety entry remains the incumbent; local diagnostics do not supply
a new official score or an accepted replacement.
