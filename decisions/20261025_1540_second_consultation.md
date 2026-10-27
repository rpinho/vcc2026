# 2026-10-25 15:40Z — second three-model consultation: verdict holds, bounded diagnostic authorized

Ricardo asked whether the local harness could be improved and whether the research corpus had been
audited. The audit surfaced four facts that were not in the first brief. The same three advisors were
re-asked, with a neutral brief, and no vote from the orchestrator.

Records: `decisions/consult_1025_harness/{brief.md,brief_short.md,astra.md,kimi.md,fable.md}`.

| advisor | verdict on the submission | on diagnostic work |
|---|---|---|
| Astra | **freeze, do not replace** | HYBRID: allow <=16 h learning-only discrepancy audit |
| Kimi | **freeze, do not replace** | HYBRID: allow <=4 h read-only Fact-2 diagnostic |
| Fable | **keep finalized version unchanged** | bounded diagnostics only; "no new version before the deadline, whatever the diagnostics show" |

**3-0 on the submission, unchanged from 15:25Z.** No model authorized a replacement.

## What changed in their reasoning

All three withdrew the stronger claim that the research space was exhausted:

- Astra: "'all useful mechanisms have been exhausted' is no longer defensible. It does not establish a
  sufficiently reliable way to select a new final entry."
- Kimi: "The rho-raising search space was **not** exhausted; any earlier implication otherwise should be
  withdrawn. But 'unimplemented and mechanistically relevant' is not evidence of improvement."
- Fable: "'never tried' is zero evidence of payoff."

All three also rejected spending the remaining 11 days on ideas #6/#7/#11/#12, because the closed
scoring server dominates: an unselectable artifact has no value at the deadline.

Astra and Kimi both flagged that Fact 2 must be *audited before it is believed*. Kimi: "Treat
officially implied rho only as an aggregate sanity check, not per-target truth." Astra: "raw loss does
not uniquely identify correlation without the amplitude and normalization assumptions." That is
correct and the job below tests it first.

## Corrections to the record

1. `learnednorm` is **VOID, not negative** — it aborted at 50/79 on its own 30 GB disk reserve
   (`run.py:49-50`; free fell 36.3 -> 26.1 GB). No result exists. Do not cite it either way.
2. The blanket "harness inverts 6/15 orderings" was stated too uniformly. Per fold: H1_rich 15/21,
   CD4_depth192 18/21, pooled 5-recipe macro 9/10. H1_rich is the 10.1%-overlap fold.
3. The 4.8x same-line figure is a *same-experiment* figure (H1), not purchasable for D/E/F. The
   measured same-line premium is 1.19x and partly contaminated by a shared control component.

## Authorized work

One bounded, learning-only discrepancy audit (`final_push/harness_reconcile_1025/`), hard-capped and
forbidden from touching the package. It exists to make the next cycle's instrument honest, not to
produce a candidate. No submission.
