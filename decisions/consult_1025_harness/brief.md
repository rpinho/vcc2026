# DECISION REQUEST #2 — same question, materially new evidence

You advised 3-0 on 2026-10-25 15:25Z: STOP AND SHIP. The orchestrator did not vote. Since then an
audit of the retained research corpus produced four new measured facts that were NOT in your first
brief. Re-answer with them. You may hold or reverse your verdict.

Answer in English. Be decisive. <=800 words. Structure:
  (A) VERDICT: [STOP AND SHIP | PURSUE, specify | HYBRID, specify]
  (B) WHAT TO DO, ranked, first action + hours
  (C) DID THE NEW FACTS CHANGE ANYTHING? If not, say so plainly and why.
  (D) RISK

## Unchanged constraints
Deadline Nov 05 23:59 UTC (11 days). Only the last final entry counts. Current last entry
`5NeeoWn96lfm23hgExcm`, SHA `ce78c81c...9dc4`, published FINAL, blind score unknown. Official
validation scoring is CLOSED (last public entry 2026-10-22; zero entries 10-23/24/25), so no new
candidate can ever be scored before the deadline. Local CPU only. No human in the loop overnight.
Agent forbidden to submit.

Best official validation score: billmA 0.16438 (rank 463). Recent attempts: v9b 0.12049,
slim_v5b 0.08833, main_a2 -0.008 vs safety (reverted). Leader 0.4716; ranks 10-50 median ~0.315.

## NEW FACT 1 — the harness failure is concentrated in one fold, not uniform
Same seven historical recipes, same scorer, two different folds:

| fold | pairs correct | Spearman |
|---|---|---|
| H1_rich_search (79 targets) | 15/21 (71%) | 0.536 |
| CD4_depth192_search | 18/21 (86%) | 0.786 |
| pooled 5-recipe macro | 9/10 (90%) | 0.900 |

H1_rich is the fold with only 8/79 = 10.1% target overlap with the official 300 panel. CD4 reaches
18/21. So most of the "inverts 6/15 orderings" claim comes from the least representative fold we hold.
This does not prove a better fold fixes it, and the pooled-macro 9/10 still fails on v4a/v4b.

## NEW FACT 2 — an unexplained ~2x gap between local and officially-implied direction correlation
Convert official scaled MSE through the measured anchors (baseline u=0.961342, replicate u=0.003359):

| entry | official MSE | implied raw u | implied rho |
|---|---:|---:|---:|
| billmA | 0.061 | 0.903 | 0.312 |
| v5b | 0.082 | 0.883 | 0.342 |
| v9a | 0.083 | 0.882 | 0.344 |
| v4a | 0.000 | 0.961 | 0.197 |

But the local H1 SEARCH replay measures per-target rho = 0.1535 and reports raw u = 0.9835, i.e. it
says our incumbent is WORSE than the official do-nothing baseline. Officially it is not (0.061 MSE
credit). So the local replay understates direction correlation by roughly 2x on the MSE component
(local 0.15 vs implied 0.31), while on PDS it overstates (local 0.754 vs official 0.680).

Interpretation is NOT settled. Candidates: the local fold's truth is genuinely harder; the local
emitter/baseline differs from the packaged one; the official val panel is easier/deeper; or the
anchor mapping does not transfer across folds. Nobody has measured which. This is the first concrete,
testable explanation of the blindness that does not require re-deriving anchors.

## NEW FACT 3 — the 13 ranked research ideas: 4 were never implemented, and they are exactly the rho-raising class
`research/IDEAS.md` (2026-10-16) ranks 13 experiments. Code search says these four have no
implementation anywhere in the private model directories:

- #6 bilinear / low-rank ridge target-to-response model using STRING + ESM + PIE features
- #7 select/combine sources by held-out RESPONSE reliability, not basal correlation
- #11 learned priors for genes with no source row (161 of the 300 official targets have no exact-target
  row in any of our 11 banks; union coverage 139/300 = 46.3%)
- #12 route between baseline and residual model by out-of-fold saturation/reliability

Everything actually run and closed — amplitude, layout, gate/deflation ablations, pool_k, alpha
re-solve, dense-bulk, eta0, learned-norm, family-quality — is an *emitter/amplitude/allocation* knob.
Under the algebra u_min = 1 - rho^2 those are all bounded by the rho we already have. The four
unimplemented ideas are the only proposed mechanisms that could RAISE rho, which the 2026-10-25
00:45 record itself names as "the only lever that moves MSE materially."

Caveat: #6/#7/#11 need features/data work and CPU training; none can be validated without the closed
official panel, and #11 is imputation (creates predictions, not evidence).

## NEW FACT 4 — the running record-only job self-aborted, so its result is void
`learnednorm` (substitute learned-response norm in the blend) died at 50/79 emissions. It hit its own
job-level 30 GB disk reserve: free fell 36.3 -> 26.1 GB as it wrote intermediates; run.py raises
`RuntimeError('30 GB free disk reserve')`. `dense` (-0.0111, refuted) and `eta0` (+0.00135) completed;
`learnednorm` has NO valid result and no exit code.

## Also available
- Direct same-line vs cross-line: 1.19x (rho 0.1435 vs 0.1205) on 94 official-overlap targets, with a
  shared input-control component (control mean alone gives 0.135) that may contaminate most of it.
- D=MOLM-13-like, E=U-87MG-like, F unresolved fibroblast-like; no same-line bank for any of them.
- 0 exact-matched local/official recipe pairs out of 67 retained records, so no transfer rule exists.

## The question
Do facts 1-4 change your verdict? Specifically:
(i) Is diagnosing the local-vs-official rho gap (Fact 2) worth doing as learning-only work, or is it the
    same "rebuild the harness" you told us to abandon?
(ii) Is any of the four unimplemented rho-raising ideas worth 11 CPU days with no validation yardstick,
    given that a failure to validate is guaranteed and only the last entry counts?
(iii) Or does STOP AND SHIP still hold unchanged?
