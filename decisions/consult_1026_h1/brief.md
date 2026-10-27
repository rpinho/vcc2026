You are deciding the next experiment for a team competing in the Virtual Cell Challenge 2026.
Be brutal and specific. We want YOUR decision, not a menu of options. We must
launch something within the hour. We cannot submit (frozen safety entry stays last).

## Our position
Official best: billmA **0.16438**, rank ~463/1635. Leader **0.47203**.
Composite = equal-weight mean of six metrics. Official scaled components for billmA:
  expr_mse_unbiased_capped_norm 0.061  (raw 0.973; 1.0 = predicting the control, so 2.7% skill)
  pds_cosine                    0.706  (raw 0.8187)
  lfc_nmae                      0.116  (raw 0.9301)
  direction_fidelity           -0.066  (raw 0.4929)
  direction_reach               0.186  (raw 0.2438)
  sig_jaccard                  -0.016  (raw 0.0245)
Leader's per-component scaled: mse 0.793, pds 0.892, nmae 0.734, fid 0.564, reach 0.495, jac 0.059.

## What we measured since the last consultation (all real, in the metric's own space)
1. A just-finished four-arm official scoring run on the **H1 in-distribution fold** (H1 is a line we
   train on; 400 cells/target; split-half reliability 0.528):
     incumbent 0.021555 | A-shared 0.008394 | B-real 0.011010 | B-shuffled 0.009601
   The incumbent wins all pairwise comparisons. Recipient-control conditioning does NOT help, even
   in-distribution. That whole line is closed.
2. On that same H1 fold the incumbent's raw `expr_mse_unbiased_capped_norm` is **1.4097** — i.e. our
   predictions are 41% WORSE than predicting the control mean. This is not a cross-line artefact; it
   reproduces in-distribution. The model's expression output is actively harmful.
3. Predicted expression delta cosine with truth is ~0.005. Global amplitude, per-gene shape, cross-line
   transfer, global LFC inflation, and call-budget fixes are all measured dead.
4. Our ranking by predicted |lfc| is genuinely informative: top-10 precision 0.752, top-50 0.705.
   Sign accuracy rises with effect size (0.512 near zero to 0.636 above |log2FC|>0.8).
   Reference DE call counts are extremely skewed (median 31, p10 3, p90 565).
5. The realization mechanism works: we can emit a chosen call set and get ~89% of them significant.

## What we just launched (so you do not duplicate it)
A within-family response-manifold experiment: per family, 5-fold over targets, SVD on the training
targets' responses (k in 8/16/32/64 selected by inner CV), ridge from 64-d STRING target embeddings
to PC scores, reconstruct held-out responses. Pass = cosine >= 0.20 and MSE ratio vs do-nothing
<= 0.90 in >=4 of 7 families, with a shuffled-feature permutation control at or below do-nothing.
This is Kimi's prior proposal, restricted to folds where the premise holds.

## Answer these, decisively
1. Given that the expression head is worse than the control mean in-distribution (1.41x), is the
   manifold experiment the right next move, or should we instead simply **replace the expression head
   with the control mean** and spend all remaining effort on the four DE metrics? Which choice has the
   higher expected composite gain, and by how much in scaled points? Show the arithmetic.
2. If the manifold experiment fails its pass gate, what is the single next experiment? Name one.
3. Is there any *legal* use of same-cell-line replicate structure that could give expression cosine
   ~0.8 without leakage? If yes, state the exact construction and the leakage test that would falsify
   it. If no, say so plainly.
4. Given 7 days and this team's evidence, what is the realistic ceiling for our composite? Give a
   number and the reasoning. We would rather hear a hard truth than an encouraging guess.

FORMAT: five labelled lines first, in this exact order and nothing else on those lines:
VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
CEILING:
Then justification, max 220 words. Cite only numbers given above; do not invent measurements.
