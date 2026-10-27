We run a computational-biology competition (Virtual Cell Challenge). Task: predict single-cell gene expression
response to CRISPR perturbation in a held-out cell line. We are ranked 451/1635 on a hidden validation panel.
Our best official score is 0.16438; the leader is 0.47203.

The official composite = mean of six rescaled components:
  expr_mse_unbiased_capped_norm   (expression MSE; LOWER raw better)          45.3% of our gap
  pds_cosine                      (perturbation discrimination; higher)       17.6%
  de_wilcoxon_lfc_nmae            (log-fold-change error; lower)              14.4%
  de_wilcoxon_direction_reach_raw (fraction of true DE genes called)          13.5%
  de_wilcoxon_direction_fidelity  (sign agreement of our calls)                7.4%
  de_wilcoxon_sig_jaccard         (overlap of significant sets)                1.9%

RAW VALUES on the hidden validation panel, ours vs leader:
  MSE       ours 0.9730   leader 0.2347
  PDS       ours 0.8187   leader 0.8919
  NMAE      ours 0.9301   leader 0.7339
  FIDELITY  ours 0.4929   leader 0.5640   (field median 0.5023 -- we are ABOVE median)
  REACH     ours 0.2438   leader 0.4946
  JACCARD   ours 0.0245   leader 0.0593

== WHAT WE ESTABLISHED TONIGHT (all measured, much of it new) ==

1. THE MSE COMPONENT IS CLIPPED AT ZERO in scored space. Uncensored refit: base(raw)=1.0453,
   span=-1.1861, R2=0.835. Any entry with raw MSE > 1.0453 scores EXACTLY 0.000 regardless of improvement;
   61.7% of the 1635-entry board is parked there. The metric's no-skill point raw=1.0 ("predict the
   control") maps to +0.0382 scaled. Ours 0.973 -> +0.061. Leader 0.2347 -> +0.793. Perfect MSE = +0.878.

2. WE MEASURED THE CONTROL REFERENCE. On our held-out CD4 fold, a submission that simply predicts the
   INPUT-CONTROL MEAN scores raw MSE 1.1087. Our incumbent scores 1.2862 (i.e. 16% WORSE than doing nothing:
   negative skill). A previously tested amplitude-calibrated arm scores 1.1104 -- it removes the excess
   error but lands exactly ON the control reference, gaining NO skill (skill -0.15%).

   => On that fold our model adds nothing over "predict the control" in expression space, and the
      amplitude fix merely cancels the model's own added noise.

   BUT on the hidden val panel our submission is 0.973, i.e. BETTER than the control reference (1.0) by
   2.7%. The local fold and the real panel disagree about whether we have positive skill. We do not know
   which regime the hidden panel will resemble next.

3. DIRECTION FIDELITY IS NOT OUR PROBLEM. Ours 0.4929 is 87% of the leader's and ABOVE the field median.
   PDS 0.8187 is 92% of the leader's. Both are close to competitive. The large gaps are MSE (45.3% of the
   deficit) and REACH (13.5%), plus NMAE (14.4%).

4. THE DE METRICS ARE COMPUTED ON RENDERED INTEGER COUNTS. Our pipeline rounds the continuous prediction to
   integer UMI counts under exact row/column total constraints (v9 layout). A Wilcoxon test on those counts
   decides which genes count as "called". So shrinking predicted amplitude pushes genes below significance
   and the call set collapses -- reach and fidelity fall with it. This is why a single global amplitude
   scalar is a coin flip: on one fold it lost 29% of fidelity (net -0.062 overall), on another it gained
   2.1% (net +0.049).

5. ALREADY TESTED AND DEAD (do not re-propose these):
   - learned low-rank donor->target response mapping: all arms worse than incumbent; a label-shuffled
     control scored as well as the real map, so it was fitting capacity, not signal.
   - per-target adaptive DE budget predicted from source-family reference DE counts: n_conf is NOT
     predictable from inference-time features (10-feature ridge, log scale, leave-batch-out, R2 = -0.43,
     worse than predicting the mean).
   - reliability-weighted pooling across source cell lines.
   - global amplitude scalar (item 4).
   - "call more genes by widening the DE eligibility gate": partially run; it realises calls but the
     evidence so far is that added calls are low-precision.

== CONSTRAINTS ==
One resource-limited CPU machine. One official end-to-end held-out fold evaluation
(render ~300 targets to integer counts + score) costs ~1-2 hours. We have ~8-10 hours of unattended compute
before a human reviews. That is 4-6 fold evaluations, total. There is NO leaderboard feedback available
(scores hidden until submission; we are not submitting tonight). Every decision must come from held-out
local folds scored with the official scorer.

== THE QUESTION WE ACTUALLY NEED ANSWERED ==
Given that we are at 2.7% expression skill on the real panel but NEGATIVE skill locally, and given that
fidelity/PDS are already near-competitive while MSE and REACH are far off:

1. What is the single highest-value experiment to run in the next 8 hours? Be concrete about the mechanism
   by which it is supposed to move the score, and state what the expected effect size is.
2. Is "get better than the control mean in expression space" even the right target, or is that a trap given
   the clipping and given that the leader is at 0.235 (which requires real per-gene structure, not shrinkage)?
3. We suspect there may be a cheaper structural win: our predicted response may have the right DIRECTION but
   wrong per-gene MAGNITUDE PROFILE (e.g. too flat, or dominated by highly-expressed genes). What diagnostic
   would distinguish "wrong amplitude" from "wrong per-gene shape", and which of those would you fix first?
4. Rank these and pick ONE: (a) per-gene magnitude/shape recalibration; (b) attack REACH by calling more
   genes with better precision; (c) attack NMAE by shrinking log-fold-changes of called genes; (d) something
   we have not considered.

Answer concisely and decisively: VERDICT, EXPERIMENT, PROTOCOL with concrete numbers, KILL TEST.
Pick one. Do not hedge across options.
