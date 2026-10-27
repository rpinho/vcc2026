We run a computational-biology competition (Virtual Cell Challenge). Predict single-cell gene expression
response to CRISPR perturbation in a held-out cell line. We are ranked 451/1635 on a hidden validation panel.
Our best official score is 0.16438; the leader is 0.47203.

The official composite averages six scaled components (mean of six z-like rescaled values):
  expr_mse_unbiased_capped_norm   (expression MSE; LOWER raw is better)   45.3% of our gap
  pds_cosine                      (perturbation-discrimination; higher)    17.6%
  de_wilcoxon_lfc_nmae            (log-fold-change error; lower)           14.4%
  de_wilcoxon_direction_reach_raw (fraction of true DE genes we call)      13.5%
  de_wilcoxon_direction_fidelity  (sign agreement of our calls)             7.4%
  de_wilcoxon_sig_jaccard         (overlap of significant sets)             1.9%

RAW VALUES, ours vs leader (validation panel):
  MSE       ours 0.9730   leader 0.2347
  PDS       ours 0.8187   leader 0.8919
  NMAE      ours 0.9301   leader 0.7339
  FIDELITY  ours 0.4929   leader 0.5640   (field median 0.5023 -- we are ABOVE median)
  REACH     ours 0.2438   leader 0.4946
  JACCARD   ours 0.0245   leader 0.0593

THREE FACTS ESTABLISHED TONIGHT (all verified against the live public leaderboard of 1635 entries, not a
local proxy):

FACT 1. The MSE component is CLIPPED AT ZERO in scaled space. Refitting on the uncensored region gives
base(raw)=1.0453, span=-1.1861, R2=0.835. Any submission with raw MSE above 1.0453 gets exactly 0.000 on
this component regardless of improvement. 61.7% of the board is parked at that floor. The no-skill point
(raw=1.0, "predict the control") maps to +0.0382, not 0. Our 0.973 -> +0.061. The leader's 0.2347 -> +0.793.
Matching the leader on MSE alone is worth +0.122 overall.

FACT 2. We tested global amplitude calibration (multiply the whole predicted response by one scalar ~0.35)
on two held-out cell lines, scored with the real official scorer. Raw MSE improved 13.7% (CD4) and 31.7%
(H1). But the official DE metrics are computed by a Wilcoxon test on the RENDERED INTEGER COUNTS, so shrinking
the response pushes genes below significance and the called-gene set collapses. Results:
   CD4: MSE gained +0.019 overall but direction fidelity fell 29% raw (=-0.081 overall) -> net -0.062
   H1:  MSE gained +0.043 overall and fidelity ROSE 2.1% (+0.006)                  -> net +0.049
The two folds disagree in sign. A single global scalar is a coin flip. Rejected as a global scalar.

FACT 3. Our fidelity (0.4929) is 87% of the leader's and ABOVE the field median. PDS (0.8187) is 92% of the
leader's. These are not our problem. Our problems are MSE (0.973 vs 0.235) and REACH (0.244 vs 0.495),
and to a lesser degree NMAE.

THE PROPOSED NEXT EXPERIMENT (we want your judgement, not agreement):
A "split-amplitude render". The insight is that MSE is driven by the MANY genes with wrong magnitude (all
~3000 genes enter a gene-averaged squared error), whereas the DE metrics only look at the FEW genes our
Wilcoxon test calls significant. A single global scalar cannot serve both -- which is exactly the tradeoff
we observed. So: shrink the amplitude of the NON-CALLED genes (buying MSE on the bulk of the matrix) while
leaving the CALLED/significant genes at full amplitude (preserving the Wilcoxon call set, hence reach and
fidelity). Both groups are rendered through the same existing integer-count emitter and v9 layout; only the
per-gene amplitude multiplier differs by group.

ALTERNATIVES WE ARE CONSIDERING:
  (a) split-amplitude render as above
  (b) abandon MSE, attack REACH (0.244 vs 0.495): call more DE genes per target, confidence-ordered
  (c) attack NMAE (0.930 vs 0.734): the log-fold-change magnitudes of called genes are too large
  (d) something else we have not thought of

WHAT WE ALREADY KNOW IS DEAD (do not re-propose):
  - learned low-rank donor->target response mapping (tested, all arms worse than incumbent; a shuffled-label
    control scored as well as the real map, i.e. capacity not signal)
  - per-target adaptive DE budget from source families (n_conf is NOT predictable from inference-time
    features: 10-feature ridge, log scale, leave-batch-out, R2 = -0.43, worse than the mean)
  - reliability-weighted pooling of source cell lines
  - global amplitude scalar (see FACT 2)

CONSTRAINTS: One resource-limited CPU machine, each official end-to-end fold evaluation takes
~1-2 hours. We have roughly 10 hours of unattended compute left before a human reviews. We can run maybe 4-6
such folds. There is no leaderboard feedback available to us (scores are hidden until submission, and we are
not submitting tonight) -- every decision must be made from held-out local folds with the official scorer.

QUESTIONS:
1. Which single experiment should we run next, and why? Be concrete about the mechanism you expect to pay.
2. For split-amplitude specifically: what shrinkage on the non-called genes, and how should it be fit without
   leaking the held-out line? What is the kill test?
3. Is there a reason to expect the called/non-called split to be stable, given that "called" is itself
   determined by the Wilcoxon test on the very counts we are modifying? (i.e. is this circular / does it
   just move the threshold problem?)
4. Rank (a)-(d) and say what you would do with 10 CPU-hours.

Answer concisely and decisively. Give a VERDICT, an EXPERIMENT, a PROTOCOL with numbers, and a KILL TEST.
Do not hedge across all options; pick one.
