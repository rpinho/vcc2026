Output the four labelled lines FIRST. Justification after, max 170 words. Pick ONE experiment.

NEW EVIDENCE, measured tonight on the OFFICIAL scorer (not a proxy). This changes the diagnosis, so read it
before answering.

The official composite `score_avg` is the mean of six scaled components. On our real out-of-family CD4 fold,
ONE component accounts for 99.4% of that mean:
  de_wilcoxon_direction_fidelity_yield_raw  scaled -31.53   (99.4% share)
  de_wilcoxon_lfc_nmae                      scaled  -0.24   (0.8%)
  de_wilcoxon_sig_jaccard                   scaled  -0.05   (0.1%)
  pds_cosine                                scaled  -0.005  (~0%)
  expr_mse_unbiased_capped_norm             scaled   0.000  (~0%)
  de_wilcoxon_direction_reach_raw           scaled  +0.10   (-0.3%)
  score_avg = -5.2871.  Identical structure on all four CD4 arms we scored (99.3-99.7%).

Raw values vs the scorer's own anchors:
  ours:      fidelity 0.1767, reach 0.2796, Jaccard 0.0162, PDS 0.4974, exprMSE 1.2862
  baseline anchor (TRIVIAL): fidelity 0.5866, reach 0.2095, Jaccard 0.0263, PDS 0.5000, exprMSE 0.9705
  replicate anchor (good):   fidelity 0.5996, reach 0.8922, Jaccard 0.2397, PDS 0.9765, exprMSE 0.0265

So: our fidelity is 3.3x WORSE than a trivial baseline, and it dominates the score.

Two mechanisms already died tonight and this explains BOTH: amplitude shrink improved exprMSE 13-23% but cut
fidelity -29% and the composite worsened; a learned low-rank map nearly DOUBLED per-gene cosine
(0.0123->0.0239) yet fidelity fell slightly (0.1767->0.1699) and the composite fell with it. So per-gene
cosine is NOT what the score rewards.

RESOURCES: CPU only, 8 public Perturb-seq families, 300-target panel, existing response banks, frozen MLP,
working official scorer. Best official validation 0.16438; leader 0.47203.

QUESTION: what single experiment most directly raises `de_wilcoxon_direction_fidelity_yield_raw`?
Do not propose amplitude work, donor reweighting, or another low-rank linear map. And say explicitly whether
we should first spend time reading the metric's definition before modelling.

VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
