Do not deliberate at length. Output the answer immediately, then at most 120 words of justification.

You are a scientific advisor for a computational-biology competition team. We are between submissions and building research toward a better model. The user's instruction: "make better models, try new ideas, research." Do not recommend closeout.

Facts:
- Best official validation 0.16438; leader 0.47203. Expression error dominates.
- Data on hand: 8 public Perturb-seq families, ~10,781 response rows, 300-target panel, response banks, frozen model source (5x-gene MLP -> residual on pooled donor direction, blended 0.3, then a count renderer). Donor coverage 298/300.
- Four untried mechanisms: (1) source-reliability weighting of donor pooling, (2) STRING/sequence/gene-text features, (3) missing-source prior, (4) saturation routing.
- Advisor A recommends (1) with held-out-family folds and rendered expression MSE. Advisor B (generic) recommends per-domain heteroscedastic weighting and distributional calibration.
- Budget: ~2 CPU-days. Must produce a defensible local signal without any official scorer.

Answer exactly four labelled lines and nothing else:
VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
