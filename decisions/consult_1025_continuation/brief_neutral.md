You are advising on a computational-biology benchmark for predicting gene-expression responses to gene knockdown from public cell-line data.

We have ~10,000 gene-perturbation measurements across 8 human cell lines. Features: gene embeddings, a protein-interaction graph, and measured pseudobulk responses. Current model maps gene features to a perturbation effect, averages compatible donor measurements, and predicts a held-out target's response. Best held-out validation score is 0.164 on a composite metric; most of the loss is expression error.

Three candidate improvements:
A. Weight donor cell lines by reliability (estimated noise/coverage) instead of equal weighting.
B. Fuse STRING graph embeddings with ESM2 protein embeddings and gene-text embeddings to predict responses for genes with no donor measurement.
C. Add a gating/routing head that chooses between the pooled donor direction and a learned residual.

Question: which is the strongest single experiment to run next on a CPU-only box in ~2 days? Give the hypothesis, the exact baseline, the leakage-safe split, the metric, and the result that would falsify it. Answer in <=230 words with four lines: VERDICT / EXPERIMENT / PROTOCOL / KILL TEST.
