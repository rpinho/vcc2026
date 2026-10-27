Do NOT deliberate at length. Output the four labelled lines FIRST, then at most 150 words.

Setting: computational-biology competition, research between submissions. We are picking the next model
experiment. Prior advisor round recommended reliability-weighted donor pooling; we ran it as a
preregistered, independently verified experiment and it FAILED (macro MSE +0.206%, paired 95% CI
[-1.47%, +1.34%], 5/8 families). It is dead and we will not revisit it.

New measured diagnostic (honest, no oracle):
- Setup: 8 public Perturb-seq families, 18,946 rows, 18,533 genes. For each held-out family, donor
  response is pooled from the OTHER families; prediction = pooled response; score = log-expression MSE
  on the measured-gene mask; the target's own gene is always masked out.
- Baseline uniform pooling, family-macro MSE 0.015061.
- Fitting ONE global nonnegative scalar on the seven training families, then applying it to the held-out
  family, gives macro MSE 0.013189: -12.43% relative, 7/8 families improve. The fitted scale is extremely
  stable across folds: 0.2547, 0.2715, 0.3140, 0.3150, 0.3211, 0.3240, 0.3448, 0.3568.
- Two reference points: per-row ORACLE scale reaches only -17.57% (so a single scalar captures most of the
  available amplitude gain), while predicting ZERO response reaches -8.20%.
- Per-gene cosine between pooled and true response is low: 0.013 (CD4), 0.033 (HCT116), 0.086 (H1),
  0.112 (K562), 0.171 (Jurkat), 0.217 (HepG2), 0.218 (RPE1). Mean |pred|/|true| = 0.70.

Our incumbent model blends a learned MLP residual at 0.3 into the pooled direction and RESCALES the
result to the pooled direction's norm, so it inherits this amplitude. Best official validation 0.16438;
leader 0.47203. The official metric is a debiased, capped expression MSE plus DE metrics (PDS cosine,
signature Jaccard, direction fidelity/reach), so a pure shrink could trade expression gain for DE loss.

Question: what is the single highest-value next experiment, and should it be this amplitude calibration?

Answer exactly four labelled lines, then optional justification:
VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
