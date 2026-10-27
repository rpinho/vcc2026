Output the four labelled lines FIRST. Justification after, max 180 words. Be concrete and pick ONE.

Setting: computational-biology competition (Virtual Cell Challenge). Research between submissions; a
separate frozen entry is safe, so nothing here risks a submission. We need the next MODEL experiment.

WHAT WE ALREADY KILLED tonight, with official-metric evidence, so do not recommend these again:
1. Reliability/source-weighted donor pooling. Rendered-MSE macro change +0.206%, paired 95% CI
   [-1.47%,+1.34%], 5/8 families. Dead.
2. Global amplitude/shrink calibration. Real out-of-family fold scored with the official cell_eval2 vcc2026
   profile: expression MSE improved 13.7%, but de_wilcoxon_direction_fidelity fell 29%, and the official
   composite score_avg went from -5.287 to -5.933 (delta -0.646, worse). A second independent
   pooling-only fold reproduced it (-5.047 -> -5.926, delta -0.878). Dead.

THE DIAGNOSIS (measured):
- We transfer perturbation response across cell lines by pooling donor screens. Per-gene cosine between
  pooled donor direction and the true held-out response is only 0.013 (CD4), 0.033 (HCT116), 0.086 (H1),
  0.112 (K562), 0.171 (Jurkat), 0.217 (HepG2), 0.218 (RPE1).
- A per-row ORACLE scale only recovers 17.6% MSE, and predicting zero already gets 8.2%. So amplitude is a
  symptom; per-gene DIRECTION is the disease.
- Official composite is dominated by DE metrics, especially direction_fidelity (incumbent 0.177) and
  signature Jaccard (0.016) — both near floor. PDS cosine is already decent (0.497-0.850).
- Best official validation 0.16438; leader 0.47203.

RESOURCES: CPU only, ~1-2 CPU-days. 8 public Perturb-seq families, 300-target panel, existing response
banks, a frozen 5x-gene MLP that predicts a residual on the pooled direction (blended 0.3), a working local
official scorer. No new large downloads. No GPU.

QUESTION: what single experiment most plausibly raises per-gene direction accuracy and DE fidelity?

VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
