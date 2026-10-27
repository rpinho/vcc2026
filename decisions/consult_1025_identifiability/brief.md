Output the four labelled lines FIRST. Justification after, max 150 words. Pick ONE option.

You are advising on a competition model. A concrete implementation blocker came up. Two cheaper mechanisms
are already dead on the official metric (source reweighting: +0.206% MSE, CI straddling zero; global
amplitude calibration: expression MSE -23% but direction_fidelity -29%, official composite -0.334). Per-gene
donor->target cosine is only 0.013-0.218, so per-gene DIRECTION is the disease. You (and two other advisors)
unanimously prescribed a learned low-rank donor->target map. The build then hit this:

THE BLOCKER. Advisor Fable specified fitting the map onto "the target line's own top perturbation-response
PCs". The held-out family is excluded from all treated responses (leakage rule). So the target line's
response basis cannot be estimated for the held-out family -> eligible own-line treated rows = 0 in all 8
folds -> the map is unidentifiable as literally specified. The job correctly refused to silently substitute
and stopped with n=0.

THE OPTIONS.
A. SHARED basis: fit ONE map and ONE response-PC basis on the seven training families only, apply it to the
   held-out family. Strict outer exclusion preserved. Loses line-specificity of the basis.
B. RECIPIENT CONTROL STATE: use the held-out family's UNPERTURBED CONTROL expression PCs (legally
   available) as the conditioning/recipient basis, and learn the map on training families. This is
   essentially Astra's cell-state-conditioned bilinear formulation, using control PCs x pooled-response PCs.
   Preserves strict exclusion AND keeps recipient-specificity, but the basis is control-state, not
   response-state.
C. Something else you specify.

Constraints: CPU only, 8 public families, 300-target panel, only CD4 and H1 have eligible frozen MLP
checkpoints (so at most 2 full-incumbent folds can be scored regardless). Official composite is dominated by
DE direction_fidelity and signature Jaccard. Do not propose another amplitude or pooling-weight change.

VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
