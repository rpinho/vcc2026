VERDICT: Option B. Condition a low-rank bilinear donor→target map on the held-out family's unperturbed control-state PCs, with a response basis pooled across training families. Run Option A only as a nested ablation inside the same experiment.

EXPERIMENT: Predicted target response = U_pool · (W0 + Σₖ cₖ·Wₖ) · donor_response.
- U_pool: response PCs pooled from the 7 training families.
- cₖ: the held-out family's top-K control-expression PCs (K≈5–10).
- W0, Wₖ: rank r∈{4, 8, 16}, ridge-penalised, fit only on training families.
- Option A is the same model with all Wₖ set to 0.

PROTOCOL:
- Outer split: leave one family out, all 8 folds. Treated rows from the held-out family are never used for fitting.
- Inner split: leave one family out across the 7 training families, used to select r, K and the ridge strength.
- Score direction_fidelity and signature Jaccard on all 8 folds.
- Score the official composite only on the CD4 and H1 incumbent folds.
- Use paired bootstrap over targets for all comparisons.

KILL TEST: Kill Option B if either condition holds:
- It fails to beat donor-copy on direction_fidelity in at least 6 of 8 folds.
- Its official-composite gain on CD4 or H1 has a 95% confidence interval that includes zero.

If B does not beat its own A ablation, adopt A, because the control-state conditioning adds nothing.

**Justification:** Option B is the only option that is both legal and recipient-specific. The diagnosed failure is that each gene's response direction differs between lines, which is the problem a recipient-conditioned map is meant to fix. Option A discards exactly that information. Control-state PCs are legally available, and they plausibly encode why direction differs between lines, since baseline state shapes how a line responds. Nesting A inside B costs almost nothing. It also directly tests whether conditioning carries signal, so the fallback to A is decided by evidence rather than by assumption. Only two folds can score the official composite, so the six proxy folds supply power on the metric that matters most, direction_fidelity. The thresholds are set in advance to prevent a third post-hoc rescue.