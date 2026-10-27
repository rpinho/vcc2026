**VERDICT:** Run a leave-one-family-out per-gene DE classifier that predicts sign and membership from target-line context, not just donor consensus. Use it to choose which genes get nonzero predicted response and in which direction.

**EXPERIMENT:** Train a gradient-boosted model (LightGBM) over (perturbation, gene) pairs to predict two things: sign(Δ), and whether the gene is in the true Wilcoxon DE set. Features:
- donor sign-agreement fraction and median donor effect
- the held-out line's control-cell correlation between the perturbed target and the gene
- the target's and the gene's baseline expression and detection rate in the held-out line
- the gene's program (PCA) loading in the held-out line's controls

Render the prediction as: gene included only if P(DE) exceeds a threshold, sign from the classifier, magnitude from the existing pooled/MLP amplitude.

**PROTOCOL:** Use 8-fold leave-one-family-out with no held-out-family responses in training. Tune the P(DE) threshold on inner folds only. Score every fold with the local official cell_eval2 vcc2026 profile against the incumbent, and bootstrap paired CIs over perturbations. CPU cost is roughly 6–12 hours.

**KILL TEST:** Kill it if any of these hold:
- direction_fidelity paired 95% CI lower bound ≤ 0
- score_avg fails to improve in at least 6 of 8 families
- an ablation without the control-correlation and baseline features retains more than 80% of the gain, meaning nothing beyond donor consensus was learned

**Justification (≈150 words):** You killed both pooling and scaling because direction, not amplitude, is wrong. Donors carry almost no per-gene direction signal for CD4 and HCT116. The only untapped signal that is available on CPU, with no downloads, is the held-out line's own control cells.

- **Why control covariance helps:** knocking down X in a given line moves the genes that co-vary with X there.
- **Why baseline expression helps:** genes that are unexpressed in that line cannot respond.

Predicting DE membership directly targets signature Jaccard (0.016), which is near floor. Choosing sign per gene with context targets direction_fidelity (0.177).

Gating to a confident subset also avoids the shrinkage trap, because it changes which genes respond rather than how much. It leaves PDS largely intact, since magnitudes reuse the existing pipeline. The ablation in the kill test stops you from crediting a disguised donor-consensus model.