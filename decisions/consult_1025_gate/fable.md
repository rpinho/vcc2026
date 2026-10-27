**VERDICT:** Replace the hard eligibility gate with a ranked fallback, and set the per-target budget from source-family data only. The gate is worth +0.43 measured. The budget is worth at most +0.09 even with the oracle. So fix the gate first, and treat the budget as a nested arm.

**EXPERIMENT:** Rewrite renderer eligibility from `confidence >= 1e6` to "rank all `keep & direction==bulkdir` genes by confidence (net first, then zs); call the top B_t." Run it end-to-end through the official cell_eval2 scorer on the CD4 outer fold.

**PROTOCOL:**
1. **Arm A (primary):** fixed B = 940, which is the existing source geometric mean.
2. **Arm B:** B_t = geometric mean of that target gene's n_calls across the 7 source families, shrunk toward 940 when it appears in fewer than 3 families. Tune the shrinkage by leave-one-source-family-out, never using CD4.
3. Signs and ranks use only source donor banks and frozen-MLP outputs. CD4 treated cells and the CD4 DE table are never read.
4. Log realized n_pred per target against B_t. This checks that the rendered cells actually reach p_adj < .05 under the Wilcoxon test.

**KILL TEST:** Kill if Arm A's official fidelity is below 0.5866 (the baseline anchor). Also kill if it falls more than 0.03 below the reimplementation's 0.6111. A third kill condition: precision below 0.35 on the 32 previously empty targets.

**Justification:**
The gate zeroes 42% of targets, and those targets hold half of all n_conf. That gate is the dominant loss. Precision where we already call genes (0.585) matches the baseline, so ranking past the gate should keep a sign quality close to that.

The real risk is rendering. The 0.6111 headroom assumes every budgeted gene becomes a called gene. Synthetic cells from low-confidence genes may not reach significance. That is why the official-scorer run plus the n_pred-vs-B_t log is the experiment, not more offline simulation.

Arm B attacks the budget without leakage. A perturbation's DE breadth is largely a property of the target gene, so the source families are a legitimate estimate of it.