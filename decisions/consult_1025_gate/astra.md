VERDICT: Test a **gate-free renderer with a source-trained, per-target call budget**; leave the frozen MLP and effect amplitudes unchanged.
EXPERIMENT: One preregistered renderer ablation: current renderer versus gate-free fixed-B=940 versus gate-free adaptive-B. The fixed-budget arm isolates the established eligibility bottleneck.
PROTOCOL: Preserve structural gene exclusions and frozen predicted signs. Remove donor-significance and donor-consensus eligibility vetoes; use donor evidence only as a soft ranking signal. Select the top B genes. Learn B exclusively from non-held-out families using nested family-held-out predictions and source DE-count labels.
KILL TEST: On untouched outer families, require adaptive-B to beat both the current renderer and gate-free B=940 on mean official fidelity. If it cannot beat B=940, reject the budget estimator—not gate removal. Verify that reduced zero-call frequency and increased *actual* n_pred accompany the gain.

**Justification**

Fit a shallow Poisson count tree to source-family n_conf, computed with the scorer’s exact exclusions and sampling protocol. Inputs available at deployment: frozen predicted-effect magnitude quantiles, control-only variance summaries, and continuous donor-support summaries. No held-out treated cells, DE tables, or target-specific reference counts enter features or fitting.

Generate training features through inner leave-family-out runs, so donor banks and any fitted preprocessing exclude the feature target’s family. Select tree depth and a budget multiplier by inner-fold **official fidelity**, not count-prediction error. Refit on permitted families; lock everything before outer scoring.

Retain the renderer’s continuous priority ordering, but remove its significance/consensus vetoes. Cap B only by the structurally allowed gene pool. Changing eligibility and B must not alter predicted effect amplitudes.

Crucially, B is a requested budget, not necessarily realized n_pred: the official Wilcoxon test still decides significance. Audit B, n_pred, zero-call rate, and fidelity together. The reported 0.6111 is compelling headroom, not a guaranteed executable-renderer result.