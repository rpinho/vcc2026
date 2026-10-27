**VERDICT:** Build the **minimum viable test of the same-cell-line principal component (PC) hypothesis** (Question 3). This is the only potential path to the leaderboard top (0.47). Do not build the per-target budget package; it cannot close the gap to the leader.

**Pass/Fail Threshold:** Compute top-50 PCs from the held-out cell line's training perturbations, predict test responses via linear projection of these PCs (coefficients predicted from target gene embeddings), and measure expression cosine to truth. **Pass if cosine > 0.2; Fail if cosine < 0.1.** If it fails, stop all modelling immediately and document the negative result for the remaining 7 days. If it passes, Ricardo submits this model at 08:00.

---

## Reasoning

### 1. Per-Target Budget Rule
**Yes, it beats a fixed budget, but it is tactically useful and strategically irrelevant.**

A per-target budget proportional to predicted activity (e.g., $B_t \propto \sum_g \mathbb{I}(|\hat{\beta}_{g,t}| > \tau)$) exploits your genuine ranking signal (Top-10 precision 0.752) and addresses the structural skew where 52.6% of targets have $n_{conf} \leq 32$. Fixed budgets of 600+ destroy fidelity on low-activity targets (fidelity $\approx k/600$ when $n_{conf}=3$). A calibrated budget targeting $B_t \approx \hat{n}_{conf}$ could push fidelity from -0.066 to $\sim$0.3–0.4 and Jaccard to $\sim$0.2.

**However**, this only gains $\sim$0.10 on the composite (fidelity +0.4, Jaccard +0.2, others flat $\Rightarrow$ +0.10 mean). Your score would rise from 0.164 to $\sim$0.26, still 0.21 short of the leader. The shipped package already calls enough genes; the deficit is **precision**, which the budget fixes, but the larger deficit is **expression space accuracy**, which the budget cannot touch.

### 2. DE Metric Headroom Ranking
**Fidelity $\gg$ Jaccard $\gg$ Reach $\gg$ NMAE.**

- **Fidelity** (current -0.066, potential $\sim$0.6–0.75): You have Top-10 precision 0.752. Fidelity equals precision when $n_{pred} \geq n_{conf}$. By restricting calls to high-confidence predictions on low-$n_{conf}$ targets (52.6% of the data), you can flip fidelity positive and push it toward your ranking precision. Headroom: **+0.8**.
- **Jaccard** (current -0.016, potential $\sim$0.2): Killed by overcalling ($n_{pred}=700$ vs $n_{conf}=31$ gives max Jaccard 0.044). Fixing the budget lifts this. Headroom: **+0.22**.
- **Reach** (current 0.186, potential $\sim$0.4): Requires finding small-effect genes ($|\beta| < 0.2$) where your sign accuracy is 0.512–0.562 (near chance). Without expression direction, you cannot recover these. Headroom limited to large-effect genes: **+0.2**.
- **NMAE** (current 0.116, potential $\sim$0.2): Requires accurate log-fold-change magnitude. You have proven that global rescaling is dead (optimum is control) and per-gene magnitude correction fails even with oracle truth. Headroom: **Minimal**.

### 3. Exploiting Same-Cell-Line Replicates
**This is the only legitimate path to the top.** The leader achieves expression cosine $\sim$0.875 (replicate level). You have cosine 0.0045 (orthogonal). Cross-cell-line transfer is dead (Point 4), but **within-cell-line** covariance structure is untapped. The held-out cell line likely has training data (other perturbations) that define a low-dimensional response manifold.

**Minimum Viable Test (6h, 16GB CPU):**
1. Load all training expression responses for the held-out cell line (matrix $X \in \mathbb{R}^{N_{train} \times G}$).
2. Mean-center and compute truncated SVD (randomized, $k=50$) to get basis $U \in \mathbb{R}^{G \times 50}$.
3. Train a ridge regression to predict PC scores $z_t = X_t U$ from target gene features (e.g., embeddings, GO terms, or one-hot).
4. For test targets, predict $\hat{z}_t$, reconstruct $\hat{x}_t = \bar{x} + U \hat{z}_t$.
5. Compute cosine between $\hat{x}_t$ and your current predictions (or validation truth if available).

**Why this isn't leakage:** You are not using the test perturbations' expression; you are using the **training set** of the same cell line to learn the biological covariance structure (the "shape" of possible responses). This is analogous to using the empirical covariance for shrinkage, but data-adaptively.

### 4. What to Build Tonight
**Do not build the per-target budget package.** It yields at best 0.26, wasting the 08:00 slot on an incremental gain.

**Build the MV test described in Q3.** It requires no GPU and fits in 16GB: 20k genes $\times$ 5k samples = 400MB float32; randomized SVD completes in <30 minutes on a Mac Mini M-series CPU.

- **Pass (cosine > 0.2):** You have a signal. Scale this to full submission: use the PC reconstruction as your expression prediction, derive DE calls from it, and submit at 08:00. This is your only shot at 0.47.
- **Fail (cosine < 0.1):** Confirms that expression-space modelling is definitively closed for your architecture. Stop immediately. Do not spend the remaining 7 days tuning thresholds. Write the post-mortem: the signal is confined to DE ranking of large effects, and the competition requires expression replication quality that you cannot achieve.

**Brutal truth:** If the PC test fails, you are ranked 463 because you fundamentally cannot predict the expression vector (cosine $\sim$0), and no post-processing of DE calls fixes that. Stop.