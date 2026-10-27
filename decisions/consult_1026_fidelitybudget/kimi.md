 **VERDICT:** The +0.053 projection is arithmetically correct but dangerously optimistic; it assumes the B=150 fidelity transfers perfectly to the leaderboard and ignores that **de_lfc_nmae is currently wrecking your composite** (your LFCs are 2× too small). The budget line is not dead, but you need **B≈75** (not 150) and **global LFC inflation** to maximize composite. Your decisive risk is realization failure, not collateral metric damage.

---

**1. Projection soundness & metric risks**
*Arithmetic:* Valid. (0.5892 − 0.4929)/0.301 ≈ 0.32 raw → /6 ≈ 0.053 composite gain vs. current. Fidelity is indeed 1/6 of the mean.

*Realization chokepoint:* The projection dies if realization < 30 %. If you only realize 15 calls instead of 150, fidelity collapses to ~0.25 and the gain evaporates.

*Other-metric risk:* **Minimal.** You are under-calling (1 vs. 30), so Jaccard and Reach will *improve* monotonically with B until B ≫ ref median. For Jaccard: at B=150, even with 20 % precision, Jaccard ≈ 0.13 vs. ~0 today. For NMAE: increasing B grows the intersection sample (denominator), but since your LFC errors are systematic (under-dispersion) not stochastic, the per-gene NMAE stays ~50 % (bad, but not worse). **No penalty** in Expression MSE or PDS cosine from changing DE thresholds alone.

**2. Optimal budget B**
**Pick B = 75.** Fidelity gains saturate (Δ+0.02 per 50 calls after B=100), while Jaccard starts decaying once B > n_conf (median 30) because the union term dominates. B=150 risks “discovering” false positives in the tail where your sign agreement drops; B=25 under-serves the long-tail targets with >100 DE genes. B=75 maximizes (Fidelity + Jaccard + Reach) – expect ~0.54 fidelity, ~0.35 Jaccard, vs. 0.167/0/0 today.

**3. Realization tactics (ranked by ROI)**
1. **Inflate predicted log2FC in rendered counts** (redistribute expected counts to your top-B predicted DE genes). *Compute:* trivial. *Failure:* Violates global library-size constraints; may trigger QC “implausible count distribution” flags if the scorer checks internal consistency (mean vs. variance).
2. **Shrink dispersion / liberalize the Wilcoxon** (use a smaller prior variance or raise α threshold to 0.10 for realisation). *Compute:* zero. *Failure:* P-values become uncalibrated; if the scorer uses a fixed pipeline, your submission is rejected or scores zero.
3. **Increase cells per target** (2–3×). *Compute:* linear. *Failure:* Square-root law means you need 4× cells to halve the standard error; expensive and hits diminishing returns.
4. **Increase sequencing depth.** *Compute:* linear. *Failure:* Technical variance is already negligible vs. biological variance; you gain almost no power per dollar.

**4. Reachable metrics without expression cosine**
*de_lfc_nmae:* **Yes.** It depends only on the intersection of DE calls, not the full expression space. You must **globally scale your predicted log2FC by 2.03×** (0.779/0.384) to match reference dispersion. This fixes the NMAE error (currently ~50 % → ~0 % on true DEs) without needing directional alignment. Do not use per-gene scaling (you proved it dead); use a single global multiplier.

*pds_cosine:* **Marginal.** If PDS is computed on the DE signature (mean profile of called DE genes), capturing the correct gene set (B-budget) lets you score >0 even if single-cell cosine is 0. However, it is coupled to LFC magnitude accuracy; combine with the global scaling above.

**5. Concrete next experiment**
**Experiment:** Submit a validation-fold prediction using **B = 75 fixed calls** (ranked by your internal p_adj) and **global LFC amplification factor 2.0×** applied to the rendered mean parameters (i.e., multiply predicted perturbed expression by 2.0 for the top-B genes, keep others unchanged, then render counts).

**Pass/Fail:** Compare the 4-metric DE sub-score (Fidelity + Jaccard + Reach + NMAE)/4.
* Pass if gain > +0.15 absolute (you moved from ~0.04 to >0.19).
* Fail if fidelity realization < 50 % (i.e., < 38 of the 75 calls survive the integer-count p_adj<0.05 filter), in which case the budget line is dead and you must abandon the post-hoc calling strategy entirely.