**Ruling, round 10**

**1. Is it worth the 9 days? Yes, but only as one gated bet with a hard stop.**

Overall is the plain mean of the six scaled metrics. I checked this on the v8a row: the six sum to 0.881, and 0.881/6 ≈ 0.147. v8b's PDS and MSE gains alone add (0.028 + 0.079)/6 ≈ **+0.018**. That is the only official datum above +0.01, and it is three times anything the call side has shown. The safety rail is label-free. Per `delta.py`, the gain keeps rising until B reaches the claim, so a cap of B ≤ 0.95·claim bounds the risk of overshooting.

**Main risk:** the two channels may not really separate.
- **Library size.** Moving counts into high-depth cells changes those cells' library sizes. That shifts every other gene's CPM, so NMAE can drift without any visible error. Per-cell library size must stay pinned, for example with zero-sum moves inside each cell.
- **Wilcoxon ranks.** Changing the shape of the distribution moves the ranks even when the mean is held, so the DE side moves too.
- **No local MSE check.** Local MSE is unreadable, so whatever ships rests on a projection.

**2. First experiment (B0, about 1.5 days, reusing the v7a dual-moment code)**

- **Target:** per perturbation, scale the pseudobulk delta by s ≈ √(0.95·claim/B). With B/claim ≈ 0.5, that is about ×1.38.
- **Method:** depth-tilted placement, with per-gene tilt γ solved so that the per-cell CPM arithmetic mean stays fixed, per-cell library size stays fixed, and the group sum hits its target.
- **Small groups:** leave the <10-cell stratum (12 targets) untouched. That is where the pooled-template tolerance failed on 10-22.

**Gate AMP-DM (all must pass):**
1. Per-(target, gene) CPM means: maximum relative error ≤ 1e-3. Local |ΔNMAE| ≤ 0.001.
2. Fidelity, reach and Jaccard each within the seed floor (±0.0028), with the R9 strata breakdown showing no stratum below −0.005.
3. Label-free B/claim ≥ 0.85 and ≤ 0.97.
4. Local PDS change ≥ 0.
5. Projected MSE = 0.079 × (B/claim achieved − 0.5)/0.5. Projected overall must be ≥ +0.010.

**Stop rule:** if the gate hasn't passed by **11-01 23:59 UTC**, B0 is final. If it passes, whether to submit is your call. If an official validation slot is still open, spend it on this arm, because that is the only real MSE readout.

**3.** Not needed. No other lever has an official datum above +0.01.
