# VCC consultation #2, 2026-10-25 ~00:15 ET. Read-only. Answer in <=1000 words, plain prose, no preamble.

## What just resolved since the last brief (measured, accept as given)
The CD4 replication finished: **REPLICATES**. Dense native layout beats v9 sparse on fidelity
(+0.096 H1 / +0.059 CD4 raw) and Jaccard (+0.039 / +0.023), loses reach (-0.051 both), all paired
intervals exclude zero. So the layout trade is real and generalizes across two recipients.

## Three new measurements that constrain what to do next
1. **Official instance of the same lever already lost.** v5b (dense) 0.1528 vs v8a (v9 sparse) 0.1470.
   Official component transfer vs local: fidelity 10x ATTENUATED, Jaccard 6x attenuated, reach 1.5x
   AMPLIFIED, MSE unmeasurable.
2. **Local MSE is structurally blind to layout.** Both arms scored exactly 0.9835370889129174
   (identical to 16 digits) because PRED_TRACE_CAP_K=1.0 saturates: v9 jk_ratio=2.036 vs native
   exactly 1.000, cap pins both to 1.0*real_tn. Any arm moving jk_ratio cannot be judged locally on MSE.
3. **The layout effect is entirely conditional on the hidden reference's call count.** Deltas by
   stratum (n_conf): `>=1000` -> fid +0.287 H1 / +0.395 CD4; `10-999` -> -0.0003 / +0.033;
   `<10` -> +0.010 but reach -0.157 H1. A conditional "dense only where reference is dense" policy
   would capture 98.5% of the fidelity gain for 2.8% of the reach cost. BUT n_conf is hidden at scoring
   time, and label-free proxies are weak: predicted magnitude Spearman 0.347; native per-target call
   count nearly constant (CD4 std 71 on mean 3073); v9's own count correlates 0.039.
4. **The official composite is the unweighted mean of six components** (verified on 7 entries, max
   residual 1.7e-4). Our best billmA 0.16438 decomposes as pds +0.1177 (71.6%), reach +0.0310,
   nmae +0.0193, mse +0.0102, jac -0.0027, fid -0.0110.
5. **Gap to the top-10 median (0.3322) by component contribution:** mse +0.0754, nmae +0.0277,
   fid +0.0243, pds +0.0208, reach +0.0193, jac +0.0060. MSE alone is 43% of the gap.

## Closed, do not propose
Amplitude/rho, pool_k, global bulk multipliers, B/claim, PC projection, learned-blend coefficient,
source breadth/fusion, per-target budget shrinkage, density sweeps, global dense layout for the
package, public MOLM-13/U-87MG acquisition (same-line transfer buys only 1.19x, and a shared
control component contaminates it).

## Facts
- Local CPU only, ~10 days, 2 validation submissions/day available, final panel D/E/F is BLIND and
  only the LAST entry counts. Safety entry 5NeeoWn96lfm23hgExcm is published with a hidden score.
- v5b zeroes 79% of the DE universe in the bulk profile (3,911 calls/target of 18,533) via
  `de_call_mask` with min_net>=1; 50 control PCs also projected out (`deflate`).
- Sidechain, the reference reproduction at 0.1657, has MSE 0.151 vs our 0.061 and identical PDS 0.678.

## Questions (answer each, concretely, in order)
1. Given MSE is 43% of the gap to top-10 and is the one component our local harness cannot measure,
   what is the single most valuable thing to do next that does NOT require a local MSE measurement?
   Be specific about the mechanism and what evidence would justify a visible-panel submission.
2. Is there a legitimate way to exploit the conditional-dense insight without knowing n_conf? Consider
   that the reference's call count may be estimable from the number of genes with a detectable effect
   in the target's own predicted response, or from cross-source agreement counts. Name the estimator
   you would test and the gate that would let it through.
3. What is your single highest-conviction concrete experiment for the next 8 hours of local CPU?
4. Name the specific way we are most likely to fool ourselves in the next 48 hours.
