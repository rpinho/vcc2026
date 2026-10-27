You are advising a team on the Virtual Cell Challenge (VCC 2026). Read and decide. Be blunt; we have
repeatedly fooled ourselves and need you to catch errors.

## Situation
Task: predict single-cell gene-expression response to CRISPR perturbation in a held-out cell line.
Our best official score 0.16438; leader 0.47203. Composite = EQUAL-WEIGHT MEAN OF SIX metrics
(each 1/6): expression MSE, PDS cosine, de direction fidelity, de direction reach, de sig jaccard,
de lfc NMAE.

## What we have already PROVEN DEAD (do not re-propose these)
Measured in the metric's own bulk_lognorm space on a CD4 fold:
- Global amplitude rescaling: 3rd independent failure. Every step toward zero improves MSE; the
  optimum is "predict the control". 0.50*ours=1.0604, 0.25*ours=1.0145, zero=1.0000.
- Per-gene / per-expression-bin magnitude ("shape") correction: DEAD EVEN WITH ORACLE knowledge of
  held-out truth (per-gene oracle MSE ratio 1.1719, 8-bin 1.2459, both worse than control 1.0000).
- Cross-cell-line transfer / "shared response axis": pure shrinkage, cos 0.00-0.03 with true mean.
- Our predicted expression delta has cos 0.0045 with truth (split-half reliability 0.317, so
  attenuation-corrected cos ~0.009). It is orthogonal, not noisy. A positive diagonal cannot rotate it.
So expression-space direction is essentially absent, and the leader's MSE 0.2325 implies cos~0.875
(i.e. replication quality) and is NOT reachable by modelling from what we have.

## NEW FINDING (this is what we want you to attack)
The failure is concentrated in the DE-call channel, where we DO have signal.
Deficit: our DE table emits median **1** significant call per target (p_adj<0.05); reference calls
median **30.5**. **42% of targets (32/77) call NOTHING**, which by the official convention
("silence is not rewarded": k=0, so score 0 when reference found anything) forces fidelity to 0.
Mean fidelity actual = 0.1676.

Cause is statistical POWER, not sign: sign agreement 0.5702 (better than chance); our log2FC sd is
0.384 vs reference 0.779 (~2x under-dispersed); fraction at p_adj<0.05 is 0.00095 vs 0.02985 (~30x
under-powered). Precision@n_conf on the few targets we do reach is ~0.73.

Leak-free fix tested: call a FIXED budget B per target, ranked by OUR OWN p_adj, keeping OUR OWN
signs. Fidelity scored exactly as official (k / max(n_pred, n_conf)):
  B=1(actual) 0.1676 | B=25 0.4643 | B=50 0.5297 | B=100 0.5700 | B=150 0.5892 | B=500 0.6112
Official fidelity anchors: baseline b~0.512, replicate r~0.813, span 0.301, weight 1/6.
Our official val fidelity is 0.4929 (scaled -0.063, i.e. slightly BELOW baseline).
Projected composite gain if val behaves like local fold at B=150: (0.5892-0.512)/0.301 = +0.256
scaled, x 1/6 = **+0.053 composite** (0.164 -> ~0.217).

## The decisive unknown
Fidelity is scored on the DE TABLE. To count on the leaderboard, calls must be REALISED in integer
counts: a gene counts only if p_adj<0.05 after rendering counts. A run measuring this "realisation
rate" is finishing now. If realisation is low, the budget gain does not transfer.

## Questions — answer these specifically
1. Is the +0.053 composite projection sound, or is there an error in the arithmetic/convention?
   Check especially: is fidelity really 1/6, is k/max(n_pred,n_conf) right, and does calling extra
   genes risk a PENALTY in any of the other five metrics (jaccard, nmae, reach)? Quantify the risk.
2. Does emitting a fixed B=150 significant calls risk hurting de_sig_jaccard or de_lfc_nmae enough to
   eat the fidelity gain? Which B would you pick to maximise the COMPOSITE, not fidelity alone?
3. How should we realise more calls in integer counts: increase cells per target, increase sequencing
   depth, inflate predicted log2FC, or shrink the dispersion used in the Wilcoxon test? Rank these by
   expected realisation per unit compute, and name the failure mode of each.
4. Is there anything in the six-metric structure we are still ignoring that is reachable WITHOUT
   expression-space cosine? Specifically comment on de_lfc_nmae and pds_cosine.
5. Give ONE concrete next experiment, with a clear pass/fail threshold.

Answer with a VERDICT first, then the reasoning. Be concise but specific. If you think the whole
budget line is a dead end, say so plainly and say why.
