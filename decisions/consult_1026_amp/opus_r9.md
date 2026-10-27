**Round 9 advice: run the detectability arm first; any submission needs a pre-registered bar that touches DE only**

**1. Experiments, in this order**

1. **`DET75_885`** (Skadi-style call selection). Keep the 885 calls per target and the existing set of effect sizes and signs. Change only which genes get the calls: rank genes 0.75 by expected reference Wilcoxon z (predicted lfc × control CPM × n_real) and 0.25 by predicted rank. This goes first because the sign audit shows precision is 0.677 on genes the reference also calls and 0.548 on genes it doesn't. That gap is the only lever that fixes signs instead of rearranging them.
2. **`RH20`** (reach head). Keep the call set exactly as it is. Positions 1–20 go to genes that are sign-certain and likely to be reference-significant. Below the threshold, p-values stay graded. By construction only reach can move, which makes this the cleanest arm to read.
3. **`B250`** (flat budget of 250). This is a bracket, not a candidate. The precision decay is shallow (0.563 to 0.520), so trimming mostly buys jaccard and risks coverage on mid targets with n_real between 250 and 999. The official v5b-vs-v9a comparison favoured density on fidelity, so this one is the weakest prior.

**2. Decision instrument: yes, with conditions**

- **Allowed:** reading the 53 weak+mid targets plus a composite reweighted to 31.5/65/3.5.
- **Not allowed:** calib3's global ratios (fidelity 9–12×, jaccard 7–9×). They come from one strong-dominated layout comparison.
- **Instead:** take raw deltas per stratum, convert them with the official slopes (3.30 fidelity, 2.88 jaccard), and treat local reach as a conservative estimate.
- **Pairing:** compare each arm against B0 with the same seeds.
- **Replication:** repeat on the CD4 recipient, whose stratum mix is 29/63/8.
- **Caution:** local weak+mid precision is 0.55 while official is about 0.49. The level does not transfer, so only trust deltas.

**3. Can this become a recommendation?** Yes, but only for a change that touches DE alone, and only if all of the following hold (fix this bar by 10-28):
- Pseudobulk sums and MSE/PDS/NMAE match B0 to within 1e-9, which bounds the downside to the DE side.
- The reweighted gain is at least +0.005 overall. That is twice the seed floor of about +0.002 to +0.003.
- The arm wins on at least 60% of weak+mid targets (sign test p < 0.05).
- CD4 gives the same sign.
- The strong targets show no loss larger than their 3.5% weight can absorb.

If an arm passes, give Ricardo "B0 + `DET75_885`" with an honest probability, not a predicted size. If none passes, B0 stays final. That would be a measured result, not giving up.
