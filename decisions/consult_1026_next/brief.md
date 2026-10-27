You are deciding the next experiment for a team competing in the Virtual Cell Challenge (VCC 2026).
We have badly fooled ourselves several times tonight. Be brutal. We want your decision, not options.

## Setup
Predict single-cell gene expression response to CRISPR perturbation, held-out cell line.
Our best official score **0.16438**, leader **0.47203**, rank ~463/1635.
Composite = EQUAL-WEIGHT MEAN OF SIX metrics: expression MSE, PDS cosine, and four DE metrics
(direction fidelity, direction reach, sig jaccard, lfc NMAE). Official anchors per metric: a
baseline b (context mean response) and a replicate r (reference vs itself); scaled s=(u-b)/(r-b).
Expression MSE is clamped [0,1]; NMAE floored at -6.

## PROVEN DEAD tonight, measured in the metric's own space (do NOT re-propose)
1. **Expression direction is absent.** Our predicted expression delta has cosine **0.0045** with
   truth (split-half reliability 0.317; attenuation-corrected ~0.009). Orthogonal, not noisy.
   A positive diagonal cannot rotate a zero vector.
2. **Global amplitude rescaling.** 3rd independent failure. Our=1.2459, 0.50x=1.0604, 0.25x=1.0145,
   zero=1.0000 (MSE ratio). Optimum is "predict the control".
3. **Per-gene / per-bin magnitude ("shape") correction.** Dead EVEN WITH ORACLE truth: per-gene
   oracle 1.1719, 8-bin oracle 1.2459, both worse than control 1.0000.
4. **Cross-cell-line transfer / shared axis.** Pure shrinkage (0.50x alone = 1.0604; axis = 1.0617).
5. **Global LFC inflation x2.03.** Measured NMAE 0.9171 vs null predictor 1.0000. No gain.
6. **Call budget.** We thought our renderer undercalled (median 1 significant call/target locally vs
   reference 30) and that a fixed budget B would fix fidelity 0.168 -> 0.59. **Baseline invalid:**
   the shipped package's official fidelity is ~0.492 raw (scaled -0.066), which is impossible if it
   called ~1 gene/target. Cross-check: v8a (known 600-700 calls/group) has fidelity -0.064, reach
   0.162; billmA has -0.066, reach 0.186. The shipped package already calls plenty.

## WHAT IS ACTUALLY TRUE AND NEW (measured tonight)
(a) **The signal is monotone in effect size.** Sign accuracy vs reference |log2FC|:
    0-0.05: 0.512 (chance) | 0.05-0.1: 0.534 | 0.1-0.2: 0.562 | 0.2-0.4: 0.599 |
    0.4-0.8: 0.628 | >0.8: 0.636.   On true DE genes: 0.5644 (n=19,719).
(b) **Precision of the top-N genes ranked by PREDICTED |lfc|**: N=10 -> **0.752**, 25 -> 0.711,
    50 -> 0.705, 100 -> 0.692, 200 -> 0.679, 885 -> 0.633. So our *ranking by predicted magnitude*
    is genuinely informative, much more than our significance calls.
(c) **Reference call counts are extremely skewed**: n_conf p10=3, p25=10, median 31, p75=153,
    p90=565. **52.6% of targets have n_conf <= 32.** A single fixed budget is structurally wrong.
(d) **Realisation mechanism works**: rendering a chosen call set to integer counts realises 89.17%
    of calls at p_adj<0.05 (62,342/69,915). So we CAN control the emitted call set.
(e) Our official component scores: pds **0.706**, reach 0.186, nmae 0.116, mse 0.061,
    fidelity -0.066, jaccard -0.016. Leader: pds 0.892, mse implied by raw 0.2325.
(f) Leader's MSE 0.2325 implies expression cosine ~0.875, i.e. **replication quality**; the
   cross-replicate ceiling we measured on the same line is 0.81-0.87.

## The strategic question
Given (1)-(6), expression-space modelling is closed. Our only real signal is on **large-effect DE
genes**, where sign accuracy reaches 0.64 and top-10 precision 0.75. Meanwhile 52.6% of targets have
<=32 reference DE genes, and fidelity = k/max(n_pred,n_conf) is *precision* once n_pred>=n_conf.

We have ~6 hours of compute, one 16GB CPU Mac Mini, and cannot submit (frozen entry stays).

## Answer these
1. Is there a **per-target budget rule** (predicted from the model's own outputs only, e.g. predicted
   response magnitude / number of predicted-significant genes) that beats a fixed budget, and does it
   survive the fact that the shipped package apparently already calls enough? Or is (d)+(c) useless
   because the real deficit is precision, not count?
2. Fidelity = precision in our regime, and we cap at ~0.71 (B=50) vs leader 0.78. Reach we already
   have 0.186 vs leader 0.892. **Which of the four DE metrics has the most real headroom for us**,
   given the numbers in (a)-(e)? Rank them and justify with the actual numbers.
3. Is there a legitimate way to exploit **same-cell-line replicate structure** (the leader appears to
   reproduce measured responses) that is not leakage and that we could build in 6h on CPU? If you
   believe this is the only path to the top of the board, say so and describe the minimum viable test.
4. Given we cannot submit, what is the single most valuable thing to BUILD tonight: a better model, a
   verified package for Ricardo to submit at 08:00, or a decisive negative result that stops us
   wasting the remaining 7 days? Pick ONE and give a pass/fail threshold.

VERDICT first, then reasoning. If you think we should stop modelling and write up, say that plainly.
