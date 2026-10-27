# VERDICT

**Build one thing tonight: a faithful local reproduction of the six official metrics, then run ONE experiment through it — n_conf-matched call budgets (oracle, then predicted). Gate: if the oracle-n_conf budget does not lift the local composite by ≥ +0.04 over the shipped package under the faithful scorer, stop modelling at 08:00 and write up.** Do not hand Ricardo a new package unless that gate passes. Do not train anything.

You are not at the point where another model helps. You are at the point where you cannot measure. Exhibit A: you "proved" the renderer called ~1 gene/target and the official score says it calls hundreds. Every number you have measured locally tonight was produced by a pipeline that failed that check. Until the local scorer reproduces, on the shipped package, (i) call counts in the hundreds and (ii) raw fidelity ≈ 0.49, nothing you measure is interpretable, and (b) above is as suspect as the call budget was.

Three more things you have not yet caught yourselves on:

1. **The precision numbers in (b) are pooled across targets and therefore dominated by the p90 targets (n_conf ≥ 565).** For a median target (n_conf = 31), top-885 precision is mathematically ≤ 3.5%. A pooled 0.633 at N=885 is only possible if the big-n targets carry it. So "top-10 precision 0.752" may be true for the 10% of targets where a highly-expressed gene is DE almost regardless, and false for the 52.6% you are trying to fix. Stratify by n_conf quartile before believing it.
2. **You have not subtracted the baseline.** Scaled fidelity is relative to the context-mean response, which raw-scores ≈ 0.5 — it already calls the generic ribosomal/stress genes. If your top-10 is the baseline's top-10, scaled gain is zero no matter how good 0.752 looks. Compute precision of top-N *after removing genes the context-mean baseline would call*.
3. **"Expression direction is absent" (cosine 0.0045) and PDS 0.706 cannot both be true in the same space,** and 0.0045 is also inconsistent with 0.636 sign accuracy on |LFC|>0.8 genes (big genes dominate a cosine). Almost certainly the 0.0045 was computed in count space, where a handful of high-count genes we predict as noise dominate; in log space direction exists. That means expression is closed *for MSE* — correct — but do not let that measurement be used to argue the DE layer has no direction. It does, weakly, and it's the only thing you have.

---

## 1. Per-target budget rule

**In principle yes, and the deficit is count-matching more than precision — but you have not measured the one number that decides it.**

With fidelity = k/max(n_pred, n_conf) and precision p(N) declining slowly (0.75 → 0.63 across two orders of magnitude of N), for a target with n_conf = n:
- N < n: fidelity = N·p(N)/n, increasing in N
- N > n: fidelity = p(N), decreasing in N
- Optimum is **N = n_conf exactly**. Same for Jaccard = k/(N + n − k), which peaks at N ≈ n with value p/(2−p) ≈ 0.4 at p=0.6.

n_conf spans 3 → 565 (p10 → p90). Any fixed B is wrong by ≥ 10× for at least a quarter of targets; for n_conf = 3 and B = 50, fidelity ≤ 0.06 even with perfect precision. The shipped package calling hundreds per target is why fidelity/jaccard sit at ~0 scaled: for 75% of targets n_pred ≫ n_conf and you pay precision at N ≈ 600. So (d)+(c) is not useless — it is the only lever on two of the six metrics — **conditional on**:

**The decisive measurement (20 minutes, run first):** on the local reference line, Spearman(predicted-response magnitude or predicted n_sig, n_conf). Also compute the *oracle* composite with N = n_conf per target. If oracle gives < +0.04 composite, the lever is dead regardless of predictability. If oracle passes and ρ ≥ 0.5, use the predicted rank-quantile mapped onto the empirical n_conf distribution (p10=3 … p90=565). If ρ < 0.3, fall back to a fixed budget near the Jaccard-optimal point, which is near the *median* (~30), not 50 — and expect the gain to be small.

Also account for the reach cost: reach wants N large, fidelity/jaccard want N = n_conf, equal weights. Shrinking from ~600 to ~n_conf will drop reach below 0.186. Net is probably positive but it has to be measured under the faithful scorer, not assumed.

## 2. Headroom ranking (scaled, from where you stand)

1. **Sig Jaccard (−0.016).** Structurally the most sensitive to count-matching, because Jaccard is punished symmetrically for over- and under-calling and your n_pred is off by 10× for most targets. Raw ceiling at matched counts ≈ p/(2−p) ≈ 0.4 vs a baseline that probably sits at 0.1–0.2 raw. Largest plausible scaled jump.
2. **Direction fidelity (−0.066).** Same lever, same reason. Caveat: verify in the metric code whether k requires sign agreement. If yes, your sign-correct precision at top-50 is ≈ 0.705 × ~0.63 ≈ 0.45 — *below* the shipped raw 0.49 — and the "cap at 0.71" is fiction. This one check can delete the whole plan; do it before anything else.
3. **Direction reach (0.186).** It's recall of the true DE set. You only find true genes efficiently at the top of your ranking; beyond N ~100 precision slides toward the base rate. Pushing reach fights 1 and 2 one-for-one. Leader's 0.892 comes from knowing the genes, not from calling more. Modest headroom, mostly tradeoff.
4. **LFC NMAE (0.116).** Measured ×2.03 inflation gave nothing. Only remaining lever is shrinking LFC on uncalled genes toward zero/typical reference magnitude; small and already partly realised.

PDS (0.706 → 0.892) and MSE require expression modelling, which you've correctly closed. Realistic ceiling if 1+2 work: composite ≈ 0.16 → 0.25–0.28. That is a large rank move and nowhere near 0.47.

## 3. Same-line replicate structure

The leader's MSE implies cosine ≈ 0.875 on a held-out line, which equals your measured cross-replicate ceiling (0.81–0.87). **Nobody models their way to replicate ceiling on a held-out line. The leader has, functionally, a measurement of the response.** Two legitimate sources: (a) public Perturb-seq for the same cell line with overlapping perturbations (Replogle-style genome-wide screens exist for several common lines), if external data is permitted; (b) any labelled perturbed cells of the held-out line present in the provided data (few-shot split), enabling same-line nearest-neighbour. Anything else is leakage.

I don't know the VCC 2026 rules or which line is held out; you do. **Minimum viable test (≤ 2 h, no GPU):** identify the held-out line; find any public dataset of that line; take ≥ 20 overlapping targets; compute the public mean response vectors; validate the *mechanism* on a line where you have truth (any training line with a public analogue): cosine ≥ 0.7 between public mean response and your training truth means lookup reproduces measured responses. Then read the rules. If allowed, this is the only path to the top of the board and everything in §1–2 is secondary. If not allowed, the top is unreachable for you legitimately and your goal for the next 7 days is rank, not podium — say so in the writeup.

## 4. What to build tonight, with threshold

**Build the faithful scorer + gated budget package.** Sequence:

1. (1 h) Local scorer that reproduces the shipped package's behaviour: call counts in the hundreds per target, raw fidelity ≈ 0.49, and raw→scaled mapping using the official anchors. Check the sign-agreement definition of fidelity in the metric code. **If you cannot reproduce call counts within 2× and fidelity within ±0.05, stop; the writeup says "local evaluation was invalid; here is the fix."**
2. (0.5 h) Stratified precision by n_