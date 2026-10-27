## Decision

**Ship B0. Stop contesting expression MSE. Monday's single line of work is making the DE caller stop choking on our rendering while leaving B0's predicted mean untouched. GPU: no.**

Reasoning, by your questions:

### 1. "Contest expression MSE" is the wrong objective. The arithmetic kills it.

B1 cut MSE by 33% (1.41 → 0.95), beat the anchor, and earned **+0.0016**. Extrapolate: to earn the +0.08 sitting in DE fidelity you'd need MSE ≈ 0.5, i.e. halve the anchor. Nobody is doing that with a model that currently sits at 1.41. The MSE term is flat near baseline by construction; it is a dead lever at any achievable effort. Fable and Kimi read "expression MSE" as the headline metric. It isn't. The headroom table is the objective function; they didn't read it.

Their other half — "stop DE work" — is *half* right, and this is the reason your three DE experiments failed. All three (effect-size gating, recipient conditioning, manifold shrinkage) changed the **predicted mean**. But look at the contradiction in your own table: direction is good enough for PDS 0.85 (far above baseline), yet DE fidelity is −0.48 scaled. You cannot have the right direction and catastrophically wrong DE gene sets *because of biology*. That combination has exactly one common cause: the DE caller is being run on rendered cells, and our rendered cells have the wrong within-perturbation distribution (variance too low / dropout structure wrong / cell count off), so every small deviation tests as significant, DE sets balloon, fidelity and Jaccard crater. Gating made fidelity *worse* (0.45 → 0.15) because it zeroed real effects while leaving the variance problem intact.

**Replacement objective: recover DE fidelity + Jaccard (≈ +0.11 available) through the renderer, not the mean.** Everything that touches the mean is now frozen.

### 2. B2 gate

PDS is worth roughly 0.35 points per raw unit (B1: −0.295 raw → −0.103). MSE gain is capped at ~+0.002. So B2 at PDS 0.80 nets **≈ −0.015** and loses. The gate is not "PDS ≥ 0.80." The gate is: **official composite > 0.021555 on the same 79, same anchors.** Nothing else. Don't ship on a PDS proxy, don't ship on a response-space ratio — response space already told you B4 > B2 on MSE and that B0's scalar is 0.10, and response space is blind to the only term that matters. My expectation: B2–B5 all lose or tie, B0 stays.

If B2 collapses, Monday is the plan in §1, specified in §3.

### 3. Amplitude-only rescaling: run it as a 30-minute diagnostic, not as Monday's headline, and not at s = 0.10.

Two cases, both bounded:
- If PDS is computed on deltas from control, cosine is exactly invariant to a positive per-perturbation scalar. PDS unchanged, MSE gain ≤ +0.002. Harmless, worthless.
- If PDS is computed on absolute pseudobulk, s = 0.10 drags every prediction to within 10% of control; the 79 become indistinguishable; PDS collapses like B1. B1 and B5 are already this experiment.

Upside ≤ +0.002 in the good case, −0.10 in the bad case. It is not ahead of anything. Its only value is telling you whether PDS is scale-sensitive, which you need to know before the rendering work. Do one arm, s = 0.5 per family, LOO on the 168, read PDS, move on.

**Monday, in order (all CPU):**
1. Read B2–B5 composites. Ship the max. Done by 10am.
2. Diagnostic: take B0's exact predicted means for the 168, run the *official* DE caller on our rendered cells vs. on true cells. Report per-perturbation: #DE predicted vs #DE true, p-value distribution, LFC scatter. If predicted #DE ≫ true (I'd bet on 5–50×), that is the entire fidelity/Jaccard story and nothing upstream needs to change.
3. Arm C — dispersion-matched rendering: render each perturbation by taking real control cells and applying B0's delta (additive in log space or multiplicative in counts — try both), preserving control variance and dropout; match true cell count. Pseudobulk mean is identical to B0 by construction, so PDS and MSE should move < 0.01. Score on the 79 in-house. Gate: fidelity up, PDS within 0.01 of 0.85.
4. Reliability weighting and Sidechain go *behind* C unless someone can show me they leave the pseudobulk mean invariant. If they touch the mean they're in the same bucket as the three failures.

### 4. Objective: (a) rank.

0.25 is not on the table: full recovery of every negative component to baseline gets you to ~0.16, and that requires beating the anchor on DE sub-metrics you're currently 0.48 scaled below. Component-chasing is what produced B1. Rank it is, and rank is won by not giving back PDS.

**STOP, effective now:** expression-MSE arms; manifold shrinkage of any F; effect-size gating; recipient conditioning; using response-space ratios as a decision criterion; any experiment whose output changes B0's predicted mean.

### 5. GPU: still no.

The mean is the one thing that works and it's frozen. Every remaining lever is post-processing of a fixed mean through a renderer and a DE caller. That's a CPU afternoon. If someone asks for a GPU on Monday, they're proposing to retrain the direction model, and the answer to that is also no.