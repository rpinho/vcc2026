**My call: `h2h2_1026` is the gate. Start nothing else while it runs. A noisy tie counts as failure.** If the manifold fails, STRING is dead and Monday belongs to Sidechain—not reliability weighting. No B0 bug hunt, and no GPU.

### 1. If the manifold loses: run Sidechain staged source progression

**Protocol.** Freeze the Sidechain source order, features, and hyperparameters before evaluation. Fit using only the 168 H1 training targets. Decide source retention and the stopping stage by target-held-out CV within those 168; do not use the 79 to choose stages. Then freeze one model, refit on all 168, predict the 79 once, and score through the official pipeline.

**Kill test.** Sidechain survives only if that single frozen model:
- lowers expression MSE relative to B0, and
- improves official total score over B0 with a 95% paired target-bootstrap CI entirely above zero.

Otherwise kill it. No retuning after seeing the 79, and no celebrating a CV-only or cosine-only improvement.

**Why reliability weighting waits.** I have evidence of donor-response **mismatch**, not evidence that donor-response **noise** is the limiting factor. Reliability weighting can confidently upweight a reproducibly wrong donor and adds another fitted set of knobs. Sidechain directly tests whether each source adds out-of-sample value and allows harmful sources to be rejected. That is the missing measurement.

### 2. If the manifold wins: do not package merely for beating B0

F3/B2 is the candidate; do not shop among arms afterward. Before building a package, it must:

- beat B0 on **official total score**, with the paired 95% target-bootstrap CI excluding zero;
- improve expression MSE without giving back the DE fidelity/reach performance;
- beat the training-mean control and add value beyond B5’s shrinkage;
- reproduce from training-only inputs in a clean inference path, including a deterministic fallback for the two targets lacking STRING coverage.

**The first invalidator is B5.** If shipped-B0-shrunk-toward-the-mean matches F3, the apparent gain is ordinary shrinkage, not manifold signal. Do not build a manifold package. Likewise, a cosine victory accompanied by a flat or worse official score is a loss.

### 3. B0’s 0.4487 direction fidelity is a transfer warning, not a bug investigation

Do **not** stop everything to investigate it.

Out-of-sample prediction can be genuinely below chance when donor response signs anti-transfer to H1, and raw sign agreement can be dominated by genes with near-zero effects. Moreover, 0.4487 versus 0.4929 is two summary numbers, not a uncertainty analysis; genes within targets are correlated. It is not evidence of a coding error without failed gene/target/sign-alignment assertions or comparable pathology on data B0 actually fit.

It also does **not** establish that the whole submission lacks DE skill: raw all-gene direction on this H1 panel is not the official DE score. B4 and B5 already provide the useful diagnosis. If either captures most of B0’s official score, then B0’s panel performance is indeed mostly basal-expression/mean-shrinkage value rather than target-specific response skill. If they do not, raw direction is simply the wrong lens.

### 4. Objective: **(a) move up from rank 463**

Expression MSE is the sole battlefield; it is not a separate vanity objective. Nothing measured yet justifies promising approximately 0.25.

Stop all DE direction/fidelity/reach projects, stop reopening the three closed lines, and stop treating component or cosine improvements as progress unless official total score moves.

### 5. GPU: **no**

Everything under consideration is cached-matrix fitting, source ablation, bootstrapping, and small-scale inference. A GPU removes no current bottleneck; it would merely enable a new model family that cannot be validated responsibly within 48 hours. Spend nothing.

The frozen safety entry stays last, and nothing is submitted without Ricardo’s authorization.