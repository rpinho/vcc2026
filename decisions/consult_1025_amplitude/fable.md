VERDICT: Yes, run amplitude calibration next, but only as a cheap test scored on the official metric. Most of the measured gain is shrinkage toward zero, not real signal, so expect a modest leaderboard effect.

EXPERIMENT: Replace "rescale to the pooled direction's norm" in the incumbent with "rescale to α × pooled norm." Fit α as one global scalar per training fold, nested leave-one-family-out.

PROTOCOL: Preregister the α grid {0.25, 0.3, 0.35, fitted}. Score each held-out family with the full official scorer and report every component: debiased capped MSE, PDS cosine, signature Jaccard, direction fidelity and reach. Compare against the incumbent and an α=0 control using paired family-level bootstrap CIs.

KILL TEST: Kill the change if the official macro score's paired 95% CI includes zero. Also kill it if any DE component drops more than its noise band, or if it does not beat α=0 by more than 2 percentage points on the expression term.

Zero response alone gives −8.2%, so the scalar adds only about 4 points beyond shrinkage. Cosine is scale-invariant, so PDS will not move. Thresholded Jaccard and reach may degrade. The 0.16 → 0.47 gap is a direction problem: per-gene cosines are 0.01–0.22. After this quick check, the next experiment should target direction quality, not amplitude.