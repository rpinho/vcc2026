VERDICT: Yes—read the metric implementation before modelling. Fidelity, not cosine or expression MSE, is the immediate bottleneck.
EXPERIMENT: Cross-family, **sign-only correction** of the frozen MLP’s predicted gene responses; preserve effect magnitudes.
PROTOCOL: Audit fidelity’s numerator, denominator, significance filtering, zero handling, and anchors. Train a CPU classifier on response-bank features to predict whether each MLP response needs its sign flipped. Use nested leave-one-family-out validation, selecting the flip-confidence threshold with the official scorer. Keep magnitudes and residual sampling fixed.
KILL TEST: Reject if held-out official fidelity and composite do not improve over the unchanged MLP, or if apparent gains reflect denominator collapse rather than better directional calls.

Justification: The supplied results identify the optimization target, but not its exact mechanism. A metric named “direction fidelity yield” need not be ordinary sign accuracy; reading its definition is essential before choosing training labels.

Use bank sign agreement, response reproducibility, gene identity, and control-expression context as classifier features, excluding held-out-family responses. Generate training predictions out of fold to avoid learning unrealistically easy corrections. Derive labels from the scorer’s actual direction rules.

This tests a specific hypothesis: transferable response information can correct wrong directions without shrinking effects. Report fidelity, reach, composite, and the audited numerator/denominator alongside a freshly scored trivial baseline; use paired resampling at the perturbation level to assess stability.