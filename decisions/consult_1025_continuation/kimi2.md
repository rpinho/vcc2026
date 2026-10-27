VERDICT: Pursue mechanism (1) per Advisor A: source-reliability weighting of donor pooling.
EXPERIMENT: Implement per-donor reliability weighting and validate via held-out-family rendered expression MSE.
PROTOCOL: Estimate donor residual variance on held-out families; apply inverse-variance weights to pooled direction; re-optimize blend factor; compute MSE.
KILL TEST: Terminate if held-out MSE improvement <5% relative or reliability variance collapses (CV<0.05) within 0.5 CPU-days.

Expression error dominates, so fixing donor pooling offers the highest leverage. Mechanism (1) is low-risk, requires minimal architectural change, and fits the 2 CPU-day budget. Held-out-family MSE provides a defensible local proxy without the official scorer. Mechanisms (2)-(4) demand more data or compute for uncertain gain. The kill test ensures early exit if donor reliability signals are uninformative.