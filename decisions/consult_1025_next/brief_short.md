# Decision request — what work, if any, should run for 11 days

Prior context: three independent advisors twice agreed 3-0 to keep one finalized, published result and
build no replacement. The scoring server for new entries is closed, so nothing new can be measured
before the deadline. That decision stands; this question is only about how to spend idle machine time.

State, verified: one finalized result published (score 0.164 on 0-1; leader 0.47). Package hash
unchanged. One CPU machine, one part-time operator, no overnight supervision. 11 days remain.

New findings from a completed 20-minute audit job:
1. The most informative next experiment — replaying historical approaches on a test fold with 3x better
   target overlap — is NOT runnable: required baseline matrices, calibration anchors and frozen
   configurations are missing from retained artifacts. So "is the local test instrument fixable" stays
   unmeasured, not answered.
2. One of four proposed untested mechanisms was already implemented earlier and rejected. Only a
   distinct feature-combination variant remains untested.
3. Feature coverage for the 161 targets with no direct source data is complete; what is missing is
   representative ground truth for them.

The audit's own suggested order: mechanism #7 (source weighting by held-out response reliability)
first, then #6 (feature-combination model), then #11 (priors for unsupported targets), then #12.

Complication you should weigh: the orchestrator produced a status report today containing fabricated
facts — a command that does not exist, a monitoring schedule that does not exist, an audit that never
ran, and a result ID that does not exist. The operator caught it. The finalized result was never at
risk and nothing was submitted, but intermediate status reporting is not trustworthy.

Question, answer decisively in <=600 words:
(A) VERDICT: what should run next
(B) WHAT TO DO, ranked, with hours
(C) WHAT SHOULD NOT RUN
(D) HOW SHOULD PROGRESS BE VERIFIED without trusting the orchestrator's summaries
