# Decision request — modeling project, 11 days left

I run a scientific modeling project with a hard deadline. One deliverable is already finalized and
published; only the most recent version counts. The scoring server for new candidate versions is now
closed, so no new version can be measured before the deadline.

Current finalized version: validated score 0.164 on a 0-1 scale (leader 0.47, top-50 median 0.315).
Recent attempts to improve it all measured worse.

A local test harness is used to rank candidate approaches. An audit just found:

1. Harness rank accuracy is fold-dependent: 15/21 correct on one fold (10% overlap with the real
   target set), 18/21 on another fold, 9/10 on a pooled 5-item macro.
2. The harness reports my current model's direction correlation as 0.15. The official score, converted
   through measured calibration anchors, implies about 0.31. So the harness may understate the
   quantity that matters most by ~2x. Cause unmeasured; candidates include fold difficulty, an
   emitter/baseline mismatch, or a bad anchor conversion.
3. Of 13 ranked research ideas from an earlier survey, 4 were never implemented — and they are exactly
   the ones that could raise direction correlation (feature-based regression; source selection by
   held-out response reliability; priors for targets with no source data; routing between models).
   Everything actually tried was an output-scaling/layout knob, which theory says cannot raise the
   correlation ceiling.
4. One experiment aborted at 63% due to a disk-space guard. Its result is void, not negative.

Constraints: one CPU machine, one part-time operator, no human monitoring overnight, and the finalized
version must not be put at risk.

Question, answer decisively in <=600 words:
(A) VERDICT: keep the finalized version unchanged and stop | keep it but do bounded diagnostic work |
    pursue a new version
(B) WHAT TO DO, ranked, with hours
(C) Did facts 1-4 change your earlier advice? If not, say so plainly.
(D) RISK
