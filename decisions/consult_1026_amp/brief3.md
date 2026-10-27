Round 3 of consult_1026_amp (short). You are a science advisor to a Virtual Cell Challenge 2026 team; the
orchestrator executes your rulings. Answer in <= 200 words, one decision per question.

## Round 2 outcome (unanimous parts, now executed)
- N1 only (A = mu168 + 1.25(B0 - mu168), rescaled per target to ||B0_i||, self-gene from B0). Built and hash-locked
  (h2h2_1026/AMENDMENT_N1.md). Opus kill rule: response-MSE ratio N1/B0 on the 79 = **1.00805** (< 1.055): ALIVE.
  Official N1 run is queued behind the disk guard (machine at 30.5 GB; launch needs 31.5 sustained).
- 79-gate for N1: score_avg > 0.021555, PDS >= 0.845, paired composite LB > 0, official raw-MSE ratio N1/B0 <= 1.03.
- **B3 paired bootstrap (registered machinery, 20k draws): composite +0.00406, 95% CI [-0.00197, +0.01048];
  PDS +0.00860, CI [+0.00097, +0.01688]; fidelity CI includes 0; raw MSE ratio 0.9991. Lower bound <= 0 =>
  B3 DROPPED** under all three of your round-2 rules.
- P1a build is running (p1a_1026/), with your constraints recorded.

## CORRECTION: the confirmation route two of you chose is infeasible
Opus and Astra chose (d), a Ricardo leaderboard submission compared against B0's leaderboard score. The
competition is in its FINAL phase: final-phase submissions are scored on the hidden panel vcc2026-test-1 and
their scores are hidden until the end; the public validation board has been frozen since the validation phase
closed (1,635 entries, unchanged since 10-22; competitors/DAILY-2026-10-26.md). A submission therefore returns
no score and cannot confirm or refute anything before the deadline. (Ricardo can still choose to submit; it is
a blind bet, not a test.)

Remaining options: (b) Fable's CD4 fold check — FEASIBLE but needs a build: mse_1025/inputs has the shipped CD4
predictions (inference_CD4.npz), the CD4 search h5ad and the official CD4 baseline; a CD4 training mean must be
taken from CD4 truth not in the CD4 search fold (availability to be checked); est. 2-4 h CPU to adapt the runner.
(a) no confirmation, promote on the 79 gate plus a margin (Fable: >= +0.008 over B0; Opus: 79 gate + raw-MSE
ratio <= 1.03).

## Rule
1. Confirmation: (b) CD4, or (a) with which margin? If (b) and the CD4 training mean cannot be formed legally,
   what then?
2. If N1 passes, what exactly goes in the recommendation to Ricardo (he decides whether to submit a blind final
   entry; we have at most ~2/day until 11-05 23:59 UTC)?
