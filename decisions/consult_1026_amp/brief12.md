Round 12 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = billmA 0.16438 on validation; final phase scores new contexts D/E/F blind; nothing is
submitted without Ricardo). <= 250 words. R11 ruled B0 final (2-1). Ricardo asked "how can we do better biology?";
four CPU diagnostics followed. Read bio_gap_1026/REPORT.md, bio_gap_1026/shared_resp/REPORT.md,
bio_gap_1026/align/REPORT.md, bio_gap_1026/loco_map/REPORT.md, h2h2_1026/scratch/mse_instr/REPORT.md.

## Findings (168 H1 training targets, LOO; 79 held-out never used)
1. No local MSE instrument: all six formulations fail C4 (clamp lifted -0.0012 vs official +0.0791).
2. B0's pseudobulk delta: E=1.64, mean cosine 0.084; 58% direction / 42% magnitude error. Covered targets carry 99.5%
   of the error, so new source coverage buys <=0.5%.
3. Every SSE gain is shrinkage. With cross = sum(bhat.b)/sum||b||^2: B0 0.1154; family correction 0.1086 (SSE -10.5%,
   spread 80%); shared mean + B0 0.064 (PDS 0.52). Under budget-limited credit (B/claim ~0.5-0.6) net error ~
   T - 2*cross, so only cross at matched spread counts.
4. In-context arms fitted on TRUE H1 deltas of other targets: kNN cross 0.423 (PDS 0.578), ridge map 0.398 (PDS
   0.662). Not usable: D/E/F have only control cells.
5. Leave-one-recipient-context-out (H1, K562, RPE1, HepG2, Jurkat): full maps fail PDS everywhere and lose cross on
   3/5. PDS-bounded blends pass 5/5. Best on H1 vs B0: cross 0.1154 -> 0.1181 (CI +0.0020, +0.0034), PDS -0.0023
   (CI -0.0052, +0.0004).

## Orchestrator arithmetic (check it)
Treat dcross as a raw-MSE gain of 2*0.0027 and convert at 0.811 scaled per raw unit: +0.0044 scaled MSE, about
+0.0007 overall. PDS -0.0023 local raw -> -0.0048 local scaled, official/local PDS ratio 0.27-0.6, about -0.0003
overall. Net about +0.0004, below the +0.0028 seed floor. Reading: big biology gains need in-context truth that the
final phase does not have.

## Rule
1. Is the arithmetic right? Does any finding justify building an arm (control_bounded or other) before 11-03?
2. Is there a biology lever these diagnostics missed that transfers to unseen contexts and is CPU-feasible in 7 days?
3. If not: close the biology line, B0 final? Say plainly.
