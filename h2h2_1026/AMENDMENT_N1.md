# Amendment N1 — frozen before any N1 official scoring (2026-10-26 ~14:35 UTC)

Authority: decisions/consult_1026_amp/ round 2, unanimous (Fable, Opus, Astra): run N1 only; E0/E3 dead.
This amendment does NOT modify PREREGISTRATION.md, its hash lock, or the B0-B5 verdict (still REJECT-bound
for manifold/STRING, pending B5 + verify).

Arm N1: A = mu168 + 1.25 (B0 - mu168); N1_i = A_i * ||B0_i||/||A_i||; self-gene value copied from B0.
mu168 = results/fitted_model.npz 'mu' (168 training targets only) in scoring order; sha256 of float64 bytes
e63ab57b0720f4deaa5b978467acf43d88e8855b88ee618f8eb3be4b69a59010. predictions_N1.npz sha256 recorded in
results/n1_build.json. Alpha 1.25 frozen on the 168 (training-only LOO) in round 1; not re-optimised.

Rendering/scoring: identical to B0 (score_run.py functions, same renderer, budget 885, seed 20260910, anchors,
reference cells, warm DE cache). Only change: sign-dependent donor confidence recomputed for N1 with the
identical formula (cached to scratch/confidences_N1.npz; the B0 confidence hash is re-asserted in the same pass).
Runner: score_n1.py (generated from score_run.py by two textual substitutions, saved as score_n1_generated.py).

Gate (all must hold, binding, from the round-2 rulings):
1. official score_avg(N1) > 0.021555 strictly
2. raw PDS(N1) >= 0.845
3. paired-bootstrap (registered analyze.py machinery, 20,000 draws, seed 20261026) composite N1-B0 lower bound > 0
4. official raw expression-MSE ratio N1/B0 <= 1.03 (Fable/Opus; local scaled MSE is clamped and cannot see it)
5. response-MSE ratio N1/B0 on the 79 < 1.055 (Opus kill rule) — measured 1.00805 at build: PASS
Passing the gate makes N1 a development-panel pass only. Confirmation protocol is under advisor round 3 (the
leaderboard route chosen by Opus/Astra is infeasible: final-phase scores are hidden, validation board closed).
