# Amendment P1A — official development scoring of the frozen P1a calibration (registered before scoring)

Authority: consult_1026_amp R1-R2 (unanimous: P1a is the primary lever, frozen before development scoring,
measured officially, composed N1 -> P1a, fitted on the 168 only). Does not modify PREREGISTRATION.md or AMENDMENT_N1*.
Config: p1a_1026/FROZEN_CONFIG.json sha256 9bc09faaa70b91206ccb53fb24ab9bd625510280a7d49568bf7817ca940fc949.
Arms (copied verbatim from p1a_1026):
- P1aB0 = P1a(B0)   results/predictions_P1aB0.npz sha256 8e825d1c720128366b6230878630ab2a8e0f28656709f8195c8bf16f76d3fd30
- P1aN1 = P1a(N1)   results/predictions_P1aN1.npz sha256 b1715467e6670a3b9c5078b54417d6e8e6cd016396380dccd03359b986121181
Rendering/scoring identical to B0; sign-dependent confidence recomputed with the identical formula.

Gate for each arm vs ITS base (P1aB0 vs B0; P1aN1 vs N1) — all binding:
1. official score_avg > 0.021555 and > its base's score_avg
2. raw PDS >= 0.845
3. paired composite (20k draws, seed 20261026) lower bound > 0 vs its base          (Fable)
4. official raw expression-MSE ratio vs its base <= 0.97                             (Fable; 0.97-1.0 = noise)
Training-side rules (Obj > 0, r < 1, rank guard <= +0.5; Opus) were met at freeze (p1a_1026/REPORT.md).
P1aN1 is the composition candidate only if N1 itself passes its gate (AMENDMENT_N1.md); otherwise P1aB0 is.
A passing arm is a development-panel pass only; confirmation follows AMENDMENT_N1_CD4.md's route (CD4).
