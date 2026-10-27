Round 11 of consult_1026_amp. You are a science advisor to a Virtual Cell Challenge 2026 team (deadline 11-05 23:59
UTC; frozen entry B0 = billmA 0.16438; final phase blind; nothing is submitted without Ricardo). <= 250 words.
Read dm_1026/REPORT.md.

## Measurement contradicts the R10 composition
I composed R10 with Astra's mechanism (per-cell integer depth multipliers: CPM vectors fixed, so DE and NMAE are
invariant). The 48 h feasibility test killed it in two hours. On 5 hash-verified B0 targets, the best threshold
weighting amplifies the pooled delta 1.03x with 0.20 orthogonal error. Even unconstrained NNLS weights aimed at
1.38x reach 1.00x with 0.45-0.49 orthogonal error. The emitted cells are sampling noise around B0's profile, and
no whole-cell reweighting moves the pool along the delta.

## What remains
Opus/Fable's depth-tilted placement: per gene, shift counts toward high-depth cells with per-cell zero-sum
compensation (borrowing from uncalled genes). The pooled sum rises and the per-cell CPM mean holds. Efficiency
1 - d_lo/d_hi: emitted depths p10/p50/p90 = 34.9k/55.6k/85.1k, so 0.36 (IQR) to 0.59 (10/90). Arithmetically
feasible. It changes within-group CPM dispersion, so Wilcoxon ranks (DE calls, reach, fidelity) may move, and the
10-22 dual-moment kernel failed its tolerance on the pooled template. About 1.5-2 days to build and score one arm
(local PDS, NMAE and DE readable; MSE projected from the official v8a->v8b slope plus label-free B/claim).

## Rule
1. Build the depth-tilt arm (one arm, gate as in R10 but with DE "within the seed floor" instead of
   bit-identical), or stop the decoupled-amplitude line now?
2. If build: which pooled gain s (1.25 or ~1.38), and should called genes be excluded from the tilt so the
   layout's DE goals stay untouched?
3. If stop: B0 final? Say plainly.
