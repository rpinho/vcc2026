# Amendment RENDER — render-seed noise floor + pure-scaling amplitude map (registered before scoring)

Authority: decisions/consult_1026_amp round 4 (2026-10-26 16:45 UTC). All three: B0 frozen as incumbent/default;
research continues bounded; pure-scaling ladder = exploration MAP (Fable, Astra: no promotion; majority);
renderer becomes the research target, first experiment = render-seed replicate (Fable + Opus).
Alphas = per-slot median of the three proposals {0.75|0.7|0.8}, {1.25|1.3|1.2}, {1.5|1.6}: **0.75, 1.25, 1.5**.
Research cutoff 11-02 23:59 UTC (Opus); after it B0 is the entry unless a gated arm passed on seed replicates.

Arms (results/predictions_<arm>.npz; sha256 in results/render_arms.json):
- B0s1, B0s2: B0's response verbatim, renderer anchor-emission seed base 20261101 / 20261102 instead of 20260910
  (per-target seed = model.atlas.seed_for('H1:'+t, base), as for B0). Everything else identical.
- A075, A125, A150: alpha * B0 (float64 multiply, stored in B0's dtype); signs identical to B0 (verified 100%),
  so sign-dependent confidence and the 885-call layout eligibility are identical; seed base 20260910.
Runner score_render.py = score_run.py confidence_all/run with textual substitutions only (cache path, arm list,
seed base per arm). Order: B0s1, B0s2 first (noise floor), then A075, A125, A150.

Readout (descriptive; NO promotion from this amendment):
- Noise floor: |official(B0sk) - official(B0)| for score_avg, PDS, fidelity, MSE; paired per-target spread.
  If seed-to-seed |dPDS| >= 0.01, verdicts with |dPDS| < 0.01 (N1's -0.0104 included) are flagged as within render noise.
- Amplitude map: official six components + paired deltas vs B0 for each alpha; raw MSE reported (local scaled MSE
  is clamped; leaderboard MSE credit 0.061 at risk for alpha > 1).
Any candidate arising from this map needs its own preregistered gate (Opus's alpha>1 gate recorded for that case:
composite LB > 0, raw MSE <= B0's, fidelity drop <= 0.01, holding across seed replicates).
