# DECISION 2026-10-25 01:15 ET — MSE is NOT rho-limited: 289 competitors at our exact PDS reach 2.4x our direction signal

## The question that decides where the next 10 days go
`20261025_0045` derived `u_min = 1 - rho^2`, i.e. MSE is a pure function of the response-direction
correlation rho, and concluded MSE was "rho-limited". If that were the whole story, no modelling work
could help and only new data could. The leaderboard tests it directly, because for every entry the
public dump carries BOTH the raw `expr_mse_unbiased_capped_norm` and the raw `pds_cosine`.

## Result (1,635 val entries; 732 with a usable interior raw MSE)
Restating the raw MSE as an implied direction correlation, `rho_implied = sqrt(1 - u)`:

| population | implied rho |
|---|---:|
| top-10 median | 0.707 |
| **entries with raw `pds_cosine` within 0.03 of ours (n=289)** | **median 0.363, max 0.621** |
| us (measured locally at TS=50000) | **0.154** |

Correlation between `pds_cosine` and implied rho across the field is **r = 0.795**, so the two are
strongly related in general — *but not at our operating point*. Two hundred and eighty-nine entries
sit at our discrimination level and carry **2.4x more usable response-direction signal than we do**,
some 4x. The best of them reaches MSE 0.390 scaled, worth **+0.055 overall** against our 0.061.

## What this refutes
`20261025_0045` said "MSE headroom is rho headroom … no modelling knob can cross this … only new data
can." The first clause is right and the last is **wrong**. rho is not a fixed property of our data;
it is a property of *how we build the response*. Competitors with the same discrimination score are
extracting far more direction from the same public information. So the MSE gap is a **modelling**
gap, and it is the largest single opening we have (43% of the distance to top-10).

## What it does NOT license
- It does not say our data contains more signal — it says our *processing* discards or corrupts
  signal that others keep. The candidates are the ones already visible in our own code:
  the `de_call_mask` that zeroes 79% of the DE universe (`model_v5/v5.py:97`, `min_net>=1`), the
  50-control-PC `deflate`, and the 0.3 learned blend. Each suppresses response structure by
  construction. A competitor with identical PDS but far higher rho is most simply explained by
  **not** zeroing most of the response.
- It does not make the local harness a decision instrument: fidelity attenuates 10x, reach amplifies
  1.5x, and MSE must be scored on the same run as the incumbent.

## The concrete next experiment (top priority)
One emitted, officially scored arm that **stops suppressing the bulk response while keeping the DE
call set exactly as it is today**: skip `np.where(mask, L-Lc, 0)` for the *expression* profile
(`model_v6/pipeline.py:185` / `rsi/v5b_loco_baseline.py:118`) and skip the 50-PC `deflate`, leaving
`de_call_mask` only as the call-set generator. Everything else — head, amplitude, integer totals,
v9 layout, B/claim — frozen. Measure all six components plus B/claim, paired with the incumbent.

This is **not** one of the closed arms: the closed "density"/"bulk multiplier" sweeps scaled an
already-masked profile, and the closed "source breadth" sweeps added donors to an already-masked
profile. None of them removed the mask from the expression leg. `mse_direction_1024` built a dense
arm only as a *layout* comparison (`layout_native`), never as a profile change, and it therefore
scored identical MSE by `ratio_of_sums` invariance.

Predicted sign: rho rises, MSE rises, PDS is unchanged by construction (same group totals), fidelity
and reach are unchanged (same call set). The risk is NMAE — the same trade as `20261025_0105`.

## Status
Recorded; screening harness written in `emitter_screen_1025/` and not yet producing a trustworthy
number (two broadcast bugs fixed, then stopped when the arm set was recognised as re-testing closed
levers). **This needs a real emitted arm, not a screen.** Safety entry `5NeeoWn96lfm23hgExcm`
unchanged; no package change, no submission.
