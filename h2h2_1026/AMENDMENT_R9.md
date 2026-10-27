# Amendment R9 — call-side arms from public research + sign audit (registered before scoring)

Authority: decisions/consult_1026_amp round 9 (2026-10-26 18:57 UTC). K250 and the reach head are unanimous;
detectability selection is a majority (Fable, Opus). K100 (Astra only) is not run. Built by build_r9.py, which is
label-free: it uses B0's response, B0's frozen confidence and the control reference mean CPM only. All arms use
B0's response verbatim, seed base 20260910, anchors and depth; sha256 in results/r9_arms.json.
- K250: B0 confidence, layout call budget 250 instead of 885.
- RH20: budget 885. B0's own top-885 eligible genes (by confidence) keep their membership. The 20 with the highest
  pct(detectability) x pct(confidence) take positions 1-20 in that order; the rest keep B0's order.
  Detectability = sqrt(control mean CPM) x |expm1(response)|.
- DET885: budget 885, eligibility unchanged; ranking among eligible genes = 0.75 pct(detectability)
  + 0.25 pct(B0 confidence) (Skadi mixs75). Median top-885 overlap with B0 is 0.13.
Runner score_r9.py = score_run.py confidence_all/run with textual substitutions only: confidence cache path; layout
budget per arm (BUDGET); render-json budget record.

Instrument (R9 Q2 composition; the local composite is reported but is NOT the reading):
- Strata from layout_cd4_1025/artifacts/per_target.csv n_conf (<10 / 10-999 / >=1000; calib3's definition).
  signaudit's p_adj count gives 9/44/26 and is reported for reconciliation.
- Per stratum: mean paired raw deltas vs B0 for fidelity, jaccard, reach and nmae.
- Official-scale conversion: fidelity x 3.30, jaccard x 2.88 (official slopes); reach and nmae at local scaled
  (conservative, Opus); PDS and MSE must be unchanged (see validity).
- Projection = sum over components of the stratum-weighted means with w = (0.315, 0.65, 0.035), divided by 6.
  Sensitivity sweep: strong in {0.02, 0.035, 0.05, 0.08, 0.15} x weak in {0.20, 0.315, 0.45}; the minimum is
  reported as "conservative" (Astra). calib3's global-ratio projection is a separate sensitivity check, not
  applied on top.
- Seed floor = the same transform applied to B0s1 and B0s2 vs B0 (max |projection|).
- Realised significant calls (including incidental ones) are reported per stratum (Astra).

Validity: rendered pseudobulk and PDS bit-identical to B0; raw MSE within 0.1%; raw nmae within ±0.001.
Screen PASS (median of the R9 bars): central projection >= +0.005 AND > 2 x seed floor (Fable). Reported but not
gating: share of weak+mid targets improved (Opus asks for >= 60%); conservative projection; global-ratio projection.
A PASS goes to CD4 replication (a separate amendment) and back to the advisors. Recommendation bar to Ricardo
(fixed now, by 10-28 per Opus): screen PASS + same sign on CD4 + DE-only validity + no strong-stratum loss larger
than its weight can absorb. Then "B0 + arm" goes to Ricardo with an honest probability, not a predicted size.
An official validation probe (Fable) is impossible: the board has been frozen since 10-22. No submission without Ricardo.
