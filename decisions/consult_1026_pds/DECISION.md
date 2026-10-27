# DECISION — 2026-10-26 08:15 EDT. Round 2: the PDS finding. All three advisors consulted.

Ricardo's rule: the advisors decide. This records what each said and what therefore runs.

## Unanimous, and it reverses Round 1
**"Contest expression MSE" is dead as a modelling objective.** All three advisors explicitly abandoned
their own Round-1 instruction once they saw the official numbers:

- **Fable:** "Kill 'contest expression MSE.' It was wrong before today and it's dead now. ... Fable and
  Kimi were reasoning about the raw metric, not the scaled composite. The arithmetic you wrote out is
  correct and it overrides them."
- **Astra:** "The advisors' objective was wrong for these measurements. ... Stop treating expression MSE
  as the objective."
- **Kimi:** "The advisors were half wrong. They were wrong to make expression MSE the objective. They
  were not thereby vindicated on DE."

The triggering evidence: B1 cut raw MSE by 0.46 (beating even the anchor baseline) and collected only
**+0.00164 points**, while its PDS fall of 0.295 cost **−0.10282 points** — a 63x bigger loss.

## Where they land, independently

| | Fable | Astra | Kimi |
|---|---|---|---|
| objective | (c) DE fidelity, PDS as hard constraint ≥0.84 | (a) rank, full composite, protect PDS | (c) defend PDS + enforce a point budget |
| Monday experiment | amplitude/shrinkage **+ two DE arms** (additive shared program; variance shaping) | **amplitude calibration**, then package gate | **amplitude calibration** |
| PDS ≥ 0.80 enough? | No — 0.84 floor | No — break-even is ~0.84533 | No — break-even ~0.845 |
| reliability weighting / Sidechain | paused | paused until scalar resolves | paused |
| GPU | no | no | no |
| 0.25 | unreachable (ceiling ~0.13) | unreachable (~0.1315) | unreachable (~0.133) |

**Verified their break-even arithmetic independently:** PDS points per raw unit = 0.3481; 0.85→0.80 costs
0.01741 pts; a B1-sized MSE gain needs PDS ≥ 0.84529 merely to break even. All three agree, and my own
computation reproduces it to 4 decimals.

## Therefore — Monday's first experiment (unanimous)
**Amplitude-only calibration of the shipped model**, ahead of reliability weighting and Sidechain:

    y_hat_alpha = mu_train + alpha * (y_B0 - mu_train)

with alpha per family, fitted and frozen by **leave-one-out on the 168 only**, never swept on the 79.
Positive scalar scaling cancels in a residual cosine, so PDS *should* be scale-invariant — but Astra
correctly notes that is NOT automatic here, because the renderer exponentiates, clips, renormalises and
emits integers, so scaling can rotate the final scored vector near zero and can move DE coverage and
significance. **PDS preservation must be measured, not assumed.**

Fable adds a companion: if a scalar > 1 raises fidelity in the LOO loop, that is free evidence of
under-calling.

## Gates (converged, stricter than Round 1)
- Package only if the **complete official component ledger nets positive** and `score_avg` strictly
  beats **0.021555**. Ties ship B0.
- PDS floor **0.84** (Fable) / break-even **0.8453** (Astra, Kimi) — treat 0.845 as the real floor.
- Positive paired-bootstrap CI lower bound on the composite, honor the existing preregistration, do not
  weaken it after seeing B2.
- Untouched confirmation fold before promoting a small gain — the 79 have already informed strategy.

## Where they disagree
Fable wants **two additional DE arms** Monday (additive shared program; variance shaping), arguing
B1's +0.014 fidelity gain while PDS was at chance proves a shared pan-perturbation program carries
fidelity/jaccard independently of the specific residual — "you never tested keeping both." Astra and
Kimi both say DE headroom is real but we have **no reliable mechanism** to collect it, and want the
scalar experiment alone.

**Resolution:** the scalar experiment runs first (unanimous). Fable's Arm 1 (additive shared program)
is cheap and directly testable, so it runs **second**, behind the scalar and behind the same PDS-budget
gate. Variance shaping waits — it is the most speculative of the three and touches the renderer.

## Immediate consequence for the running job
B2–B5 finish as **diagnostics**. No successor manifold or DE arm is authorized by this round. If B2
fails the composite gate, the manifold/STRING direction is closed permanently (all three agree).
B5 (shrink 0.5) is now doubly important: it is the immediate `alpha`-isn't-free diagnostic — if pure
shrinkage itself spends PDS/DE points, the whole calibration premise is falsified and we freeze B0.
