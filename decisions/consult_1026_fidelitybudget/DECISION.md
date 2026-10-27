# DECISION — call-budget (fidelity) line — 2026-10-26T00:25Z

Advice on record in this directory: `kimi.md`, `astra.md`, `fable.md` (fable truncated at finish_reason
`length`; its objection is resolved below). Measurement overrides advice where they conflict.

## What we measured (authoritative)

`gate_1025` Arm A finished and was ACCEPTED. Rendered at **B=885/target** on CD4:

| metric | incumbent | Arm A |
|---|---:|---:|
| fidelity (raw) | 0.1767 | **0.5706** |
| reach (raw) | 0.2796 | **0.2397** (WORSE, -14.3%) |
| jaccard (raw) | 0.01623 | 0.01641 |
| lfc_nmae (raw) | 0.93928 | 0.93925 |
| pds | 0.4974 | 0.4974 |
| expr_mse_capped_norm | 1.2862 | 1.2383 |

**Realisation = 89.17%** (62,342/69,915), newly admitted 89.79%. The correction-2 premise
(`>30%` realisation) is comfortably supported. All 32 previously-silent targets were rescued:
fidelity 0.0000 → 0.5648.

## What this does and does not license

The **mechanism** is proven: widening the call set is realisable in 9 of 10 calls and turns
structurally-zero targets into scoring targets. That is a genuine, first-of-campaign result.

The **benefit** is NOT proven. Three independent reasons:
1. Reach got **worse**, and reach is already one of our strongest members (official 0.186 scaled).
2. Local anchors are degenerate (fidelity span r-b = **0.0130**), so local `score_avg` is meaningless.
3. No held-out family was tested at the chosen B.

## Advisor positions and how they resolved

- **Kimi — B=75 + global LFC x2.03 "fixes NMAE to ~0".** Measured: NMAE at x2.03 = **0.9171**, i.e.
  the null predictor (1.0000). Global LFC inflation is **rejected**. It also sits squarely in the
  already-dead amplitude family.
- **Fable — "B=500 should crush fidelity; you are inconsistent."** Resolved: once `n_pred > n_conf`
  the denominator becomes `n_pred`, so fidelity = precision and asymptotes at the sign-agreement rate
  (~0.57). It does not collapse. No inconsistency; Fable's objection is answered, not accepted.
- **Astra — selection leak in my sweep; use B=100; NMAE gains nothing directly; reach is not recall.**
  All three correct. Sweep redone with **prediction-only selection**: B=100 → 0.5692, B=150 → 0.5878
  (curve holds). Astra's caution that "power, not sign" overstated the evidence is accepted — and
  confirmed by the fact that a *random* gene set already yields fidelity ~0.47 (our signs are
  individually right 57% of the time regardless of which genes we name).

## The unresolved scientific caveat (most important)

Our `p_adj` ranking of *which* genes are DE is near-random: a true-DE gene's median rank is
0.389 of the list (random = 0.500), and Jaccard stays ~0.007 at B=75. So the budget buys fidelity
largely through `k/max(n_pred,n_conf)` arithmetic and our above-chance sign accuracy, **not** by
finding the right genes. Jaccard (1/6 of composite) stays at its floor. This is why the budget
cannot be expected to deliver the full fidelity translation.

## DECISION — what runs next

Per Astra §5, run **one rendered budget sweep**, incumbent vs B=100 and B=150, predictor/cells/
totals fixed, selection from predictions only, scored by the exact six-metric scorer, then confirm
on one untouched held-out family.

**Pass:** >=80% realisation, >=+0.02 composite vs incumbent on the confirmation family, paired
target bootstrap 95% lower bound above zero.
**Fail:** reject the budget change.

Frozen entry `5NeeoWn96lfm23hgExcm` untouched. No submission. Nothing was pushed.

---

## ADDENDUM 00:32Z — the local "incumbent" is probably NOT a fair baseline

The Arm A delta is measured against a **local** arm that calls a median of **1** gene per target
(33/79 targets call nothing). If the **shipped** package behaved that way, official fidelity would be
`~1/30 = 0.033`, scaled `~-1.6`. Observed official billmA scaled fidelity is **-0.066** (raw ~0.492).

**Those are irreconcilable.** Either the shipped package calls far more genes than the local
incumbent arm, or the hidden panel's `n_conf` is much smaller. Either way:

> **The +0.39 raw fidelity gain (0.177 → 0.571) is measured against a pathological local baseline,
> and substantially overstates what packaging a wider call budget would buy on the leaderboard.**

This is the third time tonight that a large-looking local delta has shrunk on contact with the
official anchors. Standing lesson reinforced: **anchor every claim to the official numbers, not to
a local arm chosen for convenience.**

## What survives, and what the real target is

- Survives: **the realisation mechanism** (89.17%, and 32/32 silent targets rescued). That is real
  and reusable.
- Survives: **the renderer produces 961 calls/target with precision 0.5706** at B=885.
- Does NOT survive: the claim that this yields ~+0.05–0.22 composite.
- The actual official prize, if any, is the gap between official fidelity ~0.492 (chance) and the
  rendered 0.571 — about **+0.0437 composite** — and that requires the shipped package's DE table to
  be worse than Arm A's rendering on the same data. Unverified either way.

## Sweep status

Rendered budget sweep running (B = 50, 100, 150, 300 on CD4). It will still be informative about
*shape* (how fidelity/jaccard/nmae move with B), but per the above it cannot be read as leaderboard
headroom. Completion and results go to `gate_1025/sweep_1026/`.
