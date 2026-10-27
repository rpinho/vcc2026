# DECISION 2026-10-25 04:06 ET — expression-leg programme CLOSED: neither arm replicates

## The two-recipient table (official anchors, positive = better)
| arm | H1 Δsum | CD4 Δsum | same sign both? | verdict |
|---|---:|---:|---|---|
| `no_mask` (remove DE mask) | **−0.1060** | **+0.0121** | **no** (opposite) | rejected |
| `no_pc` (remove 50-PC deflation) | **+0.0096** | **−0.0428** | **no** (opposite) | rejected |
| `unmasked` (both) | −0.0654 | pending | — | already rejected on H1 |

Under the pre-registered both-recipient rule (`decisions/20261025_0235`, `PREREGISTRATION.md`), **no arm
passes**. The sign of both interventions **reverses** between recipients, and in both cases the reversal is
larger than the original gain:

- Mask removal: helps CD4 by +0.0121, but costs H1 −0.1060 — a **9x larger** effect in the other direction.
- PC removal: helps H1 by +0.0099, but costs CD4 −0.033 (scaled components) / −0.006 overall.

This is not "small effects with noisy signs". It is a sign flip with a large magnitude, on the same code,
same layout, same scorer, same seeds. The mechanism is not understood well enough to predict which way it
will go on a withheld context, so it cannot be shipped.

## What the whole factorial established
1. The expression-leg mask is **not** suppressing recoverable signal. Removing it moves P50-direction
   quality in opposite directions on two public contexts. On H1, where the effect is large and well
   measured, it is clearly harmful.
2. Rebalancing the direction leg changes the profile's *shape*, and the frozen `B=claim` solve silently
   converts that into an *amplitude* change (`alpha` 0.444→0.261 H1, 1.546→0.795 CD4). So these arms are
   not clean tests of "more expression information"; they are coupled shape+amplitude interventions. That
   was a design limitation of the arm, and it is now recorded rather than glossed.
3. **The pre-registered rule did its job.** A single-recipient read would have promoted `no_pc`
   (+0.0096 on H1, verified, interval excluding zero, passes the priced NAE rule). The second recipient
   killed it. This is the third time this project's one-recipient local gains failed to replicate — the
   rule exists exactly for this.

## What is now closed
The expression mask, the 50-control-PC deflation, and their interaction are **closed as package levers**.
Do not re-run them. Do not build a recipient-conditional policy ("mask on for H1-like, off for CD4-like"):
we have no label-free way to identify which regime a blind context is in, and fitting one to two observed
recipients is precisely the self-deception both consults warned about.

## Next
Both consults (Kimi, Fable — `consult_1025/kimi_0402.md`, `fable_0402.md`) independently said the same thing:
stop response-transform sweeps and spend the remaining time on **verification of the retained final package**.
That work is now done for the highest-value question: the published artifact is identified and byte-exact
(`decisions/20261025_0405`), and its 900 per-group receipts have been independently re-audited (next entry).

No package change, no submission; safety entry `5NeeoWn96lfm23hgExcm` unchanged.
