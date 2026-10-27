# DECISION 2026-10-25 01:05 ET — the amplitude lever is 4-of-6 positive; NMAE is the single blocker, and the break-even is specific

## The cleanest official experiment we have
`v8a` and `v8b` differ **only in amplitude** (OFFICIAL_SCORES.md: "v8b amplitude B=claim (candB)" vs
"v8a layout-only on exact v5b totals"; 10-19 lesson: "The amplitude change (candB vs candA) bought
MSE back"). Same totals, same layout. So every component delta is attributable to amplitude.

| component | v8a | v8b (×1.5) | Δ | contribution to overall (Δ/6) |
|---|---:|---:|---:|---:|
| pds | 0.680 | 0.708 | **+0.028** | +0.0047 |
| mse | 0.016 | 0.095 | **+0.079** | **+0.0132** |
| fid | −0.064 | −0.060 | +0.004 | +0.0007 |
| reach | 0.162 | 0.178 | +0.016 | +0.0027 |
| jac | −0.015 | −0.015 | 0.000 | 0.0000 |
| **nmae** | 0.102 | **−0.059** | **−0.161** | **−0.0268** |
| **overall** | 0.1470 | 0.1412 | **−0.0058** | (sum −0.0057) |

**Four of six components improve with more amplitude.** MSE gains +0.079 — the largest single
component move we have ever measured — and PDS gains +0.028, which matters because PDS is 71.6% of
the score. The arm loses only because NMAE collapses by −0.161, three times the MSE gain.

## The break-even, stated precisely
Holding the other five at their v8b values, amplitude ×1.5 breaks even when

```
|ΔNMAE| ≤ ΔMSE = 0.079        (i.e. NMAE ≥ +0.023 instead of −0.059)
```
and it beats v8a outright when `|ΔNMAE| < 0.079 + 0.028 + 0.020 = 0.127`.

So the question is **not** "is more amplitude good" — on five components it demonstrably is — but
**"why does NMAE fall 0.161 when amplitude rises 1.5×, and can that be decoupled?"**

## The mechanism, from the metric definition (metrics brief §7)
```
NMAE_p = Σ_{g∈S_p} |lfc_hat - lfc| / Σ_{g∈S_p} |lfc|
```
The gate `S_p` is the **reference's own significant set**, and a gene the submission does not carry
(or carries null/non-finite) is read as **predicted lfc = 0**. Therefore:

- The denominator is fixed by the reference — amplitude cannot make the task easier or harder.
- The numerator is minimised at the **correct per-gene amplitude**. Overshooting every gene by 1.5x
  adds `0.5·|lfc|` of error on every gene in the gate simultaneously, which is why the damage is so
  large and so uniform.
- Critically, NMAE is scored on the reference's significant set, **not** on the submission's calls. So
  a gene we call significant-but-tiny still enters with its (too small) lfc, and a gene we omit is
  scored as 0.

**This means NMAE and MSE want opposite things from a global amplitude, and the reason is scale, not
noise.** MSE is a `ratio_of_sums` over the whole panel and is dominated by high-expression genes;
NMAE is a per-gene relative error over the reference's DE set and is dominated by small-|lfc| genes.
A global ×1.5 helps the former and destroys the latter.

## The implied experiment (specific, and not one of the closed sweeps)
The closed "amplitude/rho sweep" moved a **global** scalar. The gap above says the fix is
**gene-dependent**: apply the larger amplitude only where it helps MSE (high-expression, high-|lfc|
genes) and keep or shrink it where NMAE is decided (the reference-DE set, especially small-|lfc|
genes). Concretely: fit an amplitude **per gene-expression decile**, holding the call set fixed, and
require in the same run that MSE rises while NMAE falls no faster than the break-even above.

This is distinguishable from every closed arm: it is not global amplitude (v8b), not amplitude vs rho
(the closed sweep), not B/claim, not call density, and not dispersion. It targets the one place where
the two best components disagree.

## Status
Recorded, not yet run. It needs an emitted + officially scored arm to be believed, and the local
harness **cannot** rank it (MSE is `ratio_of_sums`-invariant to layout but *is* responsive to bulk
amplitude, so the local MSE number is meaningful here — unlike the layout arms — but the 10x/1.5x
transfer ratios still apply). Any run must be paired with the incumbent in the same job and must
report all six components plus B/claim.

Safety entry `5NeeoWn96lfm23hgExcm` unchanged. No package change, no submission.
