# DECISION 2026-10-24 22:00 ET — same-screen signal is strong; cross-line transfer is the bottleneck

## Correction to the 21:30 note
The 21:30 entry reported "H1 bank rho_opt 0.610, not leakage". Two things were wrong and are fixed here:

1. **Coordinate.** The audit scripts used bulk target sum 1e6. The frozen official preset uses
   **50,000** (`cell_eval2/configs/vcc2026.yaml: bulk_target_sum: 50000.0`; confirmed by
   `EvalConfig.from_preset('vcc2026').bulk_target_sum == 50000.0`). All numbers below are at 50,000.
   At the correct coordinate the shipped incumbent's oracle u is 0.976862, while the *measured*
   receipt is 0.983537. The two are different quantities (the measured one includes the sampling
   correction and the cap); the earlier "exact 5-decimal match at 1e6" was coincidence, not validation.
2. **Leakage.** Kimi was right that a 0.40 correlation does not establish cell disjointness. It does not.
   Direct barcode check: **2,570 of the 38,000 task cells also appear in `public_data/H1.h5ad`**, and
   **all 2,570 carry the same target label** as in the task (median 6.8% of each 400-cell group).
   So the shipped H1 bank was partially in-sample.

## The controlled test
`ceiling_audit_1024/clean_h1_v3.py` rebuilds the H1 bank from `public_data/H1.h5ad` twice with the
bank's own matched-control definition, verified against the stored `target_count_sums` (**max relative
difference 0.00e+00** for AKT2/ACAT2/DHX36):

| bank | cells | rho_opt (TS=50k, amplitude-matched) | cos_mean |
|---|---:|---:|---:|
| all cells (= shipped bank) | 20,000 | 0.7645 | 0.4701 |
| **minus every shared barcode** | 17,430 | **0.7323** | 0.4166 |
| shipped incumbent (receipt) | — | 0.1521 | — |

**Leakage accounts for 0.032 of a 0.612 advantage.** Removing it does not remove the effect.
Same-screen, disjoint-cell transfer of the target-specific response direction is ≈0.73.

## What the effect actually is
The H1 bank and H1_rich_search come from the **same experiment** (same cell line, same screen,
same guides, disjoint cells). That is not "same cell line, different screen" and it is certainly not
cross-line. Measured cross-line transfer on the same fold is 0.116–0.144 (HCT116, HEK293T, K562, CD4).

So the honest decomposition is:

| transfer regime | rho_opt | what it means |
|---|---:|---|
| same experiment, disjoint cells | **0.73** | target-specific direction is highly reproducible |
| different cell line | 0.12–0.14 | ~5x worse |

## Consequence for the final panel
The official D/E/F panel is a **cross-line** problem for every bank we hold *unless* D/E/F is a cell
line with a public screen covering the 300 targets. That is exactly what
`context_identity_1025/` is testing. If identification succeeds and coverage is good, the MSE lever
reopens at the same-screen level. If it fails, the final panel is a cross-line problem and the earlier
"data ceiling" conclusion stands *for that regime* — with the added caveat that it is a statement about
cross-line transfer, not about the information content of perturbation screens in general.

## Claims explicitly withdrawn
- "Same-line evidence is worth 4.8x" as a general statement about the final panel. Only the
  same-*experiment* number (0.73) is measured; same-line-different-screen is untested.
- The 21:30 note's "not leakage" conclusion, replaced by the barcode measurement above.
- Any absolute comparison mixing 1e6 and 50k coordinates.

## Still standing
- Amplitude is closed (optimal global rescale worth ≈0).
- The official panel overlaps our local primary fold (H1) by only 8/79 targets (10.1%);
  HCT116/HEK293T tasks overlap it at 29.7%/31.3%. Denominated against the official 300, the H1
  *bank* covers 11 targets = 3.7%. Both numbers are correct for their own denominators.
- The incumbent's official MSE (0.061) is measured, not inferred.
