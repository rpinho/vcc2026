# calib3_1026 — local harness vs official board, component-level calibration (analysis only)

Written by the orchestrator from the calibration agent's report, because the harness refused report files from
subagents. Script: `calib3.py`; data: `contrasts.json`. Nothing was scored, rendered or submitted.

**Panels:** h2h2 and final_push H1_rich_search use the same 79 targets, the same reference and the same anchors.
Re-scaling with h2h2's anchors reproduces every stored final_push value to within 1e-9. Only the biology differs.
B0 sits in a different MSE regime: raw 1.41 against 0.985. In final_push's regime the local capped MSE doesn't move
with layout at all; in B0's it moves 0.6%.

## 1. Contrast table (scaled deltas, official vs local)

| Contrast (official vs local analog) | | pds | mse | nmae | fid | reach | jac | overall |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| C1: v5b−v9a vs v5b−v9a_885 (layout: native vs v9) | off | 0 | −0.0007 | +0.0008 | +0.0370 | −0.1040 | +0.0156 | −0.0086 |
| | loc | 0 | 0 (blind) | +0.0003 | +0.3432 | −0.0636 | +0.1132 | +0.0655 |
| | ratio | – | – | 0.34 | 9.3 | 0.61 | 7.3 | sign wrong |
| C1b: same pair, 694 surrogate | ratio | – | – | 0.45 | 11.0 | 0.62 | 8.2 | sign wrong |
| C3: v5b−v8a vs v5b−v9a_694 (official side includes v8's MSE loss) | off | 0 | +0.0663 | +0.0062 | +0.0333 | −0.0852 | +0.0143 | +0.0058 |
| | ratio | – | blind | 0.06 | 12.2 | 0.75 | 9.0 | sign right, but only because of MSE |
| C4: v8b−v8a vs v8b_like−v9a_694 (amplitude ~1.5×) | off | +0.0277 | +0.0791 | −0.1613 | +0.0038 | +0.0163 | −0.0005 | −0.0058 |
| | loc | +0.1010 | 0 | +0.0069 | +0.0517 | −0.0015 | +0.0159 | +0.0290 |
| C5: billmA−v9a vs billm-LOCO v9 885 − v9a_885 (biology) | off | +0.0257 | −0.0219 | +0.0083 | +0.0015 | +0.0055 | −0.0002 | +0.0030 |
| | loc | +0.0434 | 0 | +0.0035 | −0.0031 | +0.0067 | 0.0000 | +0.0084 |
| C6: v8a−v9a (v8 vs v9 layout code) | off | 0 | −0.0670 | −0.0053 | +0.0037 | −0.0188 | +0.0013 | −0.0144 (no local analog) |

- **Layout axis (C1, C3):** all DE signs agree. Fidelity is inflated 9–12× locally, jaccard 7–9×, and reach is
  understated at 0.6–0.75×. These ratios hold steady across C1, C1b and C3. Predicting C3's DE points from C1's
  ratios gives −0.0072 against −0.0063 official, and the reverse gives −0.0073 against −0.0086. That check is
  weak: C1 and C3 share v5b and their local legs are near-duplicates.
- **Leave-one-contrast-out over C1, C3, C4, C5:** the mapped composite gets the held-out official sign right on
  1 of 4 (2 of 4 without MSE); the raw local composite gets 2 of 4. The misses are MSE (C3, C5), the amplitude NMAE
  overshoot (C4) and the PDS ratio (C5). On the call/layout axis there is effectively one official contrast.
- **Correction to the R7 correction table:** official MSE for v5b−v9a is −0.0007, not +0.001. v9a is slightly
  better (raw 0.90912 vs 0.90982).
- **Missing pieces:** v9b has an official overall only and no local replay. v7a/v7b were never scored.

## 2. Projections vs B0 (official composite points)

**MSE assumption:** the official scaled MSE moves 0.811 per raw unit (from the exact v5b and v8a rows), and
billmA's official raw MSE is about 0.936. Half of D1770's local −0.58% raw change is assumed to carry over, giving
+0.0004 (range −0.0004 to +0.0007). The reason for halving: the 885-call v9 layout moved official raw MSE only
−0.08% against native, while the same contrast moved local MSE by exactly 0.

| Arm | Local Δ | Ratio model (DE+nmae) | Stratum model (DE) | MSE | Projected official Δ |
|---|---:|---:|---:|---:|---:|
| D1770 | +0.0345 | +0.0048 to +0.0061 | +0.0026 [5–95%: −0.0009, +0.0077] | +0.0004 | +0.003 to +0.0065 with reach; 0.000 to +0.004 without |
| B0s1 / B0s2 (seed noise) | +0.0008 / +0.0005 | +0.0011 to +0.0019 | +0.0030 / +0.0025 | 0 | noise floor about +0.002 to +0.003 |
| U885 | −0.0054 | −0.0013 | −0.0040 (78% chance ≤0) | +0.0001 | negative, about −0.001 to −0.004 |
| A125 / A150 | +0.0003 / −0.0002 | DE +0.002 / +0.007, but NMAE −0.008 to −0.06 / −0.02 to −0.18 | ≈0 / +0.006 | about −0.016 / −0.035 | negative; the harness can't see the MSE and NMAE terms that dominate (see C4) |

- **Ratio model:** maps D1770's local fid, reach and jac through the C1/C1b/C3/pooled ratios, giving official
  scaled fid +0.012 to +0.016, reach +0.011 to +0.014, jac +0.005 to +0.007.
- **Stratum model:** a mechanism-based model that splits targets by n_conf (number of truly DE genes). Fitting C1's
  exact official deltas puts the official panel at 31.5% weak (<10), 65% mid (10–999) and 3.5% strong (≥1000),
  against 15/52/33% locally. Public n_conf medians of 40–260 (`model_v8/NOTES.md`) point the same way, and the
  model reproduces C3 in sign and rough size.

**D1770's per-target pattern:**
- **Fidelity and jaccard follow the native pattern.** 112% of the fidelity gain and 103% of the jaccard gain come
  from the 26 strong targets (native: 98% and 104%). Per-target fidelity changes correlate 0.83 with native's. The
  top 5 targets give 38% of the gain (native 27%).
- **The mid targets differ.** D1770 loses fidelity on the 41 mid targets (−0.008), where native was flat. That is
  why the stratum model gives D1770's official fidelity a small negative value (−0.0055 scaled).
- **Reach moves the opposite way to native, but it is noise.** 92% of the reach gain is three extra true genes in
  three weak targets (ANTXR1, NCK2, TWF2). D1770's local reach gain (+0.0085 scaled) equals seed replicate B0s1's
  (+0.0083).

## Verdict

The calibrated projection is indeterminate, leaning slightly positive. It is not a confident official win, and it
is not an official loss either. R7's reasoning ("same axis as native, so same official loss") fails at component
level: native lost officially through reach (−0.104 scaled), and D1770 doesn't lose reach. But the positive
projection rests on two weak supports. The reach term is the same size as a seed replicate's, so it is noise. The
fidelity and jaccard terms give +0.003 to +0.004 under a single global ratio, but about 0 under the stratum model,
because on a weak-heavy panel like the one the official v5b/v9a deltas imply, mid-target losses cancel the
strong-target gains. Without the noise-level reach, D1770 projects to roughly 0.000 to +0.004 official (central
about +0.002, MSE included). That is at or below the +0.002 to +0.003 the same method assigns to pure seed noise.
The calibration rests on one official layout contrast, and its leave-one-out record off that axis is 1–2 of 4.
About 60% that D1770 is officially positive; no size should be quoted. The hidden test panel's mix of strong
targets decides the sign (break-even is about 4–5% strong). Nothing here argues for overturning "B0 final", but
nothing supports calling D1770 a known loss either.

## Proposed local scoring (not run)

One arm, about 20 minutes (CPU, `nice -n 10`, one worker): D1770 on the CD4 recipient.
- **Setup:** re-run the `layout_cd4_1025` harness with the frozen CD4 incumbent
  (`layout_cd4_1025/artifacts/frozen_head.npz`, same confidence, seeds and amplitude 1.5457), changing only the v9
  budget from 885 to 1770. The existing CD4 incumbent arm is the paired comparator.
- **Why CD4:** its stratum mix (29/63/8%) is close to the fitted official mix, so it tests the stratum model directly.
- **Decision rule:** convert the raw fidelity and jaccard changes with the official slopes (3.30 and 2.88). If the
  result is ≤0, close D1770 as officially null. If it is clearly above the seed floor, D1770 moves only to "plausible".
- **Caveat:** at 885, CD4 already selects only 243–885 columns (median 313), so 1770 may change few targets. That
  outcome would itself support "about 0 on weak panels".
- **Not needed:** a D1770 seed replicate, since reach is already treated as noise.
