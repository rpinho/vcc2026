# DECISION 2026-10-25 04:10 ET — family_quality anchor, and why the bounded test stays

## Question
The archived `family_quality` direction was built against the **post-blend** incumbent profile and emitted
without another blend. The current pipeline calibrates a pre-blend direction, blends the learned profile at
0.3, then lays out. Rebuilding the old anchor would double-apply the learned profile; using a pre-blend
anchor changes fallback rows and per-target norms slightly.

## Consulted
Kimi (`orq best`, 2026-10-25 04:09). Verbatim conclusion: use the **pre-blend calibrated direction** in A1,
apply the unchanged family fusion there, and then give **both** arms exactly one learned blend at 0.3 and the
same v9 layout. Treat the archived +0.0132 as motivation, not an expected result, and verify budget matching
on the final scored candidates because downstream transformation can break calibration matching.

Kimi also answered the priority question **yes**: spend ~70 minutes on the frozen official-score test. The
MSE-only rejection evaluated a different acceptance criterion, and CD4's `-0.0019` is a source-space rho
proxy, not an official-score result in either direction. If H1 fails `C`/`T`, close the lever and move to
package reproducibility without another anchor variant.

## Decision
- **A1 anchor:** `old = logbulk(calibrated_preblend) - Lc`; family means and the frozen
  `sqrt(min(n_cells,1000))` weights; incumbent fallback on unmeasured entries; per-target norm restored
  exactly as `mse_direction_1024/evidence.py:profile()` does (own-target genes excluded from the norm).
- **A0:** the same calibrated pre-blend direction, unchanged.
- **Both:** one learned blend at 0.3, same layout/confidence/calls/scorer/anchors.
- **Gate:** `C>0`, `T>0`, paired target-bootstrap lower bound for `C` above zero, gates, and
  `|ratio_requested(A1)-ratio_requested(A0)|<=0.02`. No MSE-only hurdle.
- **Stopping rule:** H1 failure closes the lever. H1 pass triggers the identical two arms on CD4; a CD4
  failure also closes it. No third anchor variant, no tuning round.

## Implementation note
The `family_quality` fusion must be copied from `mse_direction_1024/evidence.py:source_candidates`; it is not
re-invented. The anchor is the only deliberate difference from the archive, and the report must say so.

No submission, no package change. Safety entry `5NeeoWn96lfm23hgExcm` and package
`ce78c81c453f74676ea0e6d2e1298cccf58aa0c2d6649f930b428b62a4939dc4` are unchanged.
