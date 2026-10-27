# vcc_mse_instr_1026

**Verdict: instrument not found.** Every formulation fails C4’s sign and 0.3–3× ratio gate. No parameters were fitted. The 79 targets were used only to read contrasts.

Entries below are **local/official scaled-delta ratios**: positive means correct sign; zero means blind. The fixed ruler is the existing SEARCH anchor span, 0.95798309, without calibration to official outcomes.

| Variant | C4 +.0791 | C3 +.0663 | C6* −.0670 | C5 −.0219 | C1 −.0007 |
|---|---:|---:|---:|---:|---:|
| Current clamp | 0 | 0 | 0 | 0 | 0 |
| Capped, clamp lifted | -0.0153 | 0 | 0 | -0.0700 | 0 |
| Uncapped unbiased | -0.8147 | -0.9531 | +0.0543 | -0.3257 | +95.4708 |
| Trace cap lifted | -0.0153 | 0 | 0 | -0.0700 | 0 |
| Budget, labels | -0.0153 | 0 | 0 | -0.0700 | 0 |
| Budget, label-free | -0.0153 | 0 | 0 | -0.0700 | 0 |
| Plain bulk-delta | -0.8572 | 0 | 0 | -0.3127 | 0 |

*C6 is only the 694-versus-885 budget surrogate. Calib3 contains no genuine v8/v9 layout pair; its “v8a” uses v9 code. Thus even a favorable surrogate sign cannot validate C6.

**Leave-one-contrast-out:** no fitted mapping, so held-out signs equal direct signs. In C1/C3/C4/C5 order, uncapped unbiased gives ✓/✗/✗/✗; every other listed variant gives ✗/✗/✗/✗. Genuine C6 is unavailable; only uncapped unbiased matches its surrogate. Shared package legs make these checks dependent.

**Saturation:** the real-side jackknife cap fires on 0/79 targets in every analog. B binds native v5b and both v9 analogs, making layout invisible; 1.5× amplitude clears B but increases plug-in error more than earned credit. Raw MSE rises 0.985006→0.986167; lifting the score clamp exposes C4=−0.001212 scaled. The [0,1] clamp otherwise maps both to zero because the local baseline is 0.961342.

Label-free control resamples at n_p=400 leave claims unchanged. Per-package B/claim/ρ, exact scorer locations, formulas, provenance and reproduction commands are in [DETAILS.md](DETAILS.md) and [budgets.csv](budgets.csv).

**Verification PASS:** independent installed-kernel replay matches retained rows within 2e−12; separate budget and rational aggregation checks pass. All 79 native v5b bulk reconstructions match cached moments within 8.89e−16. Missing amplitude/billm totals were reconstructed without integer cells or layout rerenders.

**A125/A150/1.38×:** no qualified instrument, so no trustworthy official MSE prediction is issued. Next prerequisite: a faithful v8 layout replay on a representative panel. No submission, download, package, promotion, commit, or external write occurred.
