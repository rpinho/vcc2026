**Covered-target transfer dominates; missing rows explain 0.5% of error. Exact n_conf stratification remains unavailable.**

Verified frozen inference bitwise; integer-rendered 168 training targets. No 79-target outcomes entered fitting/selection. E=Σ||b̂−b||²/Σ||b||²; cosine/scalar columns are mean/median. [Per-target CSV](results/per_target.csv).

| Stratum | n | Error mass % | E | cosine | scalar* |
|---|---:|---:|---:|---:|---:|
| All | 168 | 100.0 | 1.64 | 0.084 | 0.045 |
| Direct row | 161 | 99.5 | 1.64 | 0.088 | 0.049 |
| No row | 7 | 0.5 | 1.68 | -0.013 | -0.016 |
| Essential-screen member | 55 | 54.6 | 1.39 | 0.155 | 0.151 |
| Other | 113 | 45.4 | 2.09 | 0.050 | 0.013 |
| H1 KD <50% | 7 | 2.2 | 3.23 | 0.060 | 0.015 |
| KD 50–79% | 8 | 4.1 | 1.83 | 0.076 | 0.042 |
| KD ≥80% | 153 | 93.7 | 1.62 | 0.086 | 0.045 |
| Moment proxy <10 | 0 | 0 | — | — | — |
| Proxy 10–999 | 50 | 12.3 | 13.00 | 0.020 | 0.005 |
| Proxy ≥1000 | 118 | 87.7 | 1.46 | 0.111 | 0.089 |

Scalar*=unconstrained best multiplier. Aggregate error is **58.3% direction / 41.7% magnitude**; macro target-normalized error is 6.27. Native-comparator PDS=0.856 despite cosine≈0.084. Integer conversion changes SSE by only 0.003%. Biology/calibration is deficient; official dispersion/credit effects remain unresolved.

K562 donors cover 160 targets and overlap 97.9% of error; Jurkat 58/57.1%, RPE1 56/54.8%, HepG2 43/42.1%, HEK293T 41/24.0%, HCT116 40/23.5%, CD4 37/20.6%. These overlap. One-family targets carry 28.9% of error; four-family targets 46.3%.

“Essential” means K562 essential-screen membership, not measured H1 fitness. Exact Wilcoxon n_conf needs deleted cell ranks; all 168 exact values remain missing. The separately labeled moment/BH proxy is **not** Wilcoxon. Full source-KD strata and limitations: [METHODS.md](METHODS.md), [strata.csv](results/strata.csv).

Ranked fixes; gains reduce LOO SSE, **not official MSE**:

1. **Source-family residual transfer:** **10.5%** reduction (exploratory 95% CI 8.3–12.8%); PDS +0.00385; centered budget retains 80.4%. Fit <1 minute; rendering ≈66 seconds; existing data.
2. **Magnitude calibration:** global scalar≈0.133 yields **40.0%** (CI 34–47%), PDS essentially unchanged; <1 minute. Essentiality gives 40.27%, source-KD 39.97%, response-size regression 40.23%: little incremental biology. A full [TRADE](https://www.nature.com/articles/s41588-025-02169-3) adaptation is unmeasured, 2–6 engineering hours. Global shrink leaves **1.74% of budget**, so official credit may erase/reverse its distance gain. Context-similarity weighting alone gains 0.64% and loses 0.0066 PDS.
3. **Coverage/context acquisition:** [inventory](results/source_inventory.csv), no downloads. KOLF covers 164/168 (97.6%), newly covers 7 (4.2%); VIPerturb 130/168, new 7 (filtered: 68/new 4); full CD4 156/168, new 5 and 119 additional CD4 rows. Replogle essential/GW/RPE1 cover 53/160/56; Nadig Jurkat/HepG2 58/43; all add zero union coverage. HipSci ≥3/168 and LINCS ≥9/168, zero newly confirmed; restricted metadata prevents complete fractions. Perfectly fixing all missing rows buys ≤0.5%; new-context gain **unknown**, not bounded by that ceiling. Planning: CD4 2–8 h, KOLF 4–24 h, HipSci 2–6 h, VIPerturb 4–12 h, LINCS 1–3 days; unbenchmarked acquisition/harmonization, outside this job.
4. **SLIM/rank-10:** measured ≈2 seconds each; 38% distance reductions but PDS collapses to 0.556/0.505. Reject replacements.
5. **Bill retraining:** existing weights already underpin B0. New gain unknown. Timing projects **8.7 h** for nine fits, **3.0 h** for production only, excluding caches/scoring. Commands request nice10; macOS forces inherited nice20. Old 16 h is not current measured throughput; full plan still exceeds this job.

**Test first:** render the family-transfer correction. Require refitted target-CV SSE gain, PDS noninferiority (−0.005), and improved capped MSE after `vcc_mse_instr_1026` reports. PDS is locally measurable now; official MSE gain remains unverified. No promotion.

Reproduce: `nice -n 10 python -B reproduce.py`. Primary-source checks: [PRIOR_ART.md](PRIOR_ART.md); numerical verification: [receipt](results/verification.json).
