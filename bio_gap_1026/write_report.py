import runtime as rt
rt.enforce()
import pandas as pd,numpy as np
A=rt.json.loads((rt.ROOT/'results/analysis.json').read_text());V=rt.json.loads((rt.ROOT/'results/verification.json').read_text());B=rt.json.loads((rt.ROOT/'results/benchmark_projection.json').read_text());S=pd.read_csv(rt.ROOT/'results/strata.csv');I=rt.json.loads((rt.ROOT/'results/source_inventory.json').read_text());C=rt.json.loads((rt.ROOT/'results/cpu_baselines.json').read_text())
assert V['status']=='PASS' and A['baseline_integer_emitted']
def row(kind,label,shown):
    v=S[(S.kind==kind)&(S.label==label)].iloc[0]
    if not v.n:return f'| {shown} | 0 | 0 | — | — | — |'
    return f'| {shown} | {int(v.n)} | {v.error_mass_pct:.1f} | {v.error_ratio_zero:.2f} | {v.cosine_mean:.3f} | {v.best_scalar_median:.3f} |'
table='| Stratum | n | Error mass % | E | cosine | scalar* |\n|---|---:|---:|---:|---:|---:|\n'
for kind,label,shown in [('all','all','All'),('coverage','direct','Direct row'),('coverage','none','No row'),('essential_proxy','screen_member','Essential-screen member'),('essential_proxy','other','Other'),('observed_H1_knockdown','<50%','H1 KD <50%'),('observed_H1_knockdown','50–79%','KD 50–79%'),('observed_H1_knockdown','>=80%','KD ≥80%'),('n_conf_MOMENT_PROXY','weak <10','Moment proxy <10'),('n_conf_MOMENT_PROXY','mid 10–999','Proxy 10–999'),('n_conf_MOMENT_PROXY','strong >=1000','Proxy ≥1000')]:table+=row(kind,label,shown)+'\n'
energy=A['summary']['LOO_family_transfer']['centered_energy']/A['summary']['B0']['centered_energy']
report=f'''**Covered-target transfer dominates; missing rows explain 0.5% of error. Exact n_conf stratification remains unavailable.**

Verified frozen inference bitwise; integer-rendered 168 training targets. No 79-target outcomes entered fitting/selection. E=Σ||b̂−b||²/Σ||b||²; cosine/scalar columns are mean/median. [Per-target CSV](results/per_target.csv).

{table}
Scalar*=unconstrained best multiplier. Aggregate error is **58.3% direction / 41.7% magnitude**; macro target-normalized error is {A['summary']['B0']['macro_relative_error']:.2f}. Native-comparator PDS={A['summary']['B0']['PDS168']:.3f} despite cosine≈0.084. Integer conversion changes SSE by only {100*V['emitter_relative_SSE_change']:.3f}%. Biology/calibration is deficient; official dispersion/credit effects remain unresolved.

K562 donors cover 160 targets and overlap 97.9% of error; Jurkat 58/57.1%, RPE1 56/54.8%, HepG2 43/42.1%, HEK293T 41/24.0%, HCT116 40/23.5%, CD4 37/20.6%. These overlap. One-family targets carry 28.9% of error; four-family targets 46.3%.

“Essential” means K562 essential-screen membership, not measured H1 fitness. Exact Wilcoxon n_conf needs deleted cell ranks; all 168 exact values remain missing. The separately labeled moment/BH proxy is **not** Wilcoxon. Full source-KD strata and limitations: [METHODS.md](METHODS.md), [strata.csv](results/strata.csv).

Ranked fixes; gains reduce LOO SSE, **not official MSE**:

1. **Source-family residual transfer:** **10.5%** reduction (exploratory 95% CI 8.3–12.8%); PDS +0.00385; centered budget retains {100*energy:.1f}%. Fit <1 minute; rendering ≈66 seconds; existing data.
2. **Magnitude calibration:** global scalar≈0.133 yields **40.0%** (CI 34–47%), PDS essentially unchanged; <1 minute. Essentiality gives 40.27%, source-KD 39.97%, response-size regression 40.23%: little incremental biology. A full [TRADE](https://www.nature.com/articles/s41588-025-02169-3) adaptation is unmeasured, 2–6 engineering hours. Global shrink leaves **1.74% of budget**, so official credit may erase/reverse its distance gain. Context-similarity weighting alone gains 0.64% and loses 0.0066 PDS.
3. **Coverage/context acquisition:** [inventory](results/source_inventory.csv), no downloads. KOLF covers 164/168 (97.6%), newly covers 7 (4.2%); VIPerturb 130/168, new 7 (filtered: 68/new 4); full CD4 156/168, new 5 and 119 additional CD4 rows. Replogle essential/GW/RPE1 cover 53/160/56; Nadig Jurkat/HepG2 58/43; all add zero union coverage. HipSci ≥3/168 and LINCS ≥9/168, zero newly confirmed; restricted metadata prevents complete fractions. Perfectly fixing all missing rows buys ≤0.5%; new-context gain **unknown**, not bounded by that ceiling. Planning: CD4 2–8 h, KOLF 4–24 h, HipSci 2–6 h, VIPerturb 4–12 h, LINCS 1–3 days; unbenchmarked acquisition/harmonization, outside this job.
4. **SLIM/rank-10:** measured ≈2 seconds each; 38% distance reductions but PDS collapses to 0.556/0.505. Reject replacements.
5. **Bill retraining:** existing weights already underpin B0. New gain unknown. Timing projects **{B['nine_fits_two_workers_hours']:.1f} h** for nine fits, **{B['three_seeds_two_workers_hours']:.1f} h** for production only, excluding caches/scoring. Commands request nice10; macOS forces inherited nice20. Old 16 h is not current measured throughput; full plan still exceeds this job.

**Test first:** render the family-transfer correction. Require refitted target-CV SSE gain, PDS noninferiority (−0.005), and improved capped MSE after `vcc_mse_instr_1026` reports. PDS is locally measurable now; official MSE gain remains unverified. No promotion.

Reproduce: `nice -n 10 python -B reproduce.py`. Primary-source checks: [PRIOR_ART.md](PRIOR_ART.md); numerical verification: [receipt](results/verification.json).
'''
words=len(report.split());assert words<=600,words
(rt.ROOT/'REPORT.md').write_text(report)
rt.save('results/report_check.json',dict(words_whitespace=words,maximum=600,verification=V['status'],exact_n_conf_complete=False,unresolved='Exact n_conf and H1-specific essentiality unavailable; official MSE pending instrumentation',baseline_integer_emitted=True))
rt.status(f'REPORT.md written ({words} words); exact n_conf limitation explicit; no official MSE improvement claimed.')
