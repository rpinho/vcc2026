**B0/ceiling: 12.3% for half-replicate agreement; 9.95% for the estimated full-data ideal-predictor ceiling.**

**Verification PASS:** 168 training targets; 196,483 treated cells, 38,176 unique controls. No 79-target held-out outcomes read. Five seeded obs-ID splits; one resumable expression pass.

Select bins on one half, evaluate the other, swap and average. Cosine: target mean; E: aggregate SSE/evaluation energy. ±: five-seed sample SD, not experimental-replicate uncertainty. B0 evaluates half-data truth; ratios divide means.

| Bin | n* | Replicate cos ± SD | B0 cos | Cos fraction | Replicate E ± SD | B0 E | E ratio |
|---|---:|---:|---:|---:|---:|---:|---:|
| All | 168 | .6486 ± .0006 | .0798 | 12.30% | .3296 ± .0036 | 1.5896 | 4.82 |
| Norm low | 56 | .3524 ± .0032 | .0174 | 4.93% | 1.1969 ± .0068 | 9.2530 | 7.73 |
| Norm middle | 56 | .7394 ± .0019 | .0751 | 10.15% | .4392 ± .0050 | 2.5418 | 5.79 |
| Norm high | 56 | .8540 ± .0021 | .1470 | 17.21% | .2868 ± .0051 | 1.2143 | 4.23 |
| Proxy <10 | 0 | — | — | — | — | — | — |
| Proxy 10–999 | 53.6 [53–54] | .3235 ± .0022 | .0159 | 4.92% | 1.2841 ± .0043 | 8.8094 | 6.86 |
| Proxy ≥1000 | 114.4 [114–115] | .8010 ± .0011 | .1097 | 13.70% | .3068 ± .0037 | 1.4174 | 4.62 |

n*: mean target membership [range]. Proxy ≠ Wilcoxon.

Spearman–Brown per target: R=2r/(1+r); report mean R and C=mean(√R), the ideal-predictor cosine ceiling. Parallel-half independent errors are assumed; nonlinear pooling/unequal halves make this approximate. No clipping or E correction. Full B0 cosine=.08413.

| Bin | Mean R±SD | C±SD | Full B0/C |
|---|---:|---:|---:|
| All | .7454±.0007 | .8459±.0009 | 9.95% |
| Norm low | .4839±.0033 | .6704±.0029 | 3.00% |
| Norm middle | .8385±.0016 | .9130±.0010 | 8.72% |
| Norm high | .9137±.0016 | .9542±.0010 | 16.00% |
| Proxy 10–999 | .4586±.0025 | .6540±.0024 | 2.88% |
| Proxy ≥1000 | .8798±.0009 | .9358±.0006 | 12.26% |

“Barely pointed the right way” is fair: noise correction leaves B0 near 10% of the estimated cosine ceiling.

Reproduce: `nice -n 10 python -B SCRIPT` for `stream_halves.py`, `analyze_halves.py`, `verify_halves.py`, sequentially. Verification (verification.json, kept locally); per-seed bins (kept locally); [methods](METHODS.md). [Credit/prior-art](PRIOR_ART.md): zymobyte012; local code only. Next: contextualize B0.
