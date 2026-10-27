**Verdict: transfers.** One arm worth building for D/E/F: **control_bounded** (exploratory; no official-score forecast).

Cells are **C/A/PDS**, where C is cross at plain-transfer centered spread and A=C/√Q; brackets give 2,000 paired target-bootstrap ΔC 95% CIs versus plain, rematching spread per draw. All main arms use **6,006 shared measured genes**. B/T is matched within each recipient.

|Recipient (n)|Plain|Pooled map|Control map|H1-direct|Pooled blend [CI]|Control blend [CI]|B/T|
|---|---|---|---|---|---|---|---|
|H1 (168)|0.2076/0.1882/0.8404|0.1851/0.1534/0.6365 [-0.0694,+0.0216]|0.1954/0.1604/0.6370 [-0.0601,+0.0322]|—|0.2155/0.1948/0.8388 [+0.0059,+0.0096]|0.2159/0.1951/0.8393 [+0.0065,+0.0098]|1.1729|
|K562 (6864)|0.0262/0.0568/0.5304|0.0190/0.0259/0.5130 [-0.0094,-0.0052]|0.0171/0.0296/0.5168 [-0.0109,-0.0073]|0.0095/0.0254/0.5118 [-0.0187,-0.0147]|0.0264/0.0562/0.5312 [+0.0002,+0.0003]|0.0265/0.0562/0.5304 [+0.0003,+0.0004]|0.1298|
|RPE1 (1996)|0.1712/0.3477/0.8313|0.2016/0.3940/0.7561 [+0.0248,+0.0358]|0.2056/0.3995/0.7548 [+0.0287,+0.0397]|0.1019/0.1828/0.6396 [-0.0755,-0.0628]|0.1770/0.3580/0.8315 [+0.0055,+0.0060]|0.1772/0.3584/0.8311 [+0.0057,+0.0063]|0.2203|
|HepG2 (1803)|0.2514/0.3700/0.7822|0.2799/0.3925/0.7281 [+0.0174,+0.0394]|0.2773/0.3907/0.7367 [+0.0150,+0.0369]|0.1162/0.1539/0.6494 [-0.1438,-0.1263]|0.2605/0.3804/0.7793 [+0.0084,+0.0097]|0.2603/0.3802/0.7796 [+0.0082,+0.0095]|0.3792|
|Jurkat (2116)|0.2090/0.2548/0.8413|0.1859/0.2003/0.7364 [-0.0288,-0.0174]|0.1848/0.1992/0.7375 [-0.0299,-0.0183]|0.1133/0.1282/0.6581 [-0.1026,-0.0885]|0.2107/0.2555/0.8391 [+0.0017,+0.0019]|0.2108/0.2555/0.8390 [+0.0017,+0.0020]|0.5510|

Original 18,077-gene H1 B0: C=0.115447, PDS=0.855895. On the shared axis at plain B, B0 C=0.196385, PDS=0.865163. B0+pooled_bounded: C=0.117968, ΔC CI [+0.001755,+0.003217], ΔPDS=-0.002816. B0+control_bounded: C=0.118116, ΔC CI [+0.001953,+0.003350], ΔPDS=-0.002281. These partial-gene updates preserve B0 outside the intersection; their weights use training recipients only. Its H1 B0 ΔPDS 95% CI is [-0.005204,+0.000428]; stricter PDS-CI noninferiority does not pass.

Cross-CI/PDS successes: pooled_map: 0/5; control_map: 0/5; pooled_bounded: 5/5; control_bounded: 5/5. A success requires paired ΔC CI lower>0 and mean ΔPDS≥−0.005. The same arm must pass in a majority and beat H1 B0. Full maps fail the PDS criterion; H1-direct fails everywhere. K562’s A falls despite positive ΔC. Per-arm ΔPDS intervals, centered-cross diagnostics and stricter PDS-CI flags: [all_metrics.csv](results/all_metrics.csv).

Each outer recipient is absent from labels AND every training source; inner recipient holdouts choose hyperparameters and bounded weights. Controls alone condition gene weights. The 79 held-out H1 identities are excluded from every bank; H1 uses training168 only. Pooled maps use a 256-component input approximation; retained energy and conditioning limits are in [METHODS.md](METHODS.md) and the retained fit diagnostics. The original operational STATUS.md is omitted from this export.

Verification: independent Cholesky/SVD ridge/RRR checks; exact dyadic cross, energy and centered-B certificates; scalar PDS checks; 35-arm bootstrap replay; input hashes and all 20 selectors/source-exclusion receipts audited independently. These certify local diagnostics only.

**Official MSE is not measurable locally. The final phase is blind.** No rendering, packaging, submission, upload, promotion or commit.

Reproduce sequentially: `nice -n 10 python -B analyze.py`, then `verify.py --replay`, `audit_design.py`, `write_report.py`, with the same Python/nice prefix. CPU only, one worker, BLAS=1; effective inherited nice20. Credits and public checks: [PRIOR_ART.md](PRIOR_ART.md).
