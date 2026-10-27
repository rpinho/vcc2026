**vcc_align_1026 — verified training168 LOO; no held-out79 access.**

T=Σ||b||²; C=cross/T; Q=Σ||b̂||²/T; A=C/√Q. Every arm is rescaled to B0’s centered B. B below means B/T; uncentered Q is separate.

| Arm | C at B0 B | A | PDS | B/T | ΔC 95% CI | ΔPDS 95% CI |
|---|---:|---:|---:|---:|---|---|
| B0 | 0.1154 | 0.1235 | 0.8559 | 0.8161 | [+0.0000, +0.0000] | [+0.0000, +0.0000] |
| kNN | 0.4230 | 0.3193 | 0.5782 | 0.8161 | [+0.2424, +0.3738] | [-0.3287, -0.2293] |
| Ridge / RRR | 0.3984 | 0.3951 | 0.6619 | 0.8161 | [+0.2052, +0.3557] | [-0.2420, -0.1472] |
| B0 + kNN | 0.4286 | 0.3330 | 0.6055 | 0.8161 | [+0.2542, +0.3740] | [-0.2975, -0.2035] |
| B0 + map | 0.4035 | 0.4028 | 0.6833 | 0.8161 | [+0.2165, +0.3561] | [-0.2146, -0.1328] |
| Gene variance prior | 0.0956 | 0.1004 | 0.7500 | 0.8161 | [-0.0319, -0.0070] | [-0.1371, -0.0787] |
| B0 + kNN, PDS constraint | 0.1243 | 0.1329 | 0.8503 | 0.8161 | [+0.0064, +0.0116] | [-0.0128, +0.0019] |
| B0 + map, PDS constraint | 0.1356 | 0.1450 | 0.8503 | 0.8161 | [+0.0141, +0.0265] | [-0.0156, +0.0034] |
| B0 + kNN, PDS bound | 0.1186 | 0.1269 | 0.8539 | 0.8161 | [+0.0024, +0.0041] | [-0.0061, +0.0016] |
| B0 + map, PDS bound | 0.1235 | 0.1321 | 0.8556 | 0.8161 | [+0.0057, +0.0105] | [-0.0048, +0.0043] |

**Reading:** B0 + map, PDS bound adds aligned signal at matched spread: Δcross=+0.0080, ΔPDS=-0.0003; both bootstrap criteria pass (PDS narrowly). Q: 0.8732→0.8738; centered cross: 0.1054→0.1114.

Fits exclude the outer target; hyperparameters/weights use inner LOO. Follow-up blends constrain inner ΔPDS≥−0.005, using either its mean (“constraint”) or normal 95% lower bound (“bound”). These are exploratory. Similarities use non-H1 responses.

CIs: 2000 paired target bootstraps, fixed fits/candidate panel, spread re-matched per draw; exploratory, without multiplicity correction. PDS exactly follows the preceding comparator. Rescaling includes the mean; centered-cross diagnostics: results/analysis.json.

**One arm worth rendering:** B0 + map, PDS bound; analysis only here, no rendering performed. Modal choices: w=0.02, rank40, λ=0.1 (154/168 folds).

**Official MSE is not measurable locally.** These diagnostics establish no official-score gain.

Verified: independent primal fits/all outer predictions, three complete nested grids, native PDS, bootstrap replay, exact rational cross/energy/B. [Receipt](results/verification.json), [certificate](results/exact_alignment.json).

Reproduce: `analyze.py`, `constrained.py` (also `--conservative`), `verify.py`, `write_report.py`, each with `nice -n 10 python -B`. Definitions: [METHODS.md](METHODS.md); credits: [PRIOR_ART.md](PRIOR_ART.md). CPU only, one BLAS thread, inherited nice20. No rendering, packaging, submission, upload, promotion or commits.
