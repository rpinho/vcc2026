**The shared response explains a modest part of B0’s failure.** Training168 only; no held-out access.

Projection onto each target-LOO shared mean captures **11.74% of truth energy** (a descriptive per-target amplitude oracle). The unscaled mean explains **4.44%** versus zero. B0’s mean cosine with the shared mean is **0.0418**; only **0.94%** of B0 energy lies along it. After removing that direction, B0’s residual has E=**1.759**, cosine=**0.077**. Nested scalar shrinkage (mean β=0.1128) gives orthogonal E=0.988: weak target-specific distance information.

| Arm | E | cosine | PDS | B/B0 % | ΔE [95% CI] | ΔPDS [95% CI] |
|---|---:|---:|---:|---:|---|---|
| B0 | 1.642 | 0.084 | 0.856 | 100.00 | — | — |
| (a) m | 0.956 | 0.193 | 0.476 | 0.00467 | -0.687 [-0.881, -0.535] | -0.3800 [-0.4288, -0.3320] |
| (b) m + βB0 | 0.943 | 0.208 | 0.523 | 1.46 | -0.699 [-0.890, -0.550] | -0.3329 [-0.3820, -0.2862] |
| (c) m + βB0⊥ | 0.945 | 0.208 | 0.521 | 1.25 | -0.697 [-0.889, -0.548] | -0.3348 [-0.3839, -0.2877] |
| (d) family correction | 1.468 | 0.094 | 0.860 | 80.23 | -0.174 [-0.225, -0.130] | 0.0036 [0.0005, 0.0074] |
| (e) c + family | 0.944 | 0.209 | 0.520 | 1.11 | -0.698 [-0.890, -0.549] | -0.3362 [-0.3853, -0.2895] |

Differences are versus B0; intervals use **2000 paired target bootstraps**, conditional on fitted predictions and the fixed 168-candidate panel. Lower E and higher PDS are better. All fitted means/coefficients exclude the evaluated target; stacked family fits also exclude each inner target from its shared mean and beta training. PDS exactly matches the preceding job’s native comparator, including panel-wide gene exclusion.

The table evaluates the literal vector equations, excluding each row’s own gene. The previous job’s feasibility projection is a separate [sensitivity table](results/projected_table.md); projected (d) reproduces E=1.46945, PDS=0.859745 and B/B0=80.35%. B is centered across-perturbation energy.

**Reading:** Adding the mean mainly shrinks distance error while collapsing retrieval and spread; family correction alone preserves both. The shared component does not explain most of the mismatch.

**Official MSE is not measurable: vcc_mse_instr_1026 found no instrument.** B constrains official MSE credit, so these distance gains cannot establish official gains.

Verified: nested fits, native PDS, bootstrap intervals and exact rational SSE. Reproduction, conventions and next step: [METHODS.md](METHODS.md). Credits: [PRIOR_ART.md](PRIOR_ART.md). CPU only; nice10 requested, inherited nice20; one BLAS thread. No rendering, submission, upload, packaging or promotion.
