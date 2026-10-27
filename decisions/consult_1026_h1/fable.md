VERDICT: Do not replace the expression head with the control mean — it costs 0.061/6 ≈ 0.010 composite on mse and risks collapsing pds (0.706, our only strong component). Let the manifold run finish but it is not the main bet; the DE-realization line is.
EXPERIMENT: Calibrated DE-call realization: predict a per-target call budget from the skewed reference distribution (median 31), fill it with top-ranked |lfc| genes (precision 0.752/0.705), assign sign only where |log2FC|>0.8 (accuracy 0.636), realize via the 89% mechanism as sparse nudges on called genes only.
PROTOCOL: Four arms on H1 fold with the official scorer: incumbent | calibrated budget + ranked genes | + sign gating | + sparse-nudge realization. Score all six components; log raw mse and pds for each arm.
KILL TEST: Best arm must beat incumbent 0.021555 on H1 composite by a margin exceeding fold noise (reliability 0.528 → demand replication on a second seed/fold), with raw mse not worse than the arm's incumbent 1.4097 and pds not below incumbent. Fail → the DE line is closed and we fall back to manifold-or-nothing.
CEILING: 0.21–0.23 realistic; 0.25 requires manifold passing. 0.47 is unreachable: leader mse 0.793 vs our 0.061 with delta cosine 0.005.

Arithmetic. Composite = sum(0.061, 0.706, 0.116, −0.066, 0.186, −0.016)/6 = 0.987/6 = 0.1645 ✓. Control-mean replacement: mse scaled → 0 (raw 1.0), composite −0.010; a constant prediction cannot discriminate perturbations, so pds 0.706 is also at risk (up to −0.118). Net negative. Rejected.

Manifold: pass gate requires cosine ≥0.20 and MSE ratio ≤0.90 in only 4/7 families; even a full pass moves raw mse from 0.973 toward ~0.93, roughly +0.1 scaled → +0.017 composite, discounted by pass probability. Modest.

DE realization: fidelity (−0.066) and reach (0.186) sit far below leader (0.564, 0.495). Capturing one quarter of that gap is (+0.158 + 0.077)/6 ≈ +0.039 composite — larger than the manifold’s best case — and rests on capabilities already measured (ranking precision 0.752, sign 0.636 above 0.8, 89% realization). Jaccard is compressed (leader 0.059) and is not a target.

Q2 answered above: if manifold fails, the single next experiment is the same DE-realization arm with the budget predictor as the variable.

Q3: No. A test target’s own replicate cells are exactly the response; any construction reaching cosine ~0.8 without the manifold working at 0.005 is using them. Falsifier: rerun with the held-out target’s cells removed; if cosine collapses, it was leakage.

Q4: 0.21–0.23. Hard truth: our expression model is dead, and that is 2/6 of the score.