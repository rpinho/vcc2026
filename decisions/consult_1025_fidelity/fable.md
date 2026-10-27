**VERDICT:** Yes, read the metric's definition first, but cap it at about one hour. Your own numbers show why: (0.1767 − 0.5866) / (0.5996 − 0.5866) = −31.53 exactly. Fidelity is scaled by a denominator of only 0.013, because the trivial baseline is almost as good as replicates on this metric. So the score mostly punishes falling below about 0.587.

**EXPERIMENT:** A confidence-gated sparse delta. Output the trivial baseline for every gene. Add our predicted delta only on the top-k genes per perturbation, ranked by confidence (frozen-MLP |delta| × cross-family sign agreement in the response banks).

**PROTOCOL:** On the CD4 fold, using the official scorer, sweep k ∈ {0, 5, 10, 25, 50, 100, 300}. k=0 must reproduce baseline fidelity of 0.5866; if it does not, stop and fix the pipeline. Record all six components for each k.

**KILL TEST:** Kill the experiment if no k > 0 gives both:
- fidelity ≥ 0.58
- `score_avg` greater than the k=0 result

**Justification:** Because of the scaling, every 0.013 of fidelity is worth about one full scaled unit. Simply matching the baseline (k=0) lifts `score_avg` from −5.29 to roughly −0.03 using the other components' current values. No other lever comes close.

The two failed mechanisms fit this picture. Fidelity is a Wilcoxon-call yield, not a cosine. Shrinking amplitudes or changing directions still produces many weak significant calls with the wrong sign, so yield collapses.

Gating limits calls to genes where our direction is trustworthy. The baseline's fidelity is protected while we try to win reach and Jaccard. Before running the sweep, use the reading hour to confirm the yield denominator and how ties are handled.