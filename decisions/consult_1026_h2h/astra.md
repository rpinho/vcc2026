VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
OBJECTIVE:

If manifold loses, run **response-reliability source weighting in the shipped model**. Keep architecture, data, splits, training budget, and ensembling fixed; compare weighted and unweighted retrains. Estimate source reliability from independent replicate agreement in perturbation-minus-control responses, using training targets exclusively. Freeze weights before evaluating held-out targets. Weight reproducibility, not response magnitude.

**Kill it** if independent replication cannot support those weights, or if the paired target bootstrap does not establish lower expression MSE without worsening the official composite. No rescue sweep after seeing validation results.

If manifold wins, **run a frozen head-to-head through the official evaluation pipeline before building a submission package**. Freeze the winning arm and hyperparameters. Improvements in cosine alone are insufficient. Target leakage, transductive response fitting, inconsistent controls or gene alignment, failure on the official split, or DE losses that erase expression gains invalidate it. Keep the safety entry last; submit nothing outside the review window.

The realistic objective is **(a) move up from rank 463**. Nothing supplied supports forecasting 0.25. Stop treating pooled-branch victories as shipped-model victories, reopening the two negative lines, or launching a portfolio of weakly controlled experiments. The manifold result earns a decisive test, not a new narrative.

**No GPU rental now.** The next work is controlled weighting, evaluation, and uncertainty estimation. No supplied evidence identifies compute as the constraint. CPU is sufficient for this decision sequence; renting hardware would not resolve the central uncertainty.