You are one of three scientific advisors. The operator has been told by the user: "make better models, try new ideas, research." Do not default to closeout. Be concrete.

Context:
- The final competition entry is frozen and cannot be replaced. This is research toward a better model, not a resubmission.
- Best official validation was 0.16438; leader 0.47203. The dominant gap is expression MSE.
- Public validation has had no new entries since Oct 22. Local historical harness cannot rank recipes reliably.
- A new public leave-one-cell-line-out benchmark exists (K562/RPE1/Jurkat, official metric profile) but two advisors already rejected spending 11 days on it. Do not revisit that.
- We have retained CPU data: 8 public Perturb-seq families, ~10,781 response rows, 300-target panel, existing response banks, and frozen model source. We have a 300-target donor coverage of 298/300 (ANG and DPY30 lack a target-specific row).
- Four untried mechanisms are specified: (1) response-reliability/source weighting, (2) combined STRING/sequence/gene-text features, (3) prior for missing-source targets, (4) saturation routing.
- The response model is a 5x-gene-input MLP predicting a residual on a pooled donor direction, blended at 0.3, then rendered to counts. MSE is dominated by the renderer/layout, not only direction.

Question: what is the single highest-value model experiment to build next, using retained data and at most ~2 CPU-days, that could produce a defensible local signal even without the official scorer? Specify: hypothesis, exact comparison/baseline, split/leakage rule, metric, and falsification threshold. Prefer a mechanism that could later be ported to the frozen recipe. If none is worth doing, say so plainly.

Answer in <=260 words, four labelled lines: VERDICT / EXPERIMENT / PROTOCOL / KILL TEST.
