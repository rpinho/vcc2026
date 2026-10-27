# Scope and reproduction

Training168 only; no held-out79 input/outcome files are opened. All writes stay in align/. Existing ../results/baseline.npz and ../shared_resp/results/sufficient.npz are reused; frozen networks, features and count renderers are not rerun. No package, download, GPU, submission, upload, promotion, commit or forbidden artifact access. Numerical processes have a filesystem/network audit hook, nice10 requested (the session inherits effective nice20), one BLAS thread and one worker. Free disk is checked before heavy steps (>=31.5 GiB), with an abort below 30 GiB.

From this directory:

```sh
nice -n 10 python -B analyze.py > analyze.log 2>&1
nice -n 10 python -B constrained.py > constrained.log 2>&1
nice -n 10 python -B constrained.py --conservative > conservative.log 2>&1
nice -n 10 python -B verify.py > verify.log 2>&1
nice -n 10 python -B write_report.py
```

# Coordinates and source-only similarity

Y and integer-emitted B0 P are the preceding jobs' 168 x 18077 log1p(CP50k) deltas. Each row's own target coordinate is zeroed before fitting/evaluation, exactly as before. The saved YY, PY, PP and baseline retrieval similarities are replayed against these arrays. PDS separately excludes all 168 target genes from every vector, with fixed candidate panel, n-1 rank normalization and midrank ties, as before.

The saved similarity matrices are PREDICTION-to-H1-TRUTH comparisons, not source-context similarities. They are never used as inputs. The source kernel uses the existing model_billmB/data manifest and its twelve non-H1 banks, selecting only rows with training168 target identities; H1 banks are skipped entirely. Within each of seven source families, existing measured source responses are averaged, with missing coordinates zero. The 168 target-gene coordinates are excluded from every source vector to avoid direct-knockdown signatures. The seven family vectors are concatenated conceptually, with absent family blocks zero. The Gram matrix is summed across blocks and row-normalized to obtain cosine K. This is a fixed per-row transform, without cohort fitting or H1 labels. Bank hashes must equal the existing manifest. Seven targets lack a source row.

# Nested LOO arms

For outer held-out target t, S has 167 targets. For every j in S, inner training is S without j (166). Source normalizations are per-row and need no label fitting. All label-dependent fitting, rank fitting, hyperparameter scoring and blend scaling exclude t; inner predictions also exclude j. The final 168-target spread normalization is a descriptive evaluation functional, not an inference-time fitted parameter.

**kNN:** k in {3,5,10,20}. Select the k largest source cosine similarities among training targets; clip negative weights to zero and normalize their sum to one. With no positive weights, use the uniform mean of all available training targets. Predict the weighted mean of their own-gene-masked true H1 deltas and then mask the predicted target's own coordinate. Stable target-index ordering resolves exact neighbor ties.

**Map:** through-origin ridge / ridge reduced-rank regression from the unit source features to H1 deltas. Objective is ||Y-XW||_F^2 + lambda ||W||_F^2 subject to rank(W)<=r. Ranks {5,10,20,40,80,unrestricted}; relative penalties {.01,.1,1,10,100} multiplied by mean training K diagonal. The unrestricted option is ordinary ridge. For K=U diag(d) U', use the eigensystem of D U'YY U D, D=sqrt(d/(d+lambda)), to obtain the exact ridge reduced-rank solution. Tiny numerical zero eigenvalues of source K (<max(1e-12,1e-11*dmax)) are omitted. Output rank is fit from TRAINING truth only; there is no output PCA fit on j or t. A zero source vector maps to zero. No intercept or target-identity features.

For each candidate, collect all 167 inner held-out predictions. Compute their centered spread B and the corresponding B0 spread on S. Multiply the candidate by sqrt(B0_S/B_candidate) and score its total cross against Y_S. Choose the candidate maximizing this matched-spread cross. This scale, learned from inner OOF predictions, is applied to the final outer prediction, avoiding arbitrary fold-specific shrinkage.

**Blends:** q=(1-w)P+w*scaled_source_prediction, w in {0,.05,.1,.2,.35,.5,.75,1}. Jointly select w and the source model's k or (rank,lambda) from the same strictly inner OOF candidate predictions. Apply a further inner-learned scalar sqrt(B0_S/B_blend) before scoring. Thus source hyperparameters and blend weights are all nested without needing to train a base component on the outer held-out truth. The selection objective is cross alone; PDS is an independent final acceptance criterion, not a hidden tuning filter. The w=0 option reproduces B0.

**Gene prior:** for each t, compute unbiased per-gene variance of true H1 deltas across the other 167 targets. Multiply P_t coordinatewise by this variance, then scale to P_t's own Euclidean norm. There is no tuned exponent, floor or held-out-label input. Its final global matched-B normalization is the same descriptive normalization used for every arm.

**Bounded follow-up, PDS-constrained blends:** the original cross-only selectors chose large source contributions and failed PDS. Before evaluating follow-up outputs, add two same-family blend variants with weights {0,.005,.01,.02,.03,.05,.075,.1,.15,.2,.35,.5,.75,1}. Jointly select component hyperparameters and weight by inner matched-B cross, subject to mean inner delta PDS >=-.005. PDS is computed only against the 167 outer training candidates (denominator 166), never against outer t. The same fixed training168 target-gene exclusion is used. A finite-candidate PDS constraint is not guaranteed to generalize to the outer target; both point and CI noninferiority checks are therefore reported. These additional arms are exploratory and were motivated by the original six-arm result, not independently preregistered primary hypotheses.

**Final bounded follow-up:** those point-constrained blends missed the final PDS point limit (delta=-.005596). Keep all earlier arms and add the same two selectors with inner lower normal 95% bound deltaPDS - 1.96*SD(per-target paired delta)/sqrt(167) >=-.005. The search objective remains matched-B cross, with the identical finite source-model/weight grids. This is a stronger inner constraint; the final acceptance intervals remain the requested 2000 paired percentile bootstraps. The normal bound is only a model-selection guard, not substituted for final bootstrap evidence. This final follow-up was specified before its outer outputs; no further search follows it. All four follow-up arms are exploratory, motivated by earlier outer results, so successful outcomes need independent confirmation.

# Metrics, spread matching and intervals

T=sum_i ||Y_i||^2 is the common truth denominator. For each full set of LOO predictions q, centered spread Braw=sum_i ||q_i-mean(q)||^2. Set alpha=sqrt(B0/Braw) and report C=sum_i alpha*q_i.Y_i/T, Q=sum_i ||alpha*q_i||^2/T, A=C/sqrt(Q), and B=B0/T. A is the stacked-vector cosine, i.e. aligned signal per unit total prediction spread, and is invariant to positive scalar rescaling. Centered B and uncentered Q are distinct; both are retained. The whole delta is rescaled, including its mean. To expose possible common-response amplification, also retain centered cross=sum_i (alpha*q_i-alpha*mean(q)).(Y_i-mean(Y))/T. No feasibility projection/count rendering is applied; these are literal vector-space diagnostics.

Primary acceptance: paired 95% delta-C CI lower bound >0 and mean delta-PDS >=-0.005. Also report the more conservative check that the lower delta-PDS CI is >=-0.005. No official-score forecast follows: official MSE is not measurable locally.

2000 paired target bootstraps, seed 20261026, percentile CIs. Predictions/fits and the 168-candidate PDS panel remain fixed, as before. For each resample, recompute truth denominator and both arms' centered B, then rematch spread within that replicate before calculating delta cross. This bootstraps the actual matched-spread statistic rather than treating its estimated normalization as fixed. PDS uses the mean resampled per-target rank difference. Intervals condition on fitted predictions, are exploratory and are not multiplicity-adjusted or evidence for a new destination context. A separate fixed-scale sensitivity is retained by the verifier.

# Independent verification

verify.py does not import analyze.py. It rehashes inputs, rebuilds source similarity from existing banks, independently replays outer fits, checks nested candidate optima and selected inner fits, verifies all metrics/native PDS and paired intervals, and certifies exact dyadic-rational cross/energy/spread for every primary output. Mathematical verification establishes the numerical diagnostic, not an official or biological score claim. Detailed receipts record its actual coverage.

Prior methods and literature: PRIOR_ART.md. Small kernels, coefficients, per-target summaries and certificates are retained; full challenger vectors/cell matrices are not saved.
