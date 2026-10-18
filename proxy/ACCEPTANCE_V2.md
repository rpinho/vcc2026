# Pre-registered proxy V2 acceptance — 2026-10-18

Frozen before rebuilding tasks or computing any new seven-arm calibration.

Arms: v1, v2b, v3b, v4a, v4b, v5a, v5b, without exclusions. Six official scaled components and official overall from OFFICIAL_SCORES.md. A rounded official value is not an exact measurement; ties (notably MSE zeros) stay tied.

Choose destination weights only from similarity of released control profiles, never from arm scores. Primary proxy is their fixed weighted macro. Fit each component with a one-dimensional affine proxy→official model on six arms; predict the seventh (leave one arm out). No sign inversion is allowed for a usable component: a negative training slope predicts the training mean and fails the component direction gate. Total prediction is the equal mean of the six held-out component predictions. Also report independently fitted overall for diagnosis, not selection.

Full proxy passes only if: (1) all seven arms are complete under the same task/scorer configuration; (2) held-out total Spearman >=0.80; (3) total mean absolute error <=0.03; (4) every non-tied official pairwise component direction among v4a/v5a/v5b is correct at official reporting precision (absolute difference >=0.001); (5) no component has negative held-out rank association. All tasks use the official metric/scaling functions and valid, positive-headroom anchors. Incomplete or non-identifiable components fail closed.

Component proxy: trust only components whose own held-out Spearman >=0.80, mean absolute error <=0.05, nonnegative full-fit slope, and the same v4a/v5a/v5b direction gate pass. Report all six regardless. Any constant or missing predictor fails. No total optimization using a subset is represented as a full score.

Report Kendall, per-arm held-out absolute error, raw proxy/official rank correlation, and sensitivity excluding v3b alongside (never replacing) the primary seven-arm result. No choice of destination, panel, scale, affine fit family, thresholds or aggregation based on final correlations. A pass is provisional with seven retrospective arms and requires next prospective official shot confirmation.

Task feasibility: never duplicate labeled cells to assert 400 independent real observations. Bootstrap replication does not create information or repair MSE/anchors. If no 300-target public destination has 400 independent cells per target, record the mismatch and report a separate 400-cell subset diagnostic; do not silently declare it matched.

Numerical convention (before any new seven-arm output): rank values are rounded to 12 decimal places solely to suppress floating-point tie breaks; this is far finer than official reporting precision. Fits and errors use unrounded values. Task/scorer provenance, five-split valid anchors and reference-identity checks gate any trust declaration.
