Output the four labelled lines FIRST. Justification after, max 200 words. Pick ONE experiment.

BREAKTHROUGH, measured tonight on the official scorer. This supersedes the earlier "direction is the disease"
diagnosis. Do not propose amplitude, donor reweighting, or another low-rank linear map.

THE METRIC, now reverse-engineered exactly from cell_eval2 source and verified to 1e-6 against the job's own
tidy table on all 79 CD4 targets:
    de_wilcoxon_direction_fidelity_yield_raw(t) = k / max(n_pred, n_conf)
  n_conf = # reference genes with p_adj < .05 (on-target removed); hidden, model cannot shrink it
  n_pred = # predicted genes with p_adj < .05
  k      = of those, how many the reference can adjudicate (real lfc non-null AND not exactly 0) and whose
           SIGN matches
  So it is exactly precision x min(1, coverage). Under-calling AND over-calling are both punished; calling
  nothing scores 0.0. This component is 99.4% of score_avg (-31.53 of -5.2871).

THE DIAGNOSIS (CD4 outer fold, 77 targets, exact per-target decomposition):
  mean n_pred = 189, MEDIAN n_pred = 1        mean n_conf = 255, median 30, MAX 3218
  32 of 77 targets (42%) call NOTHING -> fidelity forced to 0. They hold 9,899 of 19,658 total n_conf genes.
  mean precision where we DO call = 0.420;  on the 16 targets that call ~940 genes, precision = 0.585
  baseline anchor 0.5866, replicate anchor 0.5996, OURS 0.1790.
So we are not being out-signed (0.585 ~= baseline 0.5866). We are being out-CALLED.

TWO MECHANISMS, both in the renderer, neither in the model:
 (a) ELIGIBILITY GATE. layout.py:145:  eligible = keep & (confidence >= 1e6) & (direction == bulkdir)
     confidence = net*1e6 + zs, where net = (# donor banks significant & agreeing) - (# disagreeing).
     So confidence >= 1e6 requires AT LEAST ONE net donor bank to call that gene significant in the agreed
     direction. On cross-cell-line transfer (per-gene cosine 0.013 on CD4) often ZERO genes clear it, so the
     target calls nothing. Median eligible genes per target = 1.
 (b) BUDGET. nconf = exp(mean(log(source median calls))) = ONE scalar (~940) for every target, while true
     n_conf spans 1 .. 3218. Wrong at both extremes.

MEASURED HEADROOM (leak-free; no n_conf used to pick any parameter):
  actual                                                                  0.1790
  keep budget B=940, only stop the gate zeroing targets                   0.6111   (+0.4321)
  per-target budget = n_conf (oracle, needs hidden reference)             0.6985   (+0.5195)
  65/77 targets improve, 8 worsen. 0.6111 would EXCEED both anchors (0.5866 / 0.5996), flipping the 99.4%
  component from -31.53 to positive.

RESOURCES: CPU only, 8 public families, 300-target panel, frozen MLP, working official scorer + exact metric
reimplementation. Two cheap mechanisms already rejected on the official metric (source reweighting; global
amplitude calibration: MSE -23% but fidelity -29% because shrink pushed genes BELOW the significance bar and
destroyed coverage).

QUESTION: what single experiment most directly raises fidelity_yield_raw given the gate and budget are the
binding constraints? Be specific about how to widen eligibility without introducing leakage (the held-out
family's treated cells and its true DE table are both forbidden), and how to set a per-target budget without
seeing the hidden n_conf.

VERDICT:
EXPERIMENT:
PROTOCOL:
KILL TEST:
