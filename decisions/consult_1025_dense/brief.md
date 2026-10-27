# Independent review request — VCC dense-bulk / learned-norm ablation

You are reviewing a pre-registered experiment BEFORE its numbers are known. Do not run code. Give
a short, opinionated verdict with a ranked list. Answer in the same language as this brief (English),
and keep it under ~700 words.

## Situation

Virtual Cell Challenge 2026. Scored on 6 official components (PDS cosine, capped normalized MSE, DE
log2FC NMAE, DE direction fidelity, DE reach, DE Jaccard). Our best official validation entry is
billmA at 0.16438 overall (rank 463). Our incumbent package for the blind final panel is v5b amplitude
+ v9 emitter layout; blind score hidden.

Local harness: frozen official cell-eval2 scorer on a public 79-target H1 leave-one-family-out SEARCH
panel (400 cells/target). Known defect: this harness does NOT rank historical recipes correctly
(final_push/harness_audit_20261023/audit.json: pairwise agreement 9/10 for 5 recipes, 18/21 for 7
macro scores, 3/6 for layout/amplitude surrogates). So no local delta is certified to transfer.

## What is already measured (same harness, previous job)

| arm | overall | PDS | MSE | NMAE | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|---:|
| incumbent | -0.019164 | 0.529831 | 0.000000 | -0.037150 | -0.462828 | 0.035901 | -0.180738 |
| rho1.4 | -0.012380 | 0.576610 | 0 | -0.033342 | -0.468576 | 0.031662 | -0.180633 |
| rho2.0 | -0.003355 | 0.627119 | 0 | -0.028597 | -0.474003 | 0.036303 | -0.180953 |
| pool1 | -0.017470 | 0.663051 | 0 | -0.026172 | -0.601534 | 0.048216 | -0.188382 |
| postblend | -0.018045 | 0.544407 | 0 | -0.036987 | -0.461266 | 0.026652 | -0.181074 |
| v5b control | 0.037952 | 0.486441 | 0 | -0.040364 | -0.116550 | -0.034331 | -0.067481 |

All scaled MSE are exactly 0. rho>1 passes the local overall/NMAE rule but violates our frozen gate
G3 (requested amplitude ratio must be 0.90-1.05), so it cannot be packaged. pool1 passes the gates
with delta +0.0017, which is far inside the observed surrogate misordering scale (0.074).

## The two untested hypotheses now being run

Let `delta` be the fused cross-source log2FC response for a target, `Lc` the control log-bulk, and
`L = Lc + deflate(delta)`, where `deflate` removes the top-50 control-cell PCs from the response
(removing "control axes") and then a sign gate zeroes every gene without cross-source sign agreement
(`de_call_mask`, min_net=1). Amplitude alpha is then solved so that the emitted jackknife budget B
equals a probe-estimated claim, and a 0.3 blend with a learned MLP response preserves the pooled
donor norm. Finally the v9 layout re-selects which genes carry DE calls, at a fixed generic budget
of 885 calls per target, from a confidence ordering.

H4 (dense): because the v9 layout already chooses the call set independently of the bulk mask, the
sign gate and PC deflation may now only cost expression MSE and PDS, not help DE. Arm `dense` drops
BOTH the sign gate and the 50-PC deflation; arm `eta0` keeps the gate and drops only the deflation.

H5 (learned norm): `blend` keeps the pooled norm and discards the learned response's magnitude, even
though the MLP was trained in the scorer's own log-bulk space. Arm `learnednorm` substitutes the
learned norm.

Pre-registered pass rule (all must hold, same-run controls): overall strictly up; scaled NMAE >= a
same-run v5b floor; scaled MSE > 0; gates G1, G2, G3 (requested and native), G3e capped, G6, G7, G8;
B/claim and B/capped <= 1.05 with requested ratio in [0.90, 1.05].

## Questions — answer each explicitly

1. Given the v9 layout already picks calls independently, is H4 mechanistically coherent, and which
   measurement would falsify it most cheaply? If the 50-PC deflation and the sign gate should be
   separated, is `eta0` the right second arm?
2. H5 requires a learned response magnitude in log-bulk units. Is substituting the learned norm into
   a pooled-norm blend well-posed, or does it mix incommensurable scales (e.g. because pooled and
   learned are normalized per-target before mixing)? What is the cleanest one-line alternative?
3. Given the harness misordering evidence, what bar should a local gain clear before it is even worth
   a final-panel build? Is a same-run MSE>0 requirement sensible, or could a genuine improvement land
   at MSE exactly 0 again because the scorer's capped MSE clamps?
4. If both H4 and H5 fail, what is the single highest-value remaining experiment before the Nov 2
   package target? Rank at most three, with the one you would run first.
5. Anything in the above plan that is methodologically wrong or that you would refuse to sign off on?
