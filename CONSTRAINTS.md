# Constraints — 2026-10-18

Official validation is the fitness function. Expected gains below are UNVERIFIED. No local proxy can rank the next arms.

| Candidate constraint | Evidence | Verdict |
|---|---|---|
| Method | v5b .1528; its MSE .082 vs .442 median; v5 emitter sets per-cell CPM and bulk to effectively the same profile; pool4 limits JK. `model_v5/NOTES.md`, `model.py:mean_counts`, live leaderboard. Independent dual moments and selective scatter are not yet official-tested here. | **binding** |
| Data | Five public source families already available (~0.8 GB summary bank). Leader .470871 has proprietary data; same-PDS Bill M .2516 identifies a gap beyond discrimination but does not isolate data from method. Our same-data v4a→v5b gain .00895 demonstrates method room. Missing KOLF/HipSci/VIPerturb/LINCS. Most current banks retain437 requested rows, H1 retains614; this is not a full all-target GWPS/Orion atlas. Global source centering is limited to retained eligible rows. Bank axes/counts are in model_v7/artifacts/bank_scope.json. | **partly**; binding for reproducing Illumina, not explanation for all our gap |
| Training | Existing source transfer .1546/.1684 public baselines beat earlier neural/scratch .088/.062. No evidence a larger trained model is required for the next improvement. PIE predictions already exist at `model_v6/pie/out/ens4_{A,B,C}.npz`. | **not binding** for next two; partly for richer biological prediction |
| Tests / validation | Shape and integer/package checks work. LOCO Spearman .20 against preregistered .80 bar, PDS -.32, MSE .78; `proxy/CALIBRATION.md`. HCT116 anchor invalid; CD4 negative. | **binding** for offline candidate ranking, not for format verification |
| Fitness function | Six equally weighted scaled components; official scores available only after upload. Labels hidden; final new contexts blind, last entry counts. No scalar proxy substitutes for official validation. | **binding** |
| Compute | Resource-limited local CPU; full builds use one numerical worker and BLAS=1. Identical-size probes took 3–126 seconds under contention. Summary-transfer/emission needs no GPU. PRiMeFlow reference pretraining uses 12 GPUs and ~244 GB compressed corpus. | **binding for observed wall time**; not an architectural blocker for these two; binding for reproducing large generative pretraining today |
| Memory | Full dense validation is 6.67 billion entries, ~26.7 GB float32 before copies; streamed 400-row blocks ~60 MB float64, bounded CLI adapter already v5-verified. | **partly**, streaming handles output size; global compression/swap contention still slows execution |
| Submissions budget | 2 per UTC day. Next pair assigned 2026-10-19 00:03 UTC and following available slot; no submissions in this job. Cannot identify independent effects of every combined lever with two entries. | **binding** |
| Process / RSI gate | User reports 26 rounds no candidates; at the initial audit `rsi/NEXT_SUBMISSIONS.md` was already round55 and still CLOSED. Historical calibration is distinct from failed proxy and must not be rewritten. Loop can keep auditing while manual official-feedback candidates use the available slots. | **binding** on RSI production; not a reason to block these authorized builds |
| Knowledge | 1,569 live entries, 610 descriptions; no verified executable top-.30 recipe. Named shared families do not establish common authorship; score clusters cannot identify causal mechanisms. Public scorer has no listed issues in API, but changelog contains issue-number decisions. | **partly**; substantial implementation leads, no secret winning recipe recovered |

## Component decomposition

Snapshot `competitors/raw/lb-2026-10-18.json`, sorted by score. v5b is now rank499 with the same .1528 score. Ranks11–50 median is a component-wise synthetic comparator, not a real team's score.

| Scaled component | Our v5b | Median 11–50 | Overall points at stake (gap / 6) | Plausible lever |
|---|---:|---:|---:|---|
| pds | 0.680335 | 0.810375 | +0.021673 | Full bulk direction; deflation; wider measured-response coverage |
| mse | 0.082141 | 0.442478 | +0.060056 | Separate bulk/DE moments; selected-gene dispersion; capped B/claim |
| nmae | 0.108512 | 0.236424 | +0.021319 | Independent per-cell CPM amplitude; transferable magnitude calibration |
| fid | -0.030517 | 0.043181 | +0.012283 | Source-consensus DE power while avoiding recall collapse |
| reach | 0.076527 | 0.263016 | +0.031081 | Confidence ranking through Wilcoxon p-values, not only effect amplitude |
| jac | -0.000258 | 0.011542 | +0.001967 | DE-significance selection; PIE; small available payoff |

Total synthetic gap = 0.148379 overall. MSE is the largest single opportunity; it is not all of the gap. The top-tier simultaneous MSE/fidelity/reach result is observed, but the exact mechanism behind every private top entry is unknown.

## Process proposal — NOT APPLIED

In `rsi/supervisor/round_prompt.md`, retain the failed calibration record and HOLDOUT seal, but allow a separate manual official-feedback lane: propose one fixed label-free candidate and one interpretable contrast per UTC day, check counts/moments/format, record the official six-vector after Ricardo submits. A gate may prevent claims of proxy ranking while still allowing bounded experiments against the actual fitness function. No supervisor prompt, journal, lock, scorer or historical gate was edited by this job. Only the requested dated candidate section in NEXT_SUBMISSIONS.md will be appended after packaging.

Source-bank scope audit (public summaries, >=20-cell eligibility, before any further source-specific quality gate): GWPS269/300 official targets on8246 genes; HCT116290/300 and HEK293T294/300 on18536 genes; H1 has137 eligible bank rows but only11/300 official targets; CD4 DE291/300 on10282 genes. Most banks have437 cached rows; H1 has614. This is not complete source-screen centering/coverage; it constrains exact public-pipeline reproduction and may matter independently of emission. No source-bank or RSI metadata repair was applied in this job.
