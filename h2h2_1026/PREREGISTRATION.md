# Frozen out-of-sample H1 head-to-head

Registered 2026-10-26 after input feasibility, before any model selection or scoring.
One experiment, CPU only, one sequential worker/BLAS thread, nice >=10, >=30,000,000,000 free bytes. No network, submission, upload, package, promotion, leaderboard, or other job. Only this directory is writable. The explicit user task supersedes older mission instructions to launch jobs or build packages and older live-prior-art requirements. PRIOR_ART.md already exists and is retained.

## Cohorts and representation

Fit H1_full's 168 unique H1 targets only, all 18,077 measured genes, in inherited response coordinates. Predict all 79 official H1 SEARCH targets, with empty H1 training/test intersection asserted at every fitting scope. Other-family same-target donor responses are legal inference features; H1 recipient banks are excluded wholesale. No held-out H1 response, DE, cell count, or normalizer is a model feature. Training-only centering, response basis, feature mean/SD, and ridge intercept. Shared legal input-control template is the frozen renderer's normalisation path, never treated-cell truth.

Bank 18,080 -> scoring 18,077 genes by symbol; manifold bank 18,533 -> its 18,077 measured genes -> scoring order. Map in results/gene_alignment.json, hash cc912b655d9f987e7d028e9deaaeeba0e756492273b9202d24e368d5c9c7b1b4. All scoring genes must exist and be measured. No gene imputation. Truth diagnostic: log1p(50000*p)-log1p(50000*ctrl), each CPM row explicitly normalized to sum one (original sums ~1e6). Official six metrics use actual held-out task cells, not this diagnostic pseudobulk.

STRING missing KDR/SOX2 (77/79 covered); inherited embed() returns an explicitly documented zero vector for those two. Keep all 79 in all arms and primary scores; report coverage and their individual results. Training STRING misses KRT18/SALL4. Zero-feature fallback is declared, not an observed embedding or silent imputation.

## Frozen six arms

B0: exact pinned shipped response and unchanged shipped v9 count emitter.
B1: F1 STRING-only model; no incumbent features.
B2: F3 STRING + projected unit pooled-donor direction.
B3: exact inherited pooled() response (family-balanced donor banks), continuity comparator.
B4: full 168 training response mean, same vector for all 79.
B5: 0.5 * B0 + 0.5 * B4 in response space. Weight 0.5 fixed a priori before scoring, not selected from any target outcome. This tests a concrete shrinkage alternative, not an oracle or optimally tuned shrinkage envelope.

Reuse AST-extracted basis/ridge_path/features/split/pooled/embed from pinned manifold run.py, not its old fold predictions or incident artifacts. Select k in [8,16,32,64] and alpha in [.01,.1,1,10,100,1000,10000] independently for F1/F3 with 3-fold training-only split(targets,n_cells,3,20261126). For fold ii, basis seed 20261026+ii; all response transforms refit inside training split. Minimize summed validation expression SSE across the three folds, inherited objective and tie break (smaller k, larger alpha). No feature-arm selection: report both F1/F3. Fit final basis once on all 168 with seed 20261035; refit selected ridge once per arm. Save complete basis/mean, fitted normalization/weights, features, train/test target identities, selections/grid/scope certificates, predictions. No retuning after outcomes.

## Rendering and official scoring

Identical pinned gatecall cell_eval2 vcc2026 CPU/pdex config, num_threads=gather_threads=1; frozen model.emit anchor_single_mean + layout.emit, original 400-cell control template/depths, seed 20260910, call budget 885, own-gene protection. For every arm use the same confidence formula from legal excluded-H1 donor z scores and that arm's predicted sign; recompute sign-dependent confidence rather than freezing B0's gene signs onto competitors. B0 confidence hash, count digest, and all six raw components must match pinned historical B0 exactly (tolerance 1e-12). No adaptive DE gating or magnitude/blend selection. Frozen reference objects are input controls only; assert their identity. Same held-out reference cells and matching local H1 SEARCH baseline/replicate anchors for all arms. Bounded official reduction/DE adapters pinned; previously bitwise-validated fixtures and official package file identity must match. Save each arm's gene response NPZ, all target-level official metric vectors as NPZ and tidy CSV, raw/scaled six metrics and score_avg. Save compact rendered count sufficient statistics/DE evidence and count digests; full count matrices are regenerable from saved responses and frozen inputs and need not be retained.

## Inference and mechanical verdict

20,000 paired target-level bootstrap draws, NumPy default_rng seed 20261026, percentile 95% linear quantiles. Identical sampled indices for all arms. Aggregate official finite metric cohorts and paired-finite ratio-of-sums for expression MSE; retain fixed full-panel PDS contributions and anchors (conditional target bootstrap, not recomputed retrieval competition). Reapply official nonlinear score_one scaling to each draw, composite = six scaled components' mean. No dropping zero-call or missing-STRING targets. Save every target vector and all bootstrap primary draws. H1 is the only family tested; no multi-family extrapolation.

Primary for each B1/B2 minus B0: positive raw direction fidelity, raw PDS cosine, and scaled score_avg deltas with every paired CI lower bound >0. Secondary: B2 minus B4 and B5 on those same three tests. For a B1 adoption also demand B1 minus B4/B5 on all three, so no manifold arm escapes the mandatory artifact controls. Report B3 comparisons and B4/B5 versus B0 descriptively; they cannot rescue a failure. No multiplicity adjustment (two prespecified candidates, three conjunctive endpoints); report this limitation.

ADOPT iff at least one manifold arm passes all three primary and both controls' three tests. Otherwise REJECT and permanently close manifold/STRING in the report headline. If B2 beats B4 but not B5 (or converse), state it plainly. Failure of any CI is binding, no interpretation-based override. Missing required inputs, identity/resource/provenance failure => INPUT BLOCKED/NOT EVALUATED, with exact cause; never infer impossibility merely from disjoint cohorts.

Separate verify.py must recompute aggregation, all primary bootstrap deltas/CIs, exact binary64 gate comparisons, alignment/cohort and saved model/control identities from saved arrays and pin hashes without importing search/analysis code. Numerical reproduction is not a formal real-arithmetic or biological proof. No leaderboard/novelty claim.
