"""Local DE check for the emitter decision: SAME billmA biology (LOCO fold weights, v5b profile, 0.3 learned blend),
three emissions, scored by the official cell-eval2 code paths on a frozen RSI SEARCH task with REAL perturbed cells.

  incumbent : v5b emission (4-cell pooled template, one-moment pinned totals) at billmA's alpha   [billmA local path]
  sel_fixed : selective scatter (single-cell template off-support, 4-cell on cross-source-supported genes),
              AtlasShift dual-moment fit with both moments = the same profile, SAME alpha        [scatter only]
  sel_recal : sel_fixed with alpha re-solved so B = 0.97 x min(own probe claim, control-400 cap)   [the (b) recipe]

Reports the official six metrics (raw + scaled with the task's measured anchors), the raw normalised MSE, the emitted
B, claim and B/claim, and the control-400 cap. Writes ONLY under final_push/artifacts/de_check/<task>/<arm>/.
Never reads the sealed HOLDOUT tasks. Usage: python de_check.py --task H1_rich_search --arms incumbent,sel_fixed,sel_recal
"""
import os, sys, types, gc, time, json, argparse
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; VCC = HERE.parent
OUT = HERE / 'artifacts' / 'de_check'; OUT.mkdir(parents=True, exist_ok=True)
for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMBA_NUM_THREADS', 'POLARS_MAX_THREADS']:
    os.environ[k] = '2'
os.environ['TMPDIR'] = str(OUT / 'tmp'); os.environ['NUMBA_CACHE_DIR'] = str(OUT / 'numba_cache')
(OUT / 'tmp').mkdir(exist_ok=True); (OUT / 'numba_cache').mkdir(exist_ok=True)
import torch; torch.set_num_threads(2); torch.set_num_interop_threads(1)
sys.path.insert(0, str(VCC / 'model_billm'))
import common as job                                  # model_billm's job-local helpers (ROOT = model_billm, read-only use)
from infer import Ensemble
from network import blend


def rsi_context():
    """The frozen RSI import context, exactly as model_billm/evaluate.py builds it (no RSI status writes)."""
    rsi = VCC / 'rsi'; os.chdir(rsi); sys.path.insert(0, str(rsi))
    runtime = types.ModuleType('runtime'); runtime.ROOT = rsi; runtime.json = job.json; runtime.sha = job.sha; runtime.save = job.save
    def guard(extra=0):
        import shutil
        if shutil.disk_usage(VCC).free - extra < 25 * 1024 ** 3:
            raise RuntimeError('25 GiB disk reserve')
    runtime.guard = guard
    def forbidden(*a, **k):
        raise RuntimeError('RSI status write forbidden')
    runtime.status = forbidden; sys.modules['runtime'] = runtime
    sys.path.insert(0, str(rsi / 'vendor')); sys.path.insert(0, str(rsi / 'vendor/deps'))
    import benchmark, official_scoring as scoring, bounded_score, bounded_de, v5b_loco_baseline as baseline
    bounded_de.install(); bounded_de.rt = types.SimpleNamespace(**{**vars(runtime), 'ROOT': OUT})
    return benchmark, scoring, bounded_score, baseline, rsi


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--task', default='H1_rich_search')
    ap.add_argument('--arms', default='incumbent,sel_fixed,sel_recal'); ap.add_argument('--rho', type=float, default=0.97)
    a = ap.parse_args(); task = a.task
    benchmark, scoring, bounded_score, baseline, rsi = rsi_context()
    helpers = baseline.helpers
    from scipy import sparse
    contract = baseline.read_spec(); model, ctrl, real = benchmark.setup(task); assert model['split'] == 'search'
    fold = model['context']; ensemble = Ensemble(fold); assert fold not in {b['family'] for b in ensemble.src.banks}
    genes = ctrl.var_names.to_numpy(str); canonical = ensemble.src.genes; lookup = {g: i for i, g in enumerate(canonical)}
    assert all(g in lookup for g in genes); indices = np.array([lookup[g] for g in genes], dtype=np.int64)
    b = np.asarray(ctrl.X.sum(0, dtype=np.float64)).ravel(); b /= b.sum()
    canonical_control = np.zeros(len(canonical)); canonical_control[indices] = b
    targets = np.array(model['targets']); T = len(targets)
    learned = ensemble.predict(targets, canonical_control)['learned'][:, indices]
    del ensemble; gc.collect()
    # incumbent profile and amplitude (billmA's exact local path), plus the cross-source support mask
    q4, depth, bf, amp = baseline.calibrated_profile(ctrl, task, contract)
    _, _, L, Lc, support, pinfo = baseline.profile(ctrl, task, contract)
    new = helpers.from_logbulk(Lc + blend(helpers.logbulk(bf) - Lc, learned, .3))      # billmA incumbent profile
    with baseline.local_v4(model['context']) as v4:
        q1, _, _, _ = v4.template(ctrl, genes, pool_k=1, seed=baseline.SEED)
    w = helpers._row_weights(T, len(genes), helpers.excl_cols(targets, genes))
    # control-400 cap at this task's depth (reference controls, never treated cells)
    labels = real.obs.target_gene.astype(str).to_numpy(); R = sparse.csr_matrix(real.X[labels == 'non-targeting'])
    rng = np.random.default_rng(73591); cap_samples = []
    for _ in range(8):
        X = R[rng.choice(R.shape[0], min(400, R.shape[0]), replace=False)]
        cap_samples.append(float(helpers.jackknife_correction(X, np.zeros(X.shape[0], np.intp), 1, helpers.TS, chunk=32)[0]))
    cap_jk = float(np.mean(cap_samples)); control400_claim = cap_jk * w.sum()
    import importlib; atlas = baseline.recipes.v4.atlas

    def jk(x):
        return float(helpers.jackknife_correction(sparse.csr_matrix(x), np.zeros(x.shape[0], np.intp), 1, helpers.TS, chunk=32)[0])

    def template_for(i):
        q = q1.copy(); m = support[i]; q[:, m] = q4[:, m]; q /= q.sum(1, keepdims=True); return q

    def emit_sel(i, bvec):
        seed = atlas.seed_for(model['context'] + ':' + targets[i], baseline.SEED)
        x = atlas.dual_moment_counts(template_for(i), bvec, bvec, depths=depth, seed=seed)
        x = np.asarray(x); assert x.dtype.kind in 'iu' and np.array_equal(x.sum(1), depth); return sparse.csr_matrix(x)

    def emit_inc(i, bvec):
        return baseline.emit(q4, depth, bvec, model['context'], targets[i])

    anchors = job.json.loads((rsi / f'artifacts/scores/{task}/anchors.json').read_text())
    assert anchors['real_sha256'] == scoring.matrix_id(real)
    results = {}
    for arm in a.arms.split(','):
        t0 = time.time(); out = OUT / task / arm; out.mkdir(parents=True, exist_ok=True)
        prof = new; alpha_factor = 1.0; calib = None
        if arm == 'sel_recal':
            # re-solve the amplitude about the control for the selective template: B = rho x min(probe claim, cap)
            probe = np.array(contract['tasks'][task]['probe_indices'], dtype=int)[:8]; Lnew = helpers.logbulk(new); hist = []; f = 1.0
            for it in range(4):
                cand = helpers.from_logbulk(helpers.rescale(Lnew, f, Lc))
                jks = np.array([jk(emit_sel(int(i), cand[int(i)]).toarray()) for i in probe])
                claim = float(np.mean(np.minimum(jks, cap_jk)) * w.sum())
                f_next, B = helpers.solve_alpha(Lnew, targets, genes, a.rho * claim, center=Lc)
                hist.append(dict(iteration=it, factor=f, probe_claim=float(jks.mean() * w.sum()), claim_capped=claim, factor_next=f_next, B_next=B))
                print(arm, 'calibration', hist[-1], flush=True)
                done = abs(f_next / f - 1) < .005; f = f_next
                if done:
                    break
            prof = helpers.from_logbulk(helpers.rescale(Lnew, f, Lc)); alpha_factor = f; calib = hist
        emitter = emit_inc if arm == 'incumbent' else emit_sel
        P = np.zeros((T, len(genes))); jks = np.zeros(T)

        def iterator():
            for i, t in enumerate(targets):
                x = emitter(i, prof[i]); jks[i] = jk(x); P[i] = np.asarray(x.sum(0)).ravel()
                if i % 10 == 0:
                    print(arm, 'EMIT', i, '/', T, 'elapsed', round(time.time() - t0), flush=True)
                yield t, x
        cfg = scoring.config(cache_real=str(OUT / task / 'real_cache'))
        tidy, digest = bounded_score.score_iter(real, iterator(), cfg)
        raw = benchmark.aggregate(tidy); scaled, overall = scoring.scale(raw, anchors['baseline'], anchors['replicate'])
        Pn = P / P.sum(1, keepdims=True); B = float(helpers.budget(Pn, targets, genes))
        claim = float((w * jks).sum()); claim_capped = float((w * np.minimum(jks, cap_jk)).sum())
        rec = dict(task=task, arm=arm, fold=fold, raw=raw, scaled=scaled, overall=overall, prediction_sha256=digest,
                   B_emitted=B, claim_emitted=claim, claim_capped_emitted=claim_capped, ratio_B_over_claim=B / claim,
                   ratio_B_over_capped=B / claim_capped, control400_claim=control400_claim, control400_jk_samples=cap_samples,
                   alpha_billmA=amp['alpha'], alpha_factor=alpha_factor, calibration=calib, support_genes_mean=float(support.sum(1).mean()),
                   seconds=round(time.time() - t0), status='official scorer on a public LOCO SEARCH task; not an official challenge score')
        tidy.write_csv(out / 'tidy.csv'); job.save(out / 'score.json', rec); results[arm] = rec
        print('SCORED', task, arm, 'overall', round(overall, 5), {k: round(v, 4) for k, v in scaled.items()},
              'rawMSE', round(raw['expr_mse_unbiased_capped_norm'], 4), 'B/claim', round(B / claim, 3), 'B/cap', round(B / control400_claim, 3), flush=True)
        gc.collect()
    job.save(OUT / task / 'summary.json', results)
    print('DONE', task, flush=True)


if __name__ == '__main__':
    main()
