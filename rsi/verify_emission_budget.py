"""Synthetic, label-free reproduction of the competitor B/claim emission lever.

Uses Sidechain's pinned emitter and cell-eval2's unchanged bulk-MSE path. The
amplitude solve sees only predictions, never reference counts, DE or scores.
No biological SEARCH/HOLDOUT inputs are read and no score gain is asserted.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path('vendor/sidechain_release/src').resolve()))
from sidechain.models.count_emitters import ContextProfile, PoissonEmitter
from sidechain.models.neighbour_arm import build_neighbour_arm
from cell_eval2.prep import pseudobulk_bulk_lognorm_with_moments
from cell_eval2.metrics.delta import mse_unbiased_capped, PRED_TRACE_CAP_K, _across_pert_budget
import official_scoring


def make(blocks, labels, genes):
    return ad.AnnData(sparse.csr_matrix(np.vstack(blocks).astype(np.float32)),
                      obs=pd.DataFrame({'target_gene': np.repeat(labels, [len(x) for x in blocks])},
                                       index=[f'c{i}' for i in range(sum(map(len, blocks)))]),
                      var=pd.DataFrame(index=genes))


def bulk(data):
    return pseudobulk_bulk_lognorm_with_moments(data, 'target_gene', bulk_target_sum=50_000)


def budget_loop(means, labels, genes):
    # Separate scalar centering, including exclusion BEFORE the per-gene mean.
    answer = 0.0
    for g, name in enumerate(genes):
        values = [float(means[p, g]) for p, target in enumerate(labels) if target != name]
        if len(values) < 2:
            continue
        center = sum(values) / len(values)
        answer += sum((v-center)**2 / (len(genes)-1) for v in values)
    return answer


def neighbour_check():
    rng = np.random.default_rng(211)
    axis = np.array([f'g{i}' for i in range(8)])
    pool = list(axis[:6])
    deltas = {name: rng.normal(size=len(axis)) for name in pool}
    table = {name: rng.normal(size=5) for name in pool}
    arm = build_neighbour_arm(pool, deltas.get, axis, table, k=2, w=0.35)
    rows = []
    for name in pool:
        d = deltas[name].copy(); d[list(axis).index(name)] = 0; rows.append(d)
    mean = np.array([sum(float(row[g]) for row in rows)/len(rows) for g in range(len(axis))])
    unit = lambda v: v / (np.sqrt(sum(float(x)**2 for x in v)) + 1e-12)
    misses = []
    for target in pool:
        other = [name for name in pool if name != target]
        near = sorted(other, key=lambda name: -float(unit(table[target]) @ unit(table[name])))[:2]
        n = sum((rows[pool.index(name)] - mean for name in near)) / 2
        r = rows[pool.index(target)] - mean
        expected = mean + np.linalg.norm(r)*unit(unit(r)+0.35*unit(n))
        actual = arm.fuse(target, deltas[target])
        misses.append(float(np.max(np.abs(expected-actual))))
        np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
    return {'status': 'PASS scalar numpy fusion', 'targets': len(pool), 'max_abs_error': max(misses),
            'limitation': 'Synthetic numpy table; does not test torch deserialization or missing public STRING artifact.'}


def main():
    rng = np.random.default_rng(2026101721)
    genes = np.array([f'g{i:02}' for i in range(96)])
    targets = list(genes[:12])
    fraction = rng.uniform(0.5, 2, len(genes)); fraction /= fraction.sum()
    ctrl_depth = rng.lognormal(np.log(5000), 0.35, 600)
    controls = rng.poisson(ctrl_depth[:, None] * fraction)
    directions = rng.normal(size=(len(targets), len(genes)))
    # Fix every model input before constructing independent synthetic references.
    with tempfile.TemporaryDirectory(dir='artifacts', prefix='emission_fixture_') as tmp:
        control_file = Path(tmp)/'controls.h5ad'
        make([controls], ['non-targeting'], genes).write_h5ad(control_file)
        profile = ContextProfile.from_controls(control_file, name='synthetic')

        def emit(alpha, pinned):
            blocks = []
            for i, target in enumerate(targets):
                em = PoissonEmitter(profile, seed=300+i, lam=1, bulk_anchor='pooled')
                d = alpha*directions[i]
                block = (em.emit_dual(400, d, d, on_fail='raise') if pinned else em.emit(400, d)).toarray()
                assert np.all(block >= 0) and np.array_equal(block, np.rint(block))
                if pinned:
                    expectation = block.sum() * em._fraction(d, bulk=True)
                    assert np.max(np.abs(block.sum(0)-expectation)) <= 1.000001
                    assert em.last_dual == 'dual', em.last_dual
                blocks.append(block)
            data = make(blocks, targets, genes)
            labels, means, moments = bulk(data)
            B = budget_loop(means, labels, genes)
            exclusions = {i: list(genes).index(label) for i, label in enumerate(labels)}
            official_B = _across_pert_budget(means, exclusions, np.full(len(labels), 1/(len(genes)-1)))
            np.testing.assert_allclose(B, official_B, rtol=2e-11, atol=2e-13)
            self_claim = float(moments.jk.sum()/(len(genes)-1))
            return data, labels, means, moments, B, self_claim

        # Reference-free bisection over this emitter's own B / sum(w*jk_pred).
        low, high = 0.0, 0.05
        trials = []
        for _ in range(16):
            alpha = (low+high)/2
            result = emit(alpha, True)
            ratio = result[4]/result[5]
            trials.append({'alpha': alpha, 'prediction_B_selfclaim': ratio})
            if ratio < 1.0: low = alpha
            else: high = alpha
        alpha = high
        pinned = emit(alpha, True)
        assert 1.0 <= pinned[4]/pinned[5] < 1.02, pinned[4]/pinned[5]
        ordinary = emit(alpha, False)
        zero = emit(0.0, True)

    # Realization/metrics come AFTER the amplitude is irrevocably selected.
    real_blocks = [controls]
    for i in range(len(targets)):
        prob = fraction*np.exp2(0.02*directions[i]); prob /= prob.sum()
        depths = rng.lognormal(np.log(5000), 0.35, 400)
        real_blocks.append(rng.poisson(depths[:, None]*prob))
    real = make(real_blocks, ['non-targeting']+targets, genes)
    rlabels, rmeans, rmom = bulk(real)
    ri = {name:i for i,name in enumerate(rlabels)}
    rows = {}
    for name, result in [('pinned_zero_amplitude', zero), ('ordinary_selected_amplitude', ordinary), ('pinned_selected_amplitude', pinned)]:
        data, labels, means, moments, B, selfclaim = result
        real_jk = np.array([rmom.jk[ri[label]] for label in labels])
        claim_terms = np.minimum(moments.jk, PRED_TRACE_CAP_K*real_jk)
        claim = float(claim_terms.sum()/(len(genes)-1))
        ratio = B/claim
        retained = min(1.0, ratio)
        actual = mse_unbiased_capped(comparator='bulk_lognorm', pred_bulk=(labels,means),
            real_bulk=(rlabels,rmeans), pred_moments=moments, real_moments=rmom, genes=genes,
            pred_bulk_full=(labels,means))
        error = 0.0
        for i, label in enumerate(labels):
            squared = sum(float(means[i,g]-rmeans[ri[label],g])**2
                          for g in range(len(genes)) if genes[g]!=label)
            expected = (squared-retained*claim_terms[i]-real_jk[i])/(len(genes)-1)
            error = max(error,abs(expected-actual[label]))
            np.testing.assert_allclose(actual[label], expected, rtol=2e-11, atol=2e-13)
        rows[name] = {'B': B, 'prediction_only_claim': selfclaim, 'scorer_claim':claim,
            'B_scorer_claim': ratio, 'retained_correction':retained,
            'raw_capped_MSE':float(np.mean(list(actual.values()))), 'independent_arithmetic_max_abs_error':error}
    assert rows['pinned_selected_amplitude']['retained_correction'] == 1
    assert rows['pinned_zero_amplitude']['B_scorer_claim'] < 0.01
    receipt = {'status':'PASS synthetic mechanism only', 'utc':datetime.now(timezone.utc).isoformat(),
        'neighbour':neighbour_check(), 'selected_alpha':alpha, 'prediction_only_trials':trials,
        'arms':rows, 'official':official_scoring.provenance(),
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'release_receipt_sha256':hashlib.sha256(Path('artifacts/sidechain_release_receipt.json').read_bytes()).hexdigest(),
        'reference_used_in_amplitude_selection':False, 'biological_search_or_holdout_read':False,
        'limitation':'Synthetic fixture, not a biological prediction or scaled official gain. Self-claim is a conservative bound on min(jk_pred,k*jk_real); target B/selfclaim=1 ensures full correction here. No full .vcc, calibration repair or variant selection.'}
    Path('artifacts/emission_budget_fixture.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)


if __name__ == '__main__':
    main()
