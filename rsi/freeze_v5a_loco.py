"""Freeze the round29 v5a port contract before any new biological score."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import anndata as ad
import numpy as np
import official_scoring
import provenance
import v5a_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists(), 'Cannot overwrite a frozen recipe'
    prior = rt.json.loads(Path('artifacts/phase1_round29_audit.json').read_text())
    assert prior['status'] == 'PASS'
    for group in ['preserved_protocol_sha256', 'historical_identities', 'v5a_input_snapshot_sha256']:
        for path, digest in prior[group].items():
            assert rt.sha(path) == digest, path
    assert rt.sha('vendor/v4_model.py') == rt.sha(baseline.SNAPSHOT/'model.py')
    assert rt.sha('vendor/atlas_model.py') == rt.sha(baseline.SNAPSHOT/'vendor/atlas_model.py')
    shared = Path('../public_data/model_v4_signature_banks')
    source_names = ['K562_GWPS_raw', 'HCT116', 'CD4_Rest', 'CD4_Stim48hr',
                    'CD4_Stim8hr', 'CD4_DE_statistics']
    sources = {str(shared/(n+'.npz')): rt.sha(shared/(n+'.npz')) for n in source_names}
    sources['artifacts/cd4_banks/K562_GWPS_raw.npz'] = rt.sha('artifacts/cd4_banks/K562_GWPS_raw.npz')
    sources['artifacts/cd4_banks/HCT116.npz'] = rt.sha('artifacts/cd4_banks/HCT116.npz')
    manifest = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())
    tasks = {}
    for ci, task in enumerate(baseline.TASKS):
        contract = manifest[task]
        assert contract['split'] == 'search' and len(contract['targets']) == 79
        path = Path('artifacts/tasks')/(task+'.h5ad')
        receipt = rt.json.loads(path.with_suffix('.json').read_text())
        assert rt.sha(path) == receipt['sha256']
        a = ad.read_h5ad(path, backed='r')
        genes = a.var_names.to_numpy(str).tolist()
        a.file.close()
        # The parent has documented per-readout abstention; every target needs
        # actual >=20-cell support in at least one eligible core bank.
        banks = [shared/'K562_GWPS_raw.npz', shared/'HCT116.npz'] if ci == 0 else [
            Path('artifacts/cd4_banks/K562_GWPS_raw.npz'), Path('artifacts/cd4_banks/HCT116.npz')]
        covered = set()
        for bank in banks:
            with np.load(bank, allow_pickle=False) as data:
                covered.update(data['targets'][data['n_cells'] >= 20])
        assert set(contract['targets']) <= covered
        probe = np.sort(np.random.default_rng(ci).choice(79, 40, replace=False))
        tasks[task] = {'destination': contract['context'], 'split': 'search',
            'targets': contract['targets'], 'genes': genes, 'cells_pred': 400,
            'probe_indices': probe.tolist(), 'probe_rng_seed': ci,
            'historical_implementation_identity': provenance.identity(task)}
    tracked = ['v5a_loco_baseline.py', 'freeze_v5a_loco.py', 'verify_v5a_loco.py',
        'run_v5a_loco.py', 'runtime.py', 'recipes.py', 'provenance.py', 'vendor/v4_model.py',
        'vendor/atlas_model.py', 'vendor/adaptive_shrink.py', 'bounded_score.py',
        'bounded_de.py', 'stat_cache.py', 'benchmark.py', 'task_io.py',
        'official_scoring.py', 'baseline_protocol.py', 'audit_phase1_calibration.py']
    tracked += list(prior['v5a_input_snapshot_sha256'])
    protocol = list(prior['preserved_protocol_sha256']) + ['artifacts/phase1_round29_audit.json',
        'artifacts/promoter_pairs.csv', 'artifacts/task_io_fixture.json',
        'artifacts/bounded_score_fixture.json']
    for task in baseline.TASKS:
        protocol += [f'artifacts/tasks/{task}.json', f'artifacts/scores/{task}/anchors.json']
    protocol += list(prior['historical_identities'])
    rt.save(baseline.SPEC, {'utc': datetime.now(timezone.utc).isoformat(),
        'baseline_id': baseline.BASELINE_ID, 'phase': '1 fixed reproduction extension; no variant selection',
        'tasks': tasks, 'implementation_sha256': {p: rt.sha(p) for p in tracked},
        'protocol_sha256': {p: rt.sha(p) for p in protocol}, 'source_sha256': sources,
        'official_environment': official_scoring.provenance(),
        'recipe': {'arm': 'a', 'parent': 'v4a local port', 'rho': 1.0, 'probe_size': 40,
            'control_pool_k': 4, 'template_seed': baseline.SEED,
            'target_seed': 'atlas.seed_for(destination:target,20260910)',
            'claim': 'mean of40 own jackknifes times sum of all79 row weights',
            'centre': 'log1p(50000*input_control_pooled_probability)',
            'solve': 'Original solve_alpha lo0 hi8 iters60; abort if unattainable',
            'emitter': 'Unmodified mean_counts and repair_bulk_totals; fitted expectation totals',
            'cap1': False, 'deflation': None, 'callset_control': None},
        'adaptations': [
            'Two frozen79-target measured-axis SEARCH panels, not hidden300-target submission.',
            'Exact original v5a helper/model hashes; same local v4a family exclusions and source masks.',
            'H1/HEK293T/RPE1 responses excluded globally; CD4 response branches excluded for CD4.',
            'Input controls only. Four-cell pooled template as actual v5a builder, despite its stale module-level pool_k1 prose.',
            'Template seed20260910 in both tasks as local v4a; context probe RNG0/1; destination-based target seeds.',
            'Probe estimates pre-amplitude own claim; final actual emitted claim is recorded, not forced to rho by a new iterative solve.',
            'Original fit stops at its tolerance; exact integer repair binds fitted expected columns, not exactly400 times requested probability.',
            'Abort missing support, identity mismatch or unattainable budget; no score-selected fallback.'],
        'external_official_label': None, 'predicted_official_score': None,
        'phase2_allowed': False, 'holdout_scoring_allowed': False})
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), flush=True)


if __name__ == '__main__':
    main()
