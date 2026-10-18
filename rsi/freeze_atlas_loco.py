"""Freeze available-input AtlasShift reproduction before its biological scores."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import numpy as np
import anndata as ad
import official_scoring
import atlas_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists(), 'Existing frozen spec must not be overwritten'
    shared = Path('../public_data/model_v4_signature_banks')
    task_sources = {
        'H1_rich_search': [shared/'K562_GWPS_raw.npz', shared/'HCT116.npz'],
        'CD4_depth192_search': [Path('artifacts/cd4_banks/K562_GWPS_raw.npz'), shared/'HCT116.npz']}
    cd4 = shared/'CD4_DE_statistics.npz'
    source_files = {}
    for path in sorted({p for paths in task_sources.values() for p in paths}):
        print('Validating AtlasShift moments', path, flush=True)
        source_files[str(path)] = baseline.validate_bank(path)
    with np.load(cd4, allow_pickle=False) as data:
        for key in ['targets', 'genes', 'log2fc', 'lfcSE', 'available', 'quality_pass', 'n_cells']:
            assert key in data, key
        assert data['log2fc'].shape == data['lfcSE'].shape
        assert data['available'].shape == data['quality_pass'].shape == data['n_cells'].shape
        assert np.isfinite(data['log2fc']).all() and np.isfinite(data['lfcSE']).all()
    source_files[str(cd4)] = {'sha256': rt.sha(cd4), 'bytes': cd4.stat().st_size,
                             'exact_complete_published_source_preparation': False}
    manifest = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())
    tasks = {}
    for task in baseline.TASKS:
        contract = manifest[task]
        assert contract['split'] == 'search'
        a = ad.read_h5ad(Path('artifacts/tasks')/(task+'.h5ad'), backed='r')
        genes = a.var_names.to_numpy(str)
        a.file.close()
        sources = task_sources[task]
        coverage = {}
        for path in sources:
            with np.load(path, allow_pickle=False) as data:
                covered = set(data['targets'][data['n_cells'] >= 20])
                coverage[str(path)] = {'targets_with_at_least20_cells':
                    sum(t in covered for t in contract['targets']),
                    'measured_readout_overlap': int(np.isin(genes, data['genes'][data['measured_genes'].astype(bool)]).sum())}
        tasks[task] = {'destination': contract['context'], 'split': 'search',
                       'targets': contract['targets'], 'sources': list(map(str, sources)),
                       'weights': [2, 1], 'cd4_de': str(cd4) if contract['context'] == 'H1' else None,
                       'coverage': coverage, 'cells_pred': 400, 'readouts': len(genes)}
    tracked = ['atlas_loco_baseline.py', 'freeze_atlas_loco.py', 'verify_atlas_loco.py',
               'runtime.py', 'bounded_score.py', 'bounded_de.py', 'stat_cache.py',
               'benchmark.py', 'task_io.py', 'official_scoring.py', 'baseline_protocol.py']
    tracked += [str(baseline.RELEASE/p) for p in ['model.py', 'predict.py', 'prepare.py', 'sources.json', 'README.md']]
    assert rt.sha(baseline.RELEASE/'model.py') == rt.sha('vendor/atlas_model.py')
    protocol = ['artifacts/tasks/manifest.json', 'artifacts/benchmark_roles.json',
                'artifacts/anchor_policy.json', 'artifacts/promoter_pairs.csv']
    for task in baseline.TASKS:
        protocol += [f'artifacts/scores/{task}/anchors.json', f'artifacts/tasks/{task}.json']
    spec = {'utc': datetime.now(timezone.utc).isoformat(), 'baseline_id': baseline.BASELINE_ID,
            'phase': '1 fixed public-baseline reproduction; not variant selection',
            'seed': 20260910, 'tasks': tasks, 'source_files': source_files,
            'implementation_sha256': {p: rt.sha(p) for p in tracked},
            'protocol_sha256': {p: rt.sha(p) for p in protocol},
            'official_environment': official_scoring.provenance(),
            'published_kernel': 'https://github.com/kaipengm2/Virtual-Cell-Challenge-2026',
            'published_model_sha256': rt.sha(baseline.RELEASE/'model.py'),
            'recipe': {'prior_counts': 100000, 'source_centering': 1,
                'minimum_source_cells': 20, 'log2fc_amplitude': .6, 'bulk_delta_amplitude': .3,
                'effect_clip': 3, 'cd4_family_weight': .5, 'promoter_fraction': .15,
                'control_pool_k': 4, 'callset_control': None, 'B_claim_calibration': None},
            'missing_input_policy': 'Fail on absent required source files or moment fields. Explicit exclusions are frozen below. Unsupported target/readout pairs have zero transfer weight under the published mask; their control fallback is a model prediction, never a fabricated source observation.',
            'adaptations': [
                'Use the two frozen79-target SEARCH panels on their measured official readout intersection; not the hidden300-target panels.',
                'HEK293T and RPE1 are globally sealed; all H1 responses are excluded globally to preserve existing source isolation.',
                'CD4 is included as one DE-family vote for H1 only; both CD4 branches are excluded for CD4.',
                'Available raw GWPS/HCT116 summaries have restricted source target/readout universes and source sampling. Center over exactly these frozen available panels; no complete-source equality claim.',
                'Use the existing promoter-pair table, not an assertion of byte equality with official_pairs.csv.',
                'Control-template seed20260910 in each task; target seed uses destination:target under the published seed_for helper.',
                'Only source moments actually present in the original AtlasShift schema are required. Its release writes float32 moments; missing Sidechain float64 variances do not block this baseline.'],
            'external_official_label': None, 'exact_full_release_reproduction': False,
            'phase2_allowed': False, 'holdout_scoring_allowed': False,
            'sidechain_exact_status': 'BLOCKED on missing submitted manifest, source moments, STRING geometry/ordered842-pool/blend and scatter table/ranking/factor;126 released synthetic tests verified separately. Do not replace missing values with defaults.'}
    rt.save(baseline.SPEC, spec)
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), flush=True)


if __name__ == '__main__':
    main()
