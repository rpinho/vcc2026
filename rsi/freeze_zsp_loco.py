"""Freeze one public-method baseline before biological scoring or head audit."""
from datetime import datetime, timezone
from pathlib import Path
import runtime as rt
import numpy as np
import official_scoring
import zsp_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists(), 'Never overwrite a frozen spec'
    manifest = rt.json.loads(Path('artifacts/tasks/manifest.json').read_text())
    receipt = rt.json.loads(Path('artifacts/zsp_release_receipt.json').read_text())
    assert receipt['status'] == 'PASS pinned code identities only'
    tasks = {}
    inputs = list(baseline.SOURCES.values())
    protocol = ['artifacts/tasks/manifest.json', 'artifacts/benchmark_roles.json',
                'artifacts/anchor_policy.json', 'artifacts/calibration_gate.json']
    for task in baseline.TASKS:
        m = manifest[task]
        assert m['split'] == 'search' and m['context'] in ['H1', 'CD4']
        tasks[task] = {'targets': m['targets'], 'destination': m['context'], 'split': 'search'}
        inputs.append(f'artifacts/tasks/{task}.h5ad')
        protocol += [f'artifacts/tasks/{task}.json', f'artifacts/scores/{task}/anchors.json']
        protocol += [f'artifacts/scores/{task}/{r}_s42.json' for r in ['v1','v2b','v3b','v4a','v4b']]
    implementation = ['zsp_loco_baseline.py', 'freeze_zsp_loco.py', 'verify_zsp_loco.py',
        'run_zsp_loco.py', 'runtime.py', 'bounded_score.py', 'bounded_de.py', 'stat_cache.py',
        'benchmark.py', 'task_io.py', 'official_scoring.py', 'baseline_protocol.py',
        'artifacts/zsp_release_receipt.json']
    for p, item in receipt['files'].items():
        path = baseline.RELEASE/p
        assert rt.sha(path) == item['sha256']
        implementation.append(str(path))
    # Metadata-only coverage is descriptive; it never changes the chosen banks.
    coverage = {}
    for name, path in baseline.SOURCES.items():
        with np.load(path, allow_pickle=False) as bank:
            targets = set(bank['targets'][bank['n_cells'] >= 10])
            coverage[name] = {t: sum(g in targets for g in tasks[t]['targets']) for t in tasks}
    spec = {'utc': datetime.now(timezone.utc).isoformat(), 'baseline_id': baseline.BASELINE_ID,
        'phase': '1 fixed public-method reproduction; not configuration selection',
        'seed': 0, 'cells_pred': 400, 'tasks': tasks, 'sources': baseline.SOURCES,
        'source_target_coverage': coverage,
        'recipe': {'source_profile': 'log2(CP10k pooled raw mean +0.1) difference',
            'minimum_source_cells': 10, 'weights': 'exact equal; temperature infinity',
            'source_centering': False, 'norm': 'mean L2 norm among sources measuring target, full union axis',
            'generic_fallback': 'published per-source gene median, only wholly unseen target',
            'emission': 'released resample-rescale control cells, stochastic integer rounding',
            'B_claim_calibration': None, 'call_set_control': None},
        'adaptations': [
            'Four available source lines K562 GWPS/HepG2/Jurkat/HCT116, never an exact eight-source reproduction.',
            'RPE1 and HEK293T are globally sealed; H1 response family excluded globally. CD4 family excluded for CD4 destination.',
            'CD4 omitted also for H1: source bank lacks released Stim48hr per-gene median baseMean needed for its reconstructed pseudobulk. Do not replace with raw Rest cells or a different chemistry/condition.',
            'Existing source subsets and float32 summaries, pooled native controls and original source filters retained; not full released source preparation.',
            'Raw legacy means are used instead of quantized signatures. HCT116 target count sums and global control probabilities enter the same CP10k normalization; matched-control moments are not substituted.',
            'Unique symbol-native measured source axes replace released Ensembl-to-symbol mapping. Missing source genes remain NaN on full union before prediction; missing destination readouts receive zero predicted LFC.',
            'Exact equal weights use temperature infinity, rather than asserting byte identity to the submitted finite-temperature config.',
            'Seed0 RNG resets per separate SEARCH task and follows frozen manifest target order; original multi-context script uses one sequential RNG.',
            'No generic-response subtraction: actual NormRestoredTransfer uses raw lfc, even though SourceEffects also stores specific effects.',
            'No pinning/dual-moment emitter, source/weight/amplitude tuning, control-PC removal, call gates or benchmark changes.'],
        'release_commit': receipt['commit'],
        'published_kernel': 'https://github.com/MMatinGerami/zero-shot-perturbation',
        'implementation_sha256': {p: rt.sha(p) for p in implementation},
        'protocol_sha256': {p: rt.sha(p) for p in protocol},
        'input_sha256': {p: rt.sha(p) for p in inputs},
        'official_environment': official_scoring.provenance(),
        'external_official_label': None, 'exact_full_release_reproduction': False,
        'phase2_allowed': False, 'holdout_scoring_allowed': False}
    rt.save(baseline.SPEC, spec)
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), coverage, flush=True)


if __name__ == '__main__':
    main()
