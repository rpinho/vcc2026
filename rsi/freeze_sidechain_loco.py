"""Freeze explicitly scoped Sidechain reproduction BEFORE biological scoring."""
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import runtime as rt
import official_scoring
import atlas_loco_baseline as atlas
import sidechain_loco_baseline as baseline


def main():
    assert not baseline.SPEC.exists(), 'Cannot overwrite a frozen recipe'
    atlas_spec = atlas.read_spec()
    upstream = rt.json.loads(Path('artifacts/sidechain_release_receipt.json').read_text())
    for name, rec in upstream['files'].items():
        assert rt.sha(baseline.RELEASE/name) == rec['sha256'], name
    tasks, sources = {}, {}
    for task in baseline.TASKS:
        entry = atlas_spec['tasks'][task]
        covered = set()
        for path in entry['sources']:
            source = baseline.MomentSource(Path(path), entry['targets'])
            covered.update(source.targets)
            sources[path] = rt.sha(path)
        assert set(entry['targets']) <= covered, (task, 'No unobserved fallback allowed')
        tasks[task] = {'targets': entry['targets'], 'destination': entry['destination'],
                       'split': 'search', 'sources': entry['sources'], 'cells_pred': 400}
    tracked = ['sidechain_loco_baseline.py', 'freeze_sidechain_loco.py',
        'verify_sidechain_loco.py', 'run_sidechain_loco.py', 'runtime.py', 'bounded_score.py',
        'bounded_de.py', 'stat_cache.py', 'benchmark.py', 'task_io.py', 'official_scoring.py',
        'baseline_protocol.py', 'audit_phase1_calibration.py']
    tracked += [str(baseline.RELEASE/p) for p in upstream['files'] if p.startswith('src/')]
    protocol = list(atlas_spec['protocol_sha256']) + ['artifacts/atlas_loco_spec.json',
        'artifacts/sidechain_release_receipt.json', 'artifacts/sidechain_kernel_audit.json',
        'artifacts/calibration_gate.json']
    record = {'utc': datetime.now(timezone.utc).isoformat(), 'baseline_id': baseline.BASELINE_ID,
        'phase': '1 scoped public-kernel reproduction, not variant selection',
        'tasks': tasks, 'source_sha256': sources,
        'implementation_sha256': {p: rt.sha(p) for p in tracked},
        'protocol_sha256': {p: rt.sha(p) for p in protocol},
        'official_environment': official_scoring.provenance(),
        'published_kernel': 'https://github.com/saberhq/sidechain/tree/'+upstream['commit'],
        'recipe': {'pseudocount': 1.0, 'minimum_source_cells': 20,
            'variance': 'Observed variance of each mean; separate delta-method terms, divided by ln(2)^2',
            'source_moment_dtype': 'Original float32, promoted to float64 for arithmetic',
            'var_floor': 'none', 'numerical_variance_clamp': 1e-6,
            'shrinkage': True, 'shrink_stage': 'pooled', 'shrink_rule': 'adaptive',
            'shrink_k': 1.0, 'gamma': 1.0, 'similarity_beta': 0.0, 'coverage_tiers': None,
            'amplitude_both_heads': 1.35, 'target_self_log2fc': -2.32,
            'bulk_anchor': 'pooled', 'emit_shape': 'controls', 'emit_lambda': 0.5,
            'control_min_libsize_exclusive': 1000.0, 'libsize_quantiles': [0.0, 1.0],
            'seed': 2026101725, 'target_seed': 'little-endian first8 SHA256(seed:task:target)',
            'dual_fit_iterations': 100, 'dual_fit_tolerance': 2e-4, 'max_projection': 0.03,
            'on_fail': 'raise', 'scatter_table': None, 'neighbours': None,
            'B_claim_amplitude_selection': None},
        'missing_input_policy': 'Fail on missing required moments, uncovered targets or failed dual fit. Unmeasured readouts abstain through the original public axis remapper. No fabricated source observations or fallback moments.',
        'adaptations': [
            'NOT SER17: missing submitted manifest and exact original source caches. No external0.1657 label.',
            'Observed matched-control float32 moment banks replace original float64 single-control PseudobulkSums. No lossless-equivalence claim.',
            'Poisson variance floor omitted because exact source control library-size/sample-count moments are unavailable; use the documented public no-floor pooling path, not a guessed floor.',
            'Use allowed GWPS and HCT116 only. H1, HEK293T, RPE1 and CD4 response branches excluded globally. No destination responses enter predictions.',
            'STRING neighbour geometry/pool/weight and1500-gene scatter table omitted; their submitted values are unknown.',
            'Amplitude1.35, pooled/adaptive and control-shaped/pooled-anchor components follow archived narrative inheritance; not a recovered manifest.',
            'Minimum20 source cells and all79 SEARCH targets/readouts inherit frozen available panels. No generic H1 fallback: every target must have observed source support.',
            'Use released project control floor1000, explicit untrimmed quantiles and local hashed seeds. Abort failed fits rather than guess submitted fallback. No knob selected by SEARCH scores.'],
        'external_official_label': None, 'exact_SER17_reproduction': False,
        'phase2_allowed': False, 'holdout_scoring_allowed': False}
    rt.save(baseline.SPEC, record)
    print('FROZEN', baseline.SPEC, rt.sha(baseline.SPEC), flush=True)


if __name__ == '__main__':
    main()
