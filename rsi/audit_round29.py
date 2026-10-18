"""Close the fixed public-kernel baseline and freeze the next Phase1 input evidence.

No biological metric or prediction is recomputed; no recipe, task, scale or gate
is changed. The v5a package is hashed, not decoded or treated as local truth.
"""
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import gc
import math
import shutil

import runtime as rt
import benchmark
import official_scoring as scoring
import provenance
from audit_phase1_calibration import exact_scaled


TASKS = ['H1_rich_search', 'CD4_depth192_search']
RECIPES = ['v1', 'v2b', 'v3b', 'v4a', 'v4b']
OUT = Path('artifacts/phase1_round29_audit.json')
FROZEN = Path('artifacts/reproduction_inputs/v5a_round29')


def load(path):
    return rt.json.loads(Path(path).read_text())


def close(a, b):
    assert math.isclose(float(a), float(b), abs_tol=1e-13, rel_tol=1e-12), (a, b)


def main():
    assert not OUT.exists(), 'Retain the original round receipt; use a new name for a replay'
    tracked = ['artifacts/calibration_gate.json', 'artifacts/calibration.json',
               'artifacts/tasks/manifest.json', 'artifacts/benchmark_roles.json',
               'artifacts/anchor_policy.json', 'artifacts/phase1_calibration_audit.json']
    before = {p: rt.sha(p) for p in tracked}
    assert load(tracked[0])['status'] == 'CLOSED'
    assert load('artifacts/benchmark_roles.json')['primary_search'] == TASKS
    index = load('artifacts/frozen_metadata_index.json')
    for path, entry in index['entries'].items():
        assert rt.sha(path) == entry['sha256'], path
    assert rt.sha(index['archive']) == index['archive_sha256']
    environment = scoring.provenance()
    refs, anchors, historical = {}, {}, {}
    for task in TASKS:
        print('Checking current SEARCH identities', task, flush=True)
        contract, ctrl, real = benchmark.setup(task)
        assert contract['split'] == 'search'
        refs[task] = scoring.matrix_id(real)
        anchors[task] = load(f'artifacts/scores/{task}/anchors.json')
        assert anchors[task]['real_sha256'] == refs[task]
        assert anchors[task]['n_splits'] == 5
        assert anchors[task]['driver_sha256'] == rt.sha('bounded_score.py')
        assert anchors[task]['scorer_sha256'] == rt.sha('official_scoring.py')
        current = provenance.identity(task)
        for recipe in RECIPES:
            path = f'artifacts/scores/{task}/{recipe}_s42.json'
            record = load(path)
            assert record['implementation_sha256'] == current
            assert record['reference_sha256'] == refs[task]
            exact_scaled(record, anchors[task])
            historical[path] = rt.sha(path)
        del ctrl, real
        gc.collect()

    baselines = {}
    for stem in ['atlas_loco', 'sidechain_loco']:
        spec_path = f'artifacts/{stem}_spec.json'
        audit_path = f'artifacts/{stem}_audit.json'
        results_path = f'artifacts/{stem}_results.json'
        spec, audit, results = map(load, [spec_path, audit_path, results_path])
        assert results['status'].startswith('PASS completed')
        assert results['spec_sha256'] == rt.sha(spec_path)
        assert audit['status'] == 'PASS' and audit['spec_sha256'] == rt.sha(spec_path)
        assert audit['verifier_sha256'] == rt.sha(f'verify_{stem}.py')
        assert results['runner_sha256'] == rt.sha(f'run_{stem}.py')
        assert spec['official_environment'] == environment
        checked = {}
        for group in ['implementation_sha256', 'protocol_sha256']:
            for path, digest in spec[group].items():
                assert rt.sha(path) == digest, path
            checked[group] = len(spec[group])
        # The two already-frozen specs use different source receipt schemas.
        sources = (spec['source_sha256'] if stem == 'sidechain_loco' else
                   {p: rec['sha256'] for p, rec in spec['source_files'].items()})
        for path, digest in sources.items():
            assert rt.sha(path) == digest, path
        checked['source_sha256'] = len(sources)
        assert results['gate_unchanged_sha256'] == before['artifacts/calibration_gate.json']
        assert set(results['tasks']) == set(TASKS)
        task_scores = {}
        for task, record in results['tasks'].items():
            assert record['spec_sha256'] == rt.sha(spec_path)
            assert record['audit_sha256'] == rt.sha(audit_path)
            assert record['reference_sha256'] == refs[task]
            assert record['anchors_sha256'] == rt.sha(f'artifacts/scores/{task}/anchors.json')
            saved = list((Path('artifacts/baseline_scores') / results['baseline_id'] / task).glob('s*.json'))
            assert len(saved) == 1 and load(saved[0]) == record
            exact_scaled(record, anchors[task])
            diagnostic = record.get('diagnostics', {}).get('prediction_only_B_selfclaim')
            if diagnostic is not None:
                close(diagnostic['B'] / diagnostic['full_prediction_jackknife_claim'], diagnostic['ratio'])
            task_scores[task] = dict(raw=record['raw'], scaled=record['scaled'],
                overall=record['overall'], prediction_sha256=record['prediction_sha256'],
                prediction_only_B_selfclaim=diagnostic)
        macro = sum(Fraction(r['overall']) for r in results['tasks'].values()) / len(TASKS)
        close(macro, results['equal_context_macro'])
        assert results['external_official_label'] is None
        assert results['predicted_official_score'] is None
        assert results['holdout_scored'] is False and results['phase2_allowed'] is False
        baselines[results['baseline_id']] = dict(results_sha256=rt.sha(results_path),
            checked_bindings=checked, tasks=task_scores, equal_context_macro=float(macro),
            paired_external_score=None)

    # Freeze the actual v5a build inputs, not the mutable narrative or local ablations.
    parent = Path('../model_v5')
    build = load(parent / 'artifacts/v5a_build.json')
    package = load(parent / 'artifacts/v5a_package.json')
    for key, name in [('v5_sha256', 'v5.py'), ('model_sha256', 'model.py'),
                      ('builder_sha256', 'build_v5.py')]:
        assert rt.sha(parent / name) == build[key], (name, 'Build input changed')
    package_path = parent / 'v5a.vcc'
    assert package_path.stat().st_size == package['package_bytes']
    print('Checking v5a package identity (read only)', flush=True)
    assert rt.sha(package_path) == package['package_sha256']
    assert package['package_metadata']['n_obs'] == build['shape'][0] == 360000
    assert package['package_metadata']['n_vars'] == build['shape'][1] == 18533
    assert package['package_metadata']['nnz'] == build['nnz']
    assert build['arm'] == 'a'
    for context in build['contexts']:
        assert context['rho'] == 1.0
        close(context['B_emitted'] / context['claim_pred_emitted'], context['ratio_B_over_claim_pred'])
    names = ['v5.py', 'model.py', 'build_v5.py', 'runtime.py', 'verify_counts.py',
             'vendor/atlas_model.py', 'vendor/atlas_predict.py', 'vendor/adaptive_shrink.py',
             'vendor/ATLAS_LICENSE', 'vendor/SIDECHAIN_LICENSE',
             'artifacts/v5a_build.json', 'artifacts/v5a_package.json',
             'artifacts/v5a_structure.json', 'NOTES.md']
    snapshot = {}
    for name in names:
        source, destination = parent / name, FROZEN / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            assert rt.sha(destination) == rt.sha(source), destination
        else:
            shutil.copyfile(source, destination)
        snapshot[str(destination)] = rt.sha(destination)
    for source, name in [('../OFFICIAL_SCORES.md', 'supplied_official_scores.md'),
                         ('artifacts/prior/scorer_metrics', 'cached_arc_metrics.md'),
                         ('../data/manifest.json', 'validation_manifest.json')]:
        destination = FROZEN / name
        if destination.exists():
            assert rt.sha(destination) == rt.sha(source)
        else:
            shutil.copyfile(source, destination)
        snapshot[str(destination)] = rt.sha(destination)

    # An input-control count cannot establish equality of cell identities.
    doc = (FROZEN / 'cached_arc_metrics.md').read_text()
    assert 'held-out real control cells out of the reference' in doc
    assert 'control cells supplied with the phase are model **inputs** only' in doc
    assert load(FROZEN / 'validation_manifest.json')['per_context']['A']['control_cells'] == 18400
    official_line, = [line for line in (FROZEN / 'supplied_official_scores.md').read_text().splitlines()
                      if line.startswith('| v5a ')]
    assert float(official_line.split('|')[2].strip().replace('**', '')) == .1499

    manifest = load('artifacts/tasks/manifest.json')
    previous = load('artifacts/protocol_revisions/CD4_depth192/manifest.json.before')
    for task, contract in manifest.items():
        if contract['split'] == 'holdout':
            assert contract == previous[task]
    assert not Path('artifacts/holdout_events.jsonl').exists()
    assert not [p for p in Path('artifacts/scores').glob('*/*_s42.json') if 'holdout' in str(p).lower()]
    for path, digest in before.items():
        assert rt.sha(path) == digest
    for path, digest in historical.items():
        assert rt.sha(path) == digest

    decision = dict(
        status='RETAIN CLOSED: existing benchmark is not validated for model selection',
        task_or_scale_change=False, historical_five_recipe_gate_unchanged=True,
        reason='Historical v4a/v4b reversal remains. Scoped public-kernel baselines are unpaired with external entries; they cannot repair that gate.',
        controls='Retain disjoint input/reference controls. Cached Arc documentation describes held-out reference controls; equal counts in the public manifest do not certify identical cell identities.',
        sidechain_scope='Completed available-moment kernels only; missing SER17 submitted source/configuration/floor/neighbours/scatter prevent exact SER17 reproduction.',
        budget_interpretation='H1 B/full-self-claim >1 is sufficient for the full capped correction but does not guarantee low MSE. CD4 ratio <1 does not establish the actual retained correction because the reference-side min cap is unmeasured here.',
        next_unit='Implement and independently audit one frozen v5a available-input LOCO port using the archived budget/solve/emitter kernels and only input controls. Score both unchanged SEARCH contexts only after audit PASS. This is a Phase1 reproduction extension, never a replacement for or escape from the failed historical rank gate.',
        next_port_contract=dict(
            parent='Existing frozen local v4a port; both SEARCH tasks unchanged',
            rule='v5a build arm a: rho=1.0; 40-target prediction-only jackknife probe; solve_alpha with control-bulk centre; pooled four-cell template; original integer repair',
            exclusions='Inherit current destination-family and globally sealed HEK293T/RPE1 exclusions, including all H1 response branches excluded from v4.',
            input_interface='No real/reference argument to prediction construction or alpha selection. Reference controls may enter official scoring only.',
            seeds='Freeze explicitly before scores: template20260910 and seed_for(destination:target,20260910), following existing local v4a convention; probe RNG0 for H1 and RNG1 for CD4, fixed40 of79.',
            failure_policy='Abort missing support, infeasible fits or identity mismatch; no alternate amplitude, probe size, centring, context or fallback chosen from scores.',
            verification='Audit original helper equality, independent B/row-weight/claim arithmetic, v4a head/emission bridge, reference isolation and unchanged historical identities before scoring.',
            scope='LOCO adaptation on restricted measured axes and source banks, not byte equality with hidden-context submission. Official0.1499 is supplied/unverified server evidence, not a forecast or gate datum.',
            phase2_allowed=False, holdout_scoring_allowed=False),
        exact_external_reproduction_obstacles_retained=True,
        public_sources_refreshed=False)
    rt.save(OUT, dict(utc=datetime.now(timezone.utc).isoformat(), status='PASS',
        auditor_sha256=rt.sha(__file__), frozen_metadata_members=len(index['entries']),
        baselines=baselines, historical_identities=historical,
        preserved_protocol_sha256=before, v5a_input_snapshot_sha256=snapshot,
        v5a_package_sha256=package['package_sha256'], v5a_code_build_identity='PASS',
        v5a_package_identity='PASS SHA256/size/build metadata; not a fresh CSR decode or server attestation',
        benchmark_validity_decision=decision, holdout_scored=False, new_public_read_bytes=0))
    print('PASS baseline arithmetic/current identities; v5a evidence frozen; gate CLOSED', flush=True)


if __name__ == '__main__':
    main()
